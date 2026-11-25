// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Indexed mzML reader, backed by [`mzdata`].
//!
//! `mzdata` parses the `indexedmzML` `<indexList>` into an offset index (native-id → byte offset)
//! and decodes a `<spectrum>` on demand (base64 → optional zlib inflate → 32/64-bit floats). We map
//! its `MultiLayerSpectrum` onto the small [`Spectrum`] view that downstream scoring consumes: peak
//! list (m/z + intensity), MS level, precursor m/z + charge, scan number.
//!
//! Precision note: intensities are kept as **f64**. `BinaryArrayMap::intensities()` would downcast a
//! 64-bit intensity array to `f32`; instead we read the intensity `DataArray` and call `to_f64()`,
//! which is lossless for both 32- and 64-bit source arrays, so the O-score stays bit-identical.
//!
//! Throughput: [`MzmlReader::load_scans`] decodes a batch of scans across worker threads (each owns
//! its own `mzdata` reader on the same file, covering a contiguous slice of the offset index for
//! forward-only seeks). The single-scan [`MzmlReader::spectrum_for_scan`] is kept for tests and
//! one-off lookups.
//!
//! Spectrum lookup: a PSM's pepXML `start_scan` is matched against the scan number MSFragger derives
//! from each native ID ([`scan_number`]), so every native-ID format MSFragger reads resolves,
//! including diaTracer's `frame=N`.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use mzdata::io::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::ArrayType;

use crate::error::{PtmProphetRsError, Result};
use crate::spectrum::scan_id::scan_number;
use crate::threads::worker_count;

/// `mzdata` reader specialized to the default peak types, reading from a `File`.
type Reader = MzMLReader<File>;

/// One MS/MS spectrum, with only the fields downstream scoring consumes.
#[derive(Debug, Clone, Default)]
pub struct Spectrum {
    pub scan: i64,
    pub ms_level: u8,
    pub precursor_mz: f64,
    pub precursor_charge: i32,
    /// Peak m/z values, ascending: the reader sorts peaks that the file stores in another order.
    pub mz: Vec<f64>,
    /// Peak intensities, parallel to `mz`.
    pub intensity: Vec<f64>,
}

impl Spectrum {
    /// Smallest positive intensity among peaks (1.0 if no peak is positive) — the O-score
    /// denominator, PTMProphet's `minInt_`. Matching never lowers it: every matched peak is a peak
    /// of this spectrum, so none is weaker (see [`crate::mpx::Mpx::process`]).
    pub fn min_intensity(&self) -> f64 {
        let m = self
            .intensity
            .iter()
            .copied()
            .filter(|&v| v > 0.0)
            .fold(f64::INFINITY, f64::min);
        if m.is_finite() && m > 0.0 { m } else { 1.0 }
    }
    pub fn max_intensity(&self) -> f64 {
        self.intensity.iter().copied().fold(0.0, f64::max)
    }
    pub fn tic(&self) -> f64 {
        self.intensity.iter().sum()
    }
}

/// An opened, scan-indexed mzML file, reused across a run's spectrum queries.
pub struct MzmlReader {
    path: PathBuf,
    /// scan number → ordinal position in the file's spectrum offset index.
    scan_to_index: HashMap<i64, usize>,
    /// A reader kept open for single-scan random access (tests / one-off lookups). The batch
    /// [`load_scans`](Self::load_scans) path opens its own per-worker readers instead.
    inner: Mutex<Reader>,
}

impl MzmlReader {
    /// Open an `indexedmzML` and build its scan→ordinal map.
    pub fn open(path: &Path) -> Result<MzmlReader> {
        let reader = Reader::open_path(path.to_path_buf())
            .map_err(|e| PtmProphetRsError::Io { path: path.to_path_buf(), source: e })?;
        let (scan_to_index, repeated) = build_scan_index(&reader);
        if scan_to_index.is_empty() {
            return Err(PtmProphetRsError::Mzml(format!(
                "empty or unindexed spectrum list in {}",
                path.display()
            )));
        }
        if repeated > 0 {
            log::warn!(
                "{} spectra in {} repeat an earlier spectrum's scan number; a PSM looked up by such a scan number gets the later spectrum",
                repeated,
                path.display()
            );
        }
        Ok(MzmlReader { path: path.to_path_buf(), scan_to_index, inner: Mutex::new(reader) })
    }

    /// Number of distinct scan numbers among the indexed spectra.
    pub fn len(&self) -> usize {
        self.scan_to_index.len()
    }
    pub fn is_empty(&self) -> bool {
        self.scan_to_index.is_empty()
    }

    /// Fetch the spectrum for `scan`. Errors if the scan is not indexed or unreadable.
    pub fn spectrum_for_scan(&self, scan: i64) -> Result<Spectrum> {
        let ord = self.ordinal_for_scan(scan)?;
        let mut guard = self
            .inner
            .lock()
            .map_err(|_| PtmProphetRsError::Mzml("mzML reader lock poisoned".into()))?;
        let spec = guard.get_spectrum_by_index(ord).ok_or_else(|| {
            PtmProphetRsError::Mzml(format!("scan {scan} (index {ord}) unreadable in {}", self.path.display()))
        })?;
        Ok(convert_spectrum(&spec, scan))
    }

    fn ordinal_for_scan(&self, scan: i64) -> Result<usize> {
        self.scan_to_index.get(&scan).copied().ok_or_else(|| {
            PtmProphetRsError::Mzml(format!("scan {scan} not found in {}", self.path.display()))
        })
    }

    /// Decode a batch of scans in parallel, returning `scan → Spectrum` for those present.
    ///
    /// Scans not in the index are silently omitted (the caller tolerates skipping missing scans).
    /// Work is split into `worker_count` contiguous slices of the offset index so each worker seeks
    /// forward through its region; results are merged into one map.
    pub fn load_scans(&self, scans: &[i64]) -> Result<HashMap<i64, Spectrum>> {
        // Resolve to (ordinal, scan), drop unknown scans, sort by ordinal for forward-only seeks,
        // and dedup (a scan resolves to exactly one ordinal).
        let mut work: Vec<(usize, i64)> = scans
            .iter()
            .filter_map(|&s| self.scan_to_index.get(&s).map(|&ord| (ord, s)))
            .collect();
        work.sort_unstable();
        work.dedup();
        if work.is_empty() {
            return Ok(HashMap::new());
        }

        let n_workers = worker_count().min(work.len()).max(1);
        let chunk_len = work.len().div_ceil(n_workers);
        let path = self.path.as_path();

        let partials: Vec<Result<HashMap<i64, Spectrum>>> = std::thread::scope(|scope| {
            let handles: Vec<_> = work
                .chunks(chunk_len)
                .map(|chunk| scope.spawn(move || decode_chunk(path, chunk)))
                .collect();
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .unwrap_or_else(|_| Err(PtmProphetRsError::Mzml("spectrum worker panicked".into())))
                })
                .collect()
        });

        let mut out = HashMap::with_capacity(work.len());
        for partial in partials {
            out.extend(partial?);
        }
        Ok(out)
    }
}

/// Decode one contiguous slice of `(ordinal, scan)` work using a dedicated reader.
fn decode_chunk(path: &Path, chunk: &[(usize, i64)]) -> Result<HashMap<i64, Spectrum>> {
    let mut reader =
        Reader::open_path(path.to_path_buf()).map_err(|e| PtmProphetRsError::Io { path: path.to_path_buf(), source: e })?;
    let mut local = HashMap::with_capacity(chunk.len());
    for &(ord, scan) in chunk {
        if let Some(spec) = reader.get_spectrum_by_index(ord) {
            local.insert(scan, convert_spectrum(&spec, scan));
        }
    }
    Ok(local)
}

/// Build `scan → ordinal` from the reader's spectrum offset index (file order), numbering each
/// spectrum as MSFragger does, with its position plus 1 as the fallback. Also returns how many
/// spectra repeat an earlier spectrum's scan number; the later spectrum keeps the number.
fn build_scan_index(reader: &Reader) -> (HashMap<i64, usize>, usize) {
    let mut map = HashMap::with_capacity(reader.spectrum_index.len());
    let mut repeated = 0;
    for (ord, id) in reader.spectrum_index.keys().enumerate() {
        if map.insert(scan_number(id, ord as i64 + 1), ord).is_some() {
            repeated += 1;
        }
    }
    (map, repeated)
}

/// Map an `mzdata` spectrum onto the [`Spectrum`] view. m/z and intensity are read as f64.
fn convert_spectrum<S>(spec: &S, scan: i64) -> Spectrum
where
    S: SpectrumLike,
{
    let (mut mz, mut intensity) = match spec.raw_arrays() {
        Some(arrays) => {
            let mz = arrays.mzs().map(|c| c.into_owned()).unwrap_or_default();
            let intensity = arrays
                .get(&ArrayType::IntensityArray)
                .and_then(|da| da.to_f64().ok())
                .map(|c| c.into_owned())
                .unwrap_or_default();
            (mz, intensity)
        }
        None => (Vec::new(), Vec::new()),
    };
    // Defensive: a corrupt file could yield mismatched array lengths; downstream code indexes the
    // two arrays in lockstep, so truncate both to the shared length.
    if mz.len() != intensity.len() {
        let n = mz.len().min(intensity.len());
        mz.truncate(n);
        intensity.truncate(n);
    }
    // Peak matching binary-searches m/z, so the peaks must be in m/z order. diaTracer writes them by
    // descending intensity instead; PTMProphet's SpectraST peak list sorts by m/z before any lookup.
    if !mz.is_sorted() {
        let mut order: Vec<usize> = (0..mz.len()).collect();
        order.sort_by(|&a, &b| mz[a].total_cmp(&mz[b]));
        mz = order.iter().map(|&i| mz[i]).collect();
        intensity = order.iter().map(|&i| intensity[i]).collect();
    }

    let (precursor_mz, precursor_charge) = match spec.precursor().and_then(|p| p.ion()) {
        Some(ion) => (ion.mz, ion.charge.unwrap_or(0)),
        None => (0.0, 0),
    };

    Spectrum { scan, ms_level: spec.ms_level(), precursor_mz, precursor_charge, mz, intensity }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Standard padded base64, enough to embed peak arrays in a test mzML.
    fn base64(bytes: &[u8]) -> String {
        const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::new();
        for chunk in bytes.chunks(3) {
            let byte = |i: usize| chunk.get(i).copied().unwrap_or(0) as u32;
            let n = byte(0) << 16 | byte(1) << 8 | byte(2);
            for k in 0..4 {
                if k <= chunk.len() {
                    out.push(TABLE[(n >> (18 - 6 * k) & 63) as usize] as char);
                } else {
                    out.push('=');
                }
            }
        }
        out
    }

    /// Write an indexed mzML with one MS2 spectrum per native ID. With `base = 100 * (i + 1)` for the
    /// spectrum at 0-based position `i`, its peaks are written by descending intensity, as diaTracer
    /// writes them: m/z `base + 50, base, base + 20` with intensities `30, 20, 10`.
    fn write_indexed_mzml(path: &Path, ids: &[&str]) {
        let mut xml = String::from(concat!(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n",
            "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">\n",
            "  <mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n",
            "    <cvList count=\"1\">\n",
            "      <cv id=\"MS\" fullName=\"PSI-MS\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n",
            "    </cvList>\n",
            "    <fileDescription>\n",
            "      <fileContent>\n",
            "        <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n",
            "      </fileContent>\n",
            "    </fileDescription>\n",
            "    <softwareList count=\"1\">\n",
            "      <software id=\"test\" version=\"0\"/>\n",
            "    </softwareList>\n",
            "    <instrumentConfigurationList count=\"1\">\n",
            "      <instrumentConfiguration id=\"IC1\"/>\n",
            "    </instrumentConfigurationList>\n",
            "    <dataProcessingList count=\"1\">\n",
            "      <dataProcessing id=\"DP1\">\n",
            "        <processingMethod order=\"0\" softwareRef=\"test\"/>\n",
            "      </dataProcessing>\n",
            "    </dataProcessingList>\n",
            "    <run id=\"test\" defaultInstrumentConfigurationRef=\"IC1\">\n",
        ));
        xml.push_str(&format!(
            "      <spectrumList count=\"{}\" defaultDataProcessingRef=\"DP1\">\n",
            ids.len()
        ));
        let mut offsets = Vec::with_capacity(ids.len());
        for (i, id) in ids.iter().enumerate() {
            xml.push_str("        ");
            offsets.push((*id, xml.len()));
            let base = 100.0 * (i + 1) as f64;
            let encode = |values: [f64; 3]| base64(&values.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<u8>>());
            let mz = encode([base + 50.0, base, base + 20.0]);
            let intensity = encode([30.0, 20.0, 10.0]);
            xml.push_str(&format!(
                concat!(
                    "<spectrum index=\"{i}\" id=\"{id}\" defaultArrayLength=\"3\">\n",
                    "          <cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>\n",
                    "          <cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>\n",
                    "          <binaryDataArrayList count=\"2\">\n",
                    "            <binaryDataArray encodedLength=\"{mz_len}\">\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\" value=\"\"/>\n",
                    "              <binary>{mz}</binary>\n",
                    "            </binaryDataArray>\n",
                    "            <binaryDataArray encodedLength=\"{int_len}\">\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\" value=\"\"/>\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" value=\"\"/>\n",
                    "              <cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" value=\"\"/>\n",
                    "              <binary>{int}</binary>\n",
                    "            </binaryDataArray>\n",
                    "          </binaryDataArrayList>\n",
                    "        </spectrum>\n",
                ),
                i = i,
                id = id,
                mz_len = mz.len(),
                mz = mz,
                int_len = intensity.len(),
                int = intensity,
            ));
        }
        xml.push_str("      </spectrumList>\n    </run>\n  </mzML>\n");
        let index_offset = xml.len() + 2;
        xml.push_str("  <indexList count=\"1\">\n    <index name=\"spectrum\">\n");
        for (id, offset) in offsets {
            xml.push_str(&format!("      <offset idRef=\"{id}\">{offset}</offset>\n"));
        }
        xml.push_str("    </index>\n  </indexList>\n");
        xml.push_str(&format!("  <indexListOffset>{index_offset}</indexListOffset>\n</indexedmzML>\n"));
        std::fs::write(path, xml).expect("write test mzML");
    }

    #[test]
    fn frame_ids_resolve_to_their_scan_numbers() {
        // diaTracer and MSFragger's calibrated mzML name spectra `frame=N`; MSFragger's pepXML
        // start_scan is N. A gap in the frame numbers shows N is parsed, not taken from the position.
        let path = std::env::temp_dir()
            .join(format!("ptmprophet-rs-frame-ids-{}.mzML", std::process::id()));
        write_indexed_mzml(&path, &["frame=1", "frame=3", "frame=4"]);

        let reader = MzmlReader::open(&path).expect("open mzML with frame=N ids");
        assert_eq!(reader.len(), 3);
        let s = reader.spectrum_for_scan(3).expect("scan 3");
        assert_eq!(s.mz, vec![200.0, 220.0, 250.0], "scan 3 is the second spectrum, peaks in m/z order");
        assert_eq!(s.intensity, vec![20.0, 10.0, 30.0], "intensities move with their m/z");
        assert_eq!(s.ms_level, 2);
        assert!(reader.spectrum_for_scan(2).is_err(), "no spectrum has scan number 2");

        let batch = reader.load_scans(&[4, 1, 2]).expect("batch load");
        assert_eq!(batch.len(), 2, "scan 2 is absent and omitted");
        assert_eq!(batch[&1].mz, vec![100.0, 120.0, 150.0]);
        assert_eq!(batch[&4].mz, vec![300.0, 320.0, 350.0]);
        assert_eq!(batch[&4].intensity, vec![20.0, 10.0, 30.0]);

        drop(reader);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn read_real_diatracer_mzml() {
        let path = Path::new(
            "F:/data/diaPASEF/20200505_Evosep_200SPD_SG06-16_MLHeLa_200ng_py8_S3-A4_1_2740_diatracer.mzML",
        );
        if !path.exists() {
            eprintln!("[skip] mzML not present: {}", path.display());
            return;
        }
        let reader = MzmlReader::open(path).expect("open diaTracer mzML");
        assert_eq!(reader.len(), 114_998, "one scan number per frame=N spectrum");
        // MSFragger's pepXML for this file has start_scan="10090" spectrumNativeID="frame=10090".
        let s = reader.spectrum_for_scan(10090).expect("scan 10090");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.mz.len(), 500);
        assert_eq!(s.mz.len(), s.intensity.len());
        // diaTracer stores the peaks by descending intensity; the reader returns them by m/z.
        assert!(s.mz.windows(2).all(|w| w[0] <= w[1]), "m/z sorted");
    }

    #[test]
    fn read_real_mzml_scan_3185() {
        let path = Path::new("F:/data/ccrcc_dda_phospho/mzml/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05_calibrated.mzML");
        if !path.exists() {
            eprintln!("[skip] mzML not present: {}", path.display());
            return;
        }
        let reader = MzmlReader::open(path).expect("open mzML");
        assert!(reader.len() > 1000, "indexed spectra = {}", reader.len());
        let s = reader.spectrum_for_scan(3185).expect("scan 3185");
        assert_eq!(s.ms_level, 2, "ms level");
        assert!(!s.mz.is_empty() && s.mz.len() == s.intensity.len(), "peaks {} / {}", s.mz.len(), s.intensity.len());
        // m/z ascending
        assert!(s.mz.windows(2).all(|w| w[0] <= w[1]), "m/z sorted");
        assert!(s.min_intensity() > 0.0);

        // The batch loader returns the same peaks as the single-scan path.
        let batch = reader.load_scans(&[3185]).expect("batch load");
        let b = batch.get(&3185).expect("scan 3185 in batch");
        assert_eq!(b.mz, s.mz);
        assert_eq!(b.intensity, s.intensity);

        eprintln!(
            "[ok] scan 3185: msLevel={} peaks={} precMz={:.4} z={} minInt={:.1} maxInt={:.1}",
            s.ms_level, s.mz.len(), s.precursor_mz, s.precursor_charge, s.min_intensity(), s.max_intensity()
        );
    }
}
