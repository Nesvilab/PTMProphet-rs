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

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use mzdata::io::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::ArrayType;

use crate::error::{PtmProphetRsError, Result};
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
    /// Peak m/z values (mzML stores them ascending).
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
        let scan_to_index = build_scan_index(&reader);
        if scan_to_index.is_empty() {
            return Err(PtmProphetRsError::Mzml(format!(
                "empty or unindexed spectrum list in {}",
                path.display()
            )));
        }
        Ok(MzmlReader { path: path.to_path_buf(), scan_to_index, inner: Mutex::new(reader) })
    }

    /// Number of indexed spectra whose native-id carries a scan number.
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

/// Build `scan → ordinal` from the reader's spectrum offset index (file order).
fn build_scan_index(reader: &Reader) -> HashMap<i64, usize> {
    let mut map = HashMap::with_capacity(reader.spectrum_index.len());
    for (ord, id) in reader.spectrum_index.keys().enumerate() {
        if let Some(scan) = scan_from_idref(id) {
            map.insert(scan, ord);
        }
    }
    map
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

    let (precursor_mz, precursor_charge) = match spec.precursor().and_then(|p| p.ion()) {
        Some(ion) => (ion.mz, ion.charge.unwrap_or(0)),
        None => (0.0, 0),
    };

    Spectrum { scan, ms_level: spec.ms_level(), precursor_mz, precursor_charge, mz, intensity }
}

/// Extract the scan number from a spectrum nativeID like `... scan=3185`.
fn scan_from_idref(idref: &str) -> Option<i64> {
    if let Some(pos) = idref.find("scan=") {
        let rest = &idref[pos + 5..];
        let n: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
        return n.parse().ok();
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_idref_parse() {
        assert_eq!(scan_from_idref("controllerType=0 controllerNumber=1 scan=3185"), Some(3185));
        assert_eq!(scan_from_idref("nonsense"), None);
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
