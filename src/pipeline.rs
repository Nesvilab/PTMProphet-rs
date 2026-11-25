// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! Top-level orchestration for the file-list's input files. Each PSM is independent, so spectrum
//! decoding and per-PSM scoring fan out across all logical CPUs (see [`crate::spectrum::mzml`]).
//! Results are merged in **input order**, so each output is byte-identical regardless of worker count.
//!
//! A single invocation localizes every pepXML named by the `--filelist` (one or many), but each file
//! is localized with **its own** EM intensity model trained on **that file's** ambiguous-site
//! evidence — not one pooled model. PTMProphet is invoked once per pepXML, so its `<mixturemodel>`
//! (bandwidth + density curve) differs per file; training per file is what reproduces its output.
//!
//! Passes:
//!   A. [`scan_all_inputs`] — scan every pepXML once (read + parse + PSM collection run **in parallel
//!      across files**), then merge into one combined index (mod tables + probability hash) serially
//!      in file order.
//!   B/C. Per file, in input order: build its [`crate::mpx::Mpx`]s (parallel across the file's PSMs)
//!      and match fragments once (spectra decoded then dropped); harvest *that file's* ambiguous-site
//!      O-scores and fit its own KDModel via [`crate::em::train::compute_ptm_model`]; localize each
//!      top hit with that model and build the per-PSM output blocks (parallel across the file's PSMs);
//!      then stream every file's output pepXML **in parallel across files** (independent outputs).
//!
//! Parallelism is layered to keep all cores busy without oversubscribing: the per-file work that has
//! no internal parallelism (PASS A's parse, the output write) fans out *across files*, while the
//! heavy per-PSM work (matching, block building) fans out *within* a file processed one at a time —
//! the two layers never nest. Everything reassembles in file/input order, so each output is identical
//! regardless of worker count. Each input pepXML is read at most twice (index/PSM scan, then the
//! streaming copy) and every needed spectrum is decoded once; fragment-matching runs once per PSM.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::config::{Config, InputSpec};
use crate::em::kdmodel::KdModel;
use crate::em::train::EmEvidence;
use crate::error::{PtmProphetRsError, Result};
use crate::mpx::Mpx;
use crate::pepxml::index::InputIndex;
use crate::pepxml::psm::{for_each_psm, PsmRecord};
use crate::pepxml::roc::InfoLocs as RocInfoLocs;
use crate::spectrum::mzml::{MzmlReader, Spectrum};
use crate::threads::{par_for_each_mut, par_map, par_map_mut};

/// One eligible PSM and its localization engine. `processed` is set once a spectrum was found and
/// [`Mpx::process`] has run; un-processed items (whose spectrum is missing) contribute to neither
/// the EM nor the output.
struct WorkItem<'a> {
    psm: PsmRecord,
    mpx: Mpx<'a>,
    processed: bool,
}

/// One input pepXML's processed work, carried from the match pass through to its own output write.
/// `items` borrow `&config.mods` (not any per-file data), so every file's items share one lifetime
/// and can be held together across the per-file EM train + build loop.
struct FileWork<'a> {
    spec: &'a InputSpec,
    items: Vec<WorkItem<'a>>,
}

/// One file's finished output payload: the rewritten per-PSM blocks + ROC summary + this file's own
/// trained EM model, ready to stream. Owns its paths (no borrow of the `WorkItem`s), so the per-file
/// match state can be freed before the parallel write fan-out.
struct FileOutput {
    in_path: PathBuf,
    out_path: PathBuf,
    roc_block: String,
    /// This file's own trained intensity model (its `<mixturemodel>` goes into this file's header).
    model: KdModel,
    /// Per-`start_scan` rewritten `<modification_info>` + `<analysis_result>` block.
    blocks: HashMap<i64, String>,
}

/// Decoded spectra keyed by run-file path, then by scan number.
type SpectrumStore = HashMap<String, HashMap<i64, Spectrum>>;

/// Run the full localization pipeline over every input in [`Config::inputs`] (each file gets its own
/// per-file EM model).
pub fn run(config: &Config) -> Result<()> {
    let mut timer = PhaseTimer::new();
    let n_files = config.inputs.len();

    // PASS A — scan every input once: one combined index (merged mod tables + probability hash) and
    // each file's probability-gated PSM list, in file order.
    let (index, per_file_psms) = scan_all_inputs(config)?;
    timer.mark_n("scan_all_inputs (PASS A)", n_files);

    // Per file: build its work items (all borrow `&config.mods`) and match fragments once. Each
    // file's spectra are decoded then dropped immediately — only the per-site scores left on each
    // `Mpx` persist into that file's own EM train and output write.
    let mut files: Vec<FileWork> = Vec::with_capacity(n_files);
    let mut total_items = 0usize;
    for (i, (spec, psms)) in config.inputs.iter().zip(per_file_psms).enumerate() {
        let mut items = build_work_items(config, &index, psms);
        let store = load_store(&items)?;
        process_all(&mut items, &store);
        drop(store);
        total_items += items.len();
        if n_files > 1 {
            log::info!(
                "File {}/{}: {} — {} PSMs to localize",
                i + 1,
                n_files,
                display_name(&spec.pep_xml),
                items.len()
            );
        }
        files.push(FileWork { spec, items });
    }
    timer.mark_n("process_all (parallel match, all files)", total_items);

    // PASS B + C — per file: fit that file's OWN intensity model on its ambiguous-site evidence, then
    // localize with it and build the file's output blocks. PTMProphet is invoked once per pepXML, so
    // it trains a separate model per file (its `<mixturemodel>` bandwidth/curve differs per file);
    // pooling across files would produce one shared model and shift every probability. Each file's
    // EM harvest is serial in input order (independent of worker count); the heavy per-PSM build fans
    // out within the file. Splitting build from write lets the writes (read-input + serial XML copy —
    // otherwise a wholly-serial tail) overlap across files below, without nesting two parallel layers.
    let ptm_labels: Vec<String> = config.mods.iter().map(|m| m.ptm_label.clone()).collect();
    let mut outputs: Vec<FileOutput> = Vec::with_capacity(n_files);
    for fw in &mut files {
        // Fit this file's own model on its own ambiguous-site evidence (in input order).
        let mut model = KdModel::new("PTMProphet_Oscore");
        if config.em != 0 {
            let mut evidence = EmEvidence::default();
            harvest_into(&fw.items, &mut evidence);
            crate::em::train::compute_ptm_model(config, &mut model, &evidence)?;
        }

        let contributions = build_contributions(&mut fw.items, config, &index, &model);

        // Assemble blocks (last-wins on a duplicate scan, in input order) and the ROC accumulator
        // (prob-weighted InfoLocs per mod type, pushed in input order).
        let mut blocks: HashMap<i64, String> = HashMap::with_capacity(contributions.len());
        let mut roc_by_type: Vec<Vec<RocInfoLocs>> = vec![Vec::new(); config.mods.len()];
        for (scan, block, rocs) in contributions.into_iter().flatten() {
            for (t, r) in rocs {
                roc_by_type[t].push(r);
            }
            blocks.insert(scan, block);
        }
        let roc_block = crate::pepxml::roc::build_roc_block(&roc_by_type, &ptm_labels);
        outputs.push(FileOutput {
            in_path: fw.spec.pep_xml.clone(),
            out_path: fw.spec.out.clone(),
            roc_block,
            model,
            blocks,
        });
    }
    timer.mark_n("per-file EM train + build (PASS B+C)", n_files);

    // Free the per-file match state before the write fan-out (the blocks are self-contained Strings).
    drop(files);

    // Write every output in parallel — each re-reads its input pepXML and streams the localized copy
    // to its own path. Errors propagate in file order.
    let stat_mods = &index.mods.stat_mods;
    let written: Vec<Result<()>> = par_map(&outputs, |o| {
        let content = read_file(&o.in_path)?;
        crate::pepxml::writer::stream_output(
            config,
            stat_mods,
            &o.in_path,
            &o.out_path,
            &content,
            &o.model,
            &o.roc_block,
            &o.blocks,
        )
    });
    for result in written {
        result?;
    }
    timer.mark_n("stream_output (PASS C write, parallel)", n_files);
    timer.total();
    Ok(())
}

/// Lightweight per-phase wall-clock timer. Logs (at `info`) only when `PTMPROPHET_RS_TIMING` is set
/// (so it is free in normal runs). Helps confirm where time goes and that the parallel phases
/// scale.
struct PhaseTimer {
    on: bool,
    start: std::time::Instant,
    last: std::time::Instant,
}
impl PhaseTimer {
    fn new() -> Self {
        let on = std::env::var("PTMPROPHET_RS_TIMING").map(|v| !v.is_empty() && v != "0").unwrap_or(false);
        let now = std::time::Instant::now();
        PhaseTimer { on, start: now, last: now }
    }
    fn mark_n(&mut self, label: &str, n: usize) {
        if self.on {
            let now = std::time::Instant::now();
            log::info!("Timing {label}: {} ms ({n} items)", now.duration_since(self.last).as_millis());
            self.last = now;
        }
    }
    fn total(&self) {
        if self.on {
            log::info!("Timing total: {} ms", self.start.elapsed().as_millis());
        }
    }
}

/// Read a pepXML file into memory.
fn read_file(path: &Path) -> Result<String> {
    std::fs::read_to_string(path)
        .map_err(|e| PtmProphetRsError::Io { path: path.to_path_buf(), source: e })
}

/// Read the configured single input pepXML (used by [`read_em_evidence`] and the test/validation
/// helpers, which drive off `config.in_file`).
fn read_input(config: &Config) -> Result<String> {
    read_file(&config.in_file)
}

/// A path's file name for progress messages (falls back to the whole path).
fn display_name(path: &Path) -> String {
    path.file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// PASS A — scan every input pepXML once, returning the combined [`InputIndex`] (mod tables +
/// probability hash merged across all files) and, per file in input order, its probability-gated
/// eligible PSMs with their mzML paths remapped onto any file-list overrides.
///
/// Files are independent, so the read + XML parse + PSM collection runs **in parallel across files**
/// (one of the few wholly-serial phases otherwise). The per-file results are then merged serially in
/// file order, keeping the merge deterministic and the static-mod conflict check intact.
fn scan_all_inputs(config: &Config) -> Result<(InputIndex, Vec<Vec<PsmRecord>>)> {
    // Per file (parallel): read, parse the index, collect eligible PSMs, remap mzML paths.
    let scanned: Vec<Result<(InputIndex, Vec<PsmRecord>)>> = par_map(&config.inputs, |spec| {
        let content = read_file(&spec.pep_xml)?;
        // Gate this file's PSMs by this file's own probabilities (the keys are unique per fraction,
        // so the merged hash would work too — but the per-file index keeps it self-contained).
        let file_index = crate::pepxml::index::parse_index_from_content(&content)?;
        let mut psms = collect_eligible_psms(config, &file_index, &content);
        // The spectrum file comes solely from the file-list; stamp it onto every PSM (the pepXML's
        // own base_name/raw_data is never consulted).
        let mzml = spec.mzml.to_string_lossy().into_owned();
        for p in &mut psms {
            p.data_file = mzml.clone();
        }
        Ok((file_index, psms))
    });

    // Merge serially in file order (deterministic; the static-mod conflict check fires here).
    let mut index = InputIndex::default();
    let mut per_file = Vec::with_capacity(config.inputs.len());
    for result in scanned {
        let (file_index, psms) = result?;
        crate::pepxml::index::merge_index(&mut index, file_index)?;
        per_file.push(psms);
    }
    Ok((index, per_file))
}

/// Every rank-1 modified PSM passing the `--min-prob` gate, in file order.
fn collect_eligible_psms(config: &Config, index: &InputIndex, content: &str) -> Vec<PsmRecord> {
    let mut psms = Vec::new();
    for_each_psm(content, |psm| {
        let prob = index.spec_prob.get(&psm.prob_key()).copied().unwrap_or(-1.0);
        if prob >= config.min_prob {
            psms.push(psm.clone());
        }
    });
    psms
}

/// Build an [`Mpx`] for each eligible PSM, keeping only those with ambiguous work, in input order.
///
/// `Mpx::build` (peptide parse + candidate-site enumeration) is non-trivial per PSM, so it runs in
/// parallel across the file's PSMs. Each PSM is moved out of `psms` via `mem::take` (no clone); the
/// returned items preserve input order (load-bearing for the in-order EM harvest and ROC tables).
fn build_work_items<'a>(
    config: &'a Config,
    index: &InputIndex,
    mut psms: Vec<PsmRecord>,
) -> Vec<WorkItem<'a>> {
    par_map_mut(&mut psms, |slot| {
        let psm = std::mem::take(slot);
        let mpx = Mpx::build(
            config,
            &config.mods,
            &index.mods,
            &psm.pep_seq,
            &psm.pep_mods,
            psm.charge,
            psm.calc_neut_mass,
        );
        mpx.has_work().then_some(WorkItem { psm, mpx, processed: false })
    })
    .into_iter()
    .flatten()
    .collect()
}

/// Decode every needed spectrum once, grouped by run file (each file decoded in parallel). A missing
/// `_calibrated.mzML` is replaced as [`resolve_spectrum_file`] describes. A run file that cannot be
/// opened is a fatal error (fail-fast), but an individual scan missing from an opened file is
/// tolerated (omitted, then skipped downstream).
fn load_store(items: &[WorkItem]) -> Result<SpectrumStore> {
    let mut needed: HashMap<&str, Vec<(i64, &str)>> = HashMap::new();
    for it in items {
        needed.entry(it.psm.data_file.as_str()).or_default().push((it.psm.scan, it.psm.native_id.as_str()));
    }
    let mut store = SpectrumStore::with_capacity(needed.len());
    for (path, requests) in needed {
        let file = resolve_spectrum_file(Path::new(path))?;
        let reader = MzmlReader::open(&file)?;
        let spectra = reader.load_scans(&requests)?;
        let missing = requests.iter().filter(|(scan, _)| !spectra.contains_key(scan)).count();
        if missing > 0 {
            log::warn!(
                "{} of {} PSMs have no matching spectrum in {}; they are written without localization",
                missing,
                requests.len(),
                file.display()
            );
        }
        store.insert(path.to_string(), spectra);
    }
    Ok(store)
}

/// The suffix FragPipe gives the mzML of MSFragger's calibrated spectra, matched in any case.
const CALIBRATED_SUFFIX: &str = "_calibrated.mzML";

/// The mzML to read for the file-list's spectrum file `listed`.
///
/// MSFragger writes no `<base>_calibrated.mzML` for a run with too few PSMs to calibrate, yet
/// FragPipe lists one for every run when calibration is on. When such a listed file does not exist,
/// this returns `<base>.mzML` or else `<base>_uncalibrated.mzML`, the files FragPipe's pepXML
/// rewriting (`RewritePepxml`, `PercolatorOutputToPepXML`) falls back to in the same order, and
/// fails if neither exists. Any other path is returned as listed.
fn resolve_spectrum_file(listed: &Path) -> Result<PathBuf> {
    if listed.exists() {
        return Ok(listed.to_path_buf());
    }
    let Some(base) = listed.to_str().and_then(|s| {
        let cut = s.len().checked_sub(CALIBRATED_SUFFIX.len())?;
        s.get(cut..).filter(|tail| tail.eq_ignore_ascii_case(CALIBRATED_SUFFIX)).map(|_| &s[..cut])
    }) else {
        return Ok(listed.to_path_buf());
    };
    let fallbacks = [PathBuf::from(format!("{base}.mzML")), PathBuf::from(format!("{base}_uncalibrated.mzML"))];
    match fallbacks.iter().find(|p| p.exists()) {
        Some(found) => {
            log::warn!(
                "Calibrated mzML {} does not exist (MSFragger writes none for a run with too few PSMs to calibrate); reading spectra from {} instead",
                listed.display(),
                found.display()
            );
            Ok(found.clone())
        }
        None => Err(PtmProphetRsError::Mzml(format!(
            "{} does not exist, and neither does {} or {}",
            listed.display(),
            fallbacks[0].display(),
            fallbacks[1].display()
        ))),
    }
}

/// Match each PSM's fragment ladder to its spectrum exactly once (parallel). Items whose spectrum is
/// absent are left `processed = false`.
fn process_all(items: &mut [WorkItem], store: &SpectrumStore) {
    par_for_each_mut(items, |it| {
        if let Some(spec) = store.get(it.psm.data_file.as_str()).and_then(|m| m.get(&it.psm.scan)) {
            it.mpx.process(spec);
            it.processed = true;
        }
    });
}

/// Append one file's ambiguous-site O-scores into the shared EM evidence, serially in input order so
/// the accumulation is independent of the worker count. Called once per file to pool across a batch.
fn harvest_into(items: &[WorkItem], ev: &mut EmEvidence) {
    for it in items {
        if it.processed {
            it.mpx.harvest_evidence(ev);
        }
    }
}

/// Harvest one file's evidence into a fresh [`EmEvidence`] (single-file convenience for
/// [`read_em_evidence`] and tests).
fn harvest_in_order(items: &[WorkItem]) -> EmEvidence {
    let mut ev = EmEvidence::default();
    harvest_into(items, &mut ev);
    ev
}

/// Per-PSM output contribution: `(scan, block, roc-contributions)` or `None` when the PSM does not
/// localize. Computed in parallel; the returned `Vec` preserves input order.
type Contribution = Option<(i64, String, Vec<(usize, RocInfoLocs)>)>;

/// Apply the trained model to every processed PSM and build its output block, in parallel; the
/// returned `Vec` preserves input order.
fn build_contributions(
    items: &mut [WorkItem],
    config: &Config,
    index: &InputIndex,
    model: &KdModel,
) -> Vec<Contribution> {
    par_map_mut(items, |it| {
        if !it.processed {
            return None;
        }
        it.mpx.compute_probs(model);
        crate::pepxml::writer::build_psm_block(config, index, &it.psm, &it.mpx)
            .map(|(block, rocs)| (it.psm.scan, block, rocs))
    })
}

/// First-pass EM evidence collector — process every eligible PSM once and harvest in order.
/// Equivalent to passes A→B of [`run`]; retained for the validation harness and tests.
pub fn read_em_evidence(config: &Config, index: &InputIndex) -> Result<EmEvidence> {
    let content = read_input(config)?;
    let mut psms = collect_eligible_psms(config, index, &content);
    // Single-file convenience: the spectrum file comes from the (sole) file-list entry, mirroring the
    // batch path in `scan_all_inputs`.
    if let Some(spec) = config.inputs.first() {
        let mzml = spec.mzml.to_string_lossy().into_owned();
        for p in &mut psms {
            p.data_file = mzml.clone();
        }
    }
    let mut items = build_work_items(config, index, psms);
    let store = load_store(&items)?;
    process_all(&mut items, &store);
    Ok(harvest_in_order(&items))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::masses::IonKind;
    use crate::model::modspec::ModType;
    use std::path::PathBuf;

    fn f05_config() -> Config {
        let pep_xml = PathBuf::from("I:/test_ptmprophet_rs/ccrcc_phospho/ptmprophet/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/interact-01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05.pep.xml");
        let mzml = PathBuf::from("F:/data/ccrcc_dda_phospho/mzml/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05_calibrated.mzML");
        Config {
            mods: vec![
                ModType { aas: vec!['M'], shift: 15.9949, ptm_label: "M:15.9949".into() },
                ModType { aas: vec!['S', 'T', 'Y'], shift: 79.96633, ptm_label: "STY:79.96633".into() },
            ],
            mod_string: "M:15.9949,STY:79.96633".into(),
            opt_string: String::new(),
            frag_tol: 10.0, em: 1, min_prob: 0.5,
            mod_prec: 0, min_o: 0, max_frag_z: 0, threads: None, nostack: true,
            static_tol: true,
            nions: vec![IonKind::B], cions: vec![IonKind::Y],
            inputs: vec![InputSpec { pep_xml: pep_xml.clone(), mzml, out: PathBuf::new() }],
            in_file: pep_xml,
            out_file: PathBuf::new(),
        }
    }

    /// End-to-end EM validation: train on all of f05, then re-score scan 3185 and compare the final
    /// localization probability to the oracle (0.3032 / 0.6968). Heavy (reads the full mzML); skips
    /// if the test data is absent.
    #[test]
    fn em_probability_scan3185_vs_oracle() {
        let cfg = f05_config();
        if !cfg.in_file.exists() || !cfg.inputs[0].mzml.exists() {
            eprintln!("[skip] test data not present");
            return;
        }
        let index = crate::pepxml::index::parse_spec_probs(&cfg).unwrap();
        let ev = read_em_evidence(&cfg, &index).unwrap();
        eprintln!("[em] harvested {} ambiguous sites", ev.len());
        assert!(ev.len() > 100, "expected many ambiguous sites, got {}", ev.len());

        let mut model = KdModel::new("PTMProphet_Oscore");
        crate::em::train::compute_ptm_model(&cfg, &mut model, &ev).unwrap();
        assert!(model.ready);

        // Re-score scan 3185.
        let content = std::fs::read_to_string(&cfg.in_file).unwrap();
        let mzml = cfg.inputs[0].mzml.to_string_lossy().into_owned();
        let mut rec: Option<PsmRecord> = None;
        for_each_psm(&content, |p| {
            if p.scan == 3185 {
                let mut p = p.clone();
                p.data_file = mzml.clone();
                rec = Some(p);
            }
        });
        let psm = rec.expect("scan 3185");
        let reader = MzmlReader::open(Path::new(&psm.data_file)).unwrap();
        let spec = reader.spectrum_for_scan(3185).unwrap();
        let mut mpx = Mpx::build(&cfg, &cfg.mods, &index.mods, &psm.pep_seq, &psm.pep_mods, psm.charge, psm.calc_neut_mass);
        mpx.process(&spec);
        mpx.compute_probs(&model);

        let sty = &mpx.types[1];
        let p5 = sty.scores[&4].prob;
        let p7 = sty.scores[&6].prob;
        eprintln!(
            "[scan3185] prob pos5={p5:.4} pos7={p7:.4} (oracle 0.3032 / 0.6968; Δ={:.4})",
            (p5 - 0.3032).abs()
        );
        // Probabilities sum to nmods (=1), best-site call matches, and the values match the oracle
        // to < 1e-3 (the numerically-faithful target).
        assert!((p5 + p7 - 1.0).abs() < 1e-6);
        assert!(p7 > p5 && p7 > 0.5, "best site should be pos7");
        assert!((p5 - 0.3032).abs() < 1e-3, "pos5 prob {p5} (oracle 0.3032)");
        assert!((p7 - 0.6968).abs() < 1e-3, "pos7 prob {p7} (oracle 0.6968)");
    }

    /// A fresh, empty directory for one test's files.
    fn empty_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("ptmprophet-rs-{name}-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).expect("create test directory");
        dir
    }

    fn touch(dir: &Path, name: &str) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, "").expect("create test file");
        path
    }

    #[test]
    fn missing_calibrated_mzml_falls_back_like_fragpipe() {
        let dir = empty_dir("resolve-spectrum-file");

        // A listed file that exists is read as listed.
        let present = touch(&dir, "a_calibrated.mzML");
        touch(&dir, "a.mzML");
        assert_eq!(resolve_spectrum_file(&present).unwrap(), present);

        // A missing calibrated file: <base>.mzML first, then <base>_uncalibrated.mzML.
        let both = touch(&dir, "b.mzML");
        touch(&dir, "b_uncalibrated.mzML");
        assert_eq!(resolve_spectrum_file(&dir.join("b_calibrated.mzML")).unwrap(), both);
        let uncalibrated = touch(&dir, "c_uncalibrated.mzML");
        assert_eq!(resolve_spectrum_file(&dir.join("c_calibrated.mzML")).unwrap(), uncalibrated);
        let other_case = touch(&dir, "d.mzML");
        assert_eq!(resolve_spectrum_file(&dir.join("d_CALIBRATED.MZML")).unwrap(), other_case, "suffix in any case");

        // No fallback exists: the error names every path tried.
        let listed = dir.join("e_calibrated.mzML");
        let err = resolve_spectrum_file(&listed).unwrap_err().to_string();
        for tried in [listed.clone(), dir.join("e.mzML"), dir.join("e_uncalibrated.mzML")] {
            assert!(err.contains(&tried.display().to_string()), "{err} names {}", tried.display());
        }

        // Any other missing file is returned as listed, so opening it fails with its own path.
        touch(&dir, "f_uncalibrated.mzML");
        let plain = dir.join("f.mzML");
        assert_eq!(resolve_spectrum_file(&plain).unwrap(), plain);
        assert_eq!(resolve_spectrum_file(&dir.join("calibrated.mzML")).unwrap(), dir.join("calibrated.mzML"));

        std::fs::remove_dir_all(&dir).ok();
    }
}
