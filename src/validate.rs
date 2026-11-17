// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Validation helpers — parse a reference `interact-*.mod.pep.xml` and compare PTMProphet-rs's
//! per-site scores against it. Used by the diff harness and the diagnostic tests.

use std::collections::HashMap;

use crate::pepxml::index::{atof, atoi, for_each_tag};

/// One reference site result parsed from `<mod_aminoacid_probability>`.
#[derive(Debug, Clone, Copy, Default)]
pub struct OracleSite {
    pub probability: f64,
    pub oscore: f64,
    pub mscore: f64,
    pub cterm_score: f64,
    pub nterm_score: f64,
}

/// key = (scan, ptm-label, 1-based position) → site scores.
pub type OracleMap = HashMap<(i64, String, usize), OracleSite>;

/// Parse a reference `.mod.pep.xml` into per-(scan, ptm, position) site scores.
pub fn parse_oracle(content: &str) -> OracleMap {
    let mut map = OracleMap::new();
    let mut scan: i64 = 0;
    let mut ptm = String::new();
    for_each_tag(content, |tag| match tag.name.as_str() {
        "spectrum_query" if tag.start => {
            scan = tag.get("start_scan").map(atoi).unwrap_or(0);
        }
        "ptmprophet_result" if tag.start => {
            ptm = tag.get("ptm").unwrap_or("").to_string();
        }
        "mod_aminoacid_probability" => {
            if let Some(pos) = tag.get("position").map(atoi) {
                let site = OracleSite {
                    probability: tag.get("probability").map(atof).unwrap_or(0.0),
                    oscore: tag.get("oscore").map(atof).unwrap_or(0.0),
                    mscore: tag.get("mscore").map(atof).unwrap_or(0.0),
                    cterm_score: tag.get("cterm_score").map(atof).unwrap_or(0.0),
                    nterm_score: tag.get("nterm_score").map(atof).unwrap_or(0.0),
                };
                map.insert((scan, ptm.clone(), pos as usize), site);
            }
        }
        _ => {}
    });
    map
}

/// Running comparison statistics for one scalar attribute.
#[derive(Debug, Default, Clone)]
pub struct DiffStats {
    pub n: usize,
    pub within_1e3: usize,
    pub within_1e2: usize,
    pub max_abs: f64,
    pub sum_abs: f64,
    pub worst: Option<(i64, usize, f64, f64)>, // (scan, pos, mine, oracle)
}

impl DiffStats {
    pub fn add(&mut self, scan: i64, pos: usize, mine: f64, oracle: f64) {
        let d = (mine - oracle).abs();
        self.n += 1;
        if d < 1e-3 {
            self.within_1e3 += 1;
        }
        if d < 1e-2 {
            self.within_1e2 += 1;
        }
        self.sum_abs += d;
        if d > self.max_abs {
            self.max_abs = d;
            self.worst = Some((scan, pos, mine, oracle));
        }
    }
    pub fn mean_abs(&self) -> f64 {
        if self.n == 0 {
            0.0
        } else {
            self.sum_abs / self.n as f64
        }
    }
    pub fn summary(&self, label: &str) -> String {
        format!(
            "{label}: n={} within1e-3={:.1}% within1e-2={:.1}% mean={:.5} max={:.5} worst={:?}",
            self.n,
            100.0 * self.within_1e3 as f64 / self.n.max(1) as f64,
            100.0 * self.within_1e2 as f64 / self.n.max(1) as f64,
            self.mean_abs(),
            self.max_abs,
            self.worst,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{Config, InputSpec};
    use crate::em::kdmodel::KdModel;
    use crate::model::modspec::ModType;
    use crate::mpx::Mpx;
    use crate::pepxml::psm::for_each_psm;
    use crate::model::masses::IonKind;
    use crate::spectrum::mzml::MzmlReader;
    use std::path::{Path, PathBuf};

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

    /// Diagnostic: compare PTMProphet-rs's per-site O-score and probability to the oracle across
    /// ALL of f05. Localizes the accuracy gap (task #12): if O-scores match everywhere, the gap is
    /// purely EM.
    #[test]
    fn diff_all_psms_vs_oracle() {
        let cfg = f05_config();
        let oracle_path = PathBuf::from(
            cfg.in_file.to_string_lossy().replace(".pep.xml", ".mod.pep.xml"),
        );
        if !cfg.in_file.exists() || !oracle_path.exists() || !cfg.inputs[0].mzml.exists() {
            eprintln!("[skip] test data not present");
            return;
        }
        let oracle = parse_oracle(&std::fs::read_to_string(&oracle_path).unwrap());
        eprintln!("[oracle] {} site records", oracle.len());

        let index = crate::pepxml::index::parse_spec_probs(&cfg).unwrap();
        // Train the EM model on all PSMs.
        let ev = crate::pipeline::read_em_evidence(&cfg, &index).unwrap();
        let mut model = KdModel::new("PTMProphet_Oscore");
        crate::em::train::compute_ptm_model(&cfg, &mut model, &ev).unwrap();

        // Second pass: score every PSM and compare to the oracle.
        let content = std::fs::read_to_string(&cfg.in_file).unwrap();
        let mzml = cfg.inputs[0].mzml.to_string_lossy().into_owned();
        let mut psms = Vec::new();
        for_each_psm(&content, |p| {
            let prob = index.spec_prob.get(&p.prob_key()).copied().unwrap_or(-1.0);
            if prob >= cfg.min_prob {
                let mut p = p.clone();
                p.data_file = mzml.clone();
                psms.push(p);
            }
        });

        let mut reader_cache: Option<(String, MzmlReader)> = None;
        let mut osc = DiffStats::default();
        let mut prb = DiffStats::default();
        let mut psm_done = 0usize;

        for psm in &psms {
            let mut mpx = Mpx::build(&cfg, &cfg.mods, &index.mods, &psm.pep_seq, &psm.pep_mods, psm.charge, psm.calc_neut_mass);
            if !mpx.has_work() {
                continue;
            }
            if reader_cache.as_ref().map(|(p, _)| p != &psm.data_file).unwrap_or(true) {
                match MzmlReader::open(Path::new(&psm.data_file)) {
                    Ok(r) => reader_cache = Some((psm.data_file.clone(), r)),
                    Err(_) => continue,
                }
            }
            let spec = match reader_cache.as_ref().unwrap().1.spectrum_for_scan(psm.scan) {
                Ok(s) => s,
                Err(_) => continue,
            };
            mpx.process(&spec);
            mpx.compute_probs(&model);
            psm_done += 1;

            for (t, mt) in cfg.mods.iter().enumerate() {
                for (&site, sc) in &mpx.types[t].scores {
                    let key = (psm.scan, mt.ptm_label.clone(), site + 1);
                    if let Some(o) = oracle.get(&key) {
                        osc.add(psm.scan, site + 1, sc.oscore, o.oscore);
                        prb.add(psm.scan, site + 1, sc.prob, o.probability);
                    }
                }
            }
        }

        eprintln!("[diff] scored {psm_done} PSMs");
        eprintln!("{}", osc.summary("OSCORE"));
        eprintln!("{}", prb.summary("PROB  "));
        // Acceptance: O-score and probability match the reference within 1e-3 on >98% of all
        // sites (the residual <2% are adjacent-site ties / boundary peaks).
        assert!(osc.n > 1000, "compared {} sites", osc.n);
        let osc_frac = osc.within_1e3 as f64 / osc.n as f64;
        let prb_frac = prb.within_1e3 as f64 / prb.n as f64;
        assert!(osc_frac > 0.98, "O-score <1e-3 fraction {osc_frac:.3}");
        assert!(prb_frac > 0.98, "probability <1e-3 fraction {prb_frac:.3}");
        assert!(prb.mean_abs() < 1e-3, "mean probability diff {}", prb.mean_abs());
    }
}
