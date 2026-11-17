// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! Streaming PSM reader — shared by the EM read pass and the localizing write pass.
//!
//! Walks the pepXML linearly and yields one [`PsmRecord`] per `spectrum_query` whose rank-1
//! `search_hit` carries a `<modification_info>`. Each record carries everything needed to build an
//! [`crate::mpx::Mpx`]: the stripped peptide, the position→absolute-mass map (`pep_mods_hash`
//! convention: key 0 = N-term, `len+1` = C-term, 1..len = residues), charge, calc neutral mass,
//! scan number, and the probability-hash key. The mzML path is assigned later from the `--filelist`,
//! never read from the pepXML.

use std::collections::BTreeMap;

use crate::pepxml::index::{atof, atoi, for_each_tag};

/// One PSM (rank-1 hit) ready for localization.
#[derive(Debug, Clone, Default)]
pub struct PsmRecord {
    pub exp_label: String,
    pub spectrum: String,
    pub scan: i64,
    pub charge: i32,
    pub calc_neut_mass: f64,
    pub pep_seq: String,
    /// position → absolute residue+mod mass (0 = N-term, len+1 = C-term).
    pub pep_mods: BTreeMap<usize, f64>,
    /// mzML path for this PSM, assigned from the `--filelist` (never read from the pepXML).
    pub data_file: String,
}

impl PsmRecord {
    /// The `spec_prob` hash key (`experiment_label + spectrum`).
    pub fn prob_key(&self) -> String {
        format!("{}{}", self.exp_label, self.spectrum)
    }
}

/// Iterate every rank-1, modified PSM in `content`, calling `f` for each.
pub fn for_each_psm(content: &str, mut f: impl FnMut(&PsmRecord)) {
    let mut cur = PsmRecord::default();
    let mut in_query = false;
    let mut hit_rank = 0i64;
    let mut sh_go = false; // inside the rank-1 search_hit
    let mut in_mod = false;
    let mut has_mod = false;

    for_each_tag(content, |tag| {
        let name = tag.name.as_str();
        match name {
            "spectrum_query" if tag.start => {
                cur = PsmRecord::default();
                cur.exp_label = tag.get("experiment_label").unwrap_or("").to_string();
                cur.spectrum = tag.get("spectrum").unwrap_or("").to_string();
                cur.charge = tag.get("assumed_charge").map(atoi).unwrap_or(0) as i32;
                cur.scan = tag.get("start_scan").map(atoi).unwrap_or(0);
                in_query = true;
                hit_rank = 0;
                sh_go = false;
                in_mod = false;
                has_mod = false;
            }
            "search_hit" if tag.start => {
                hit_rank = tag.get("hit_rank").map(atoi).unwrap_or(1);
                if hit_rank == 1 {
                    sh_go = true;
                    cur.pep_seq = tag.get("peptide").unwrap_or("").to_string();
                    cur.calc_neut_mass = tag.get("calc_neutral_pep_mass").map(atof).unwrap_or(0.0);
                    cur.pep_mods.clear();
                } else {
                    sh_go = false;
                }
            }
            "modification_info" if tag.start && sh_go => {
                in_mod = true;
                has_mod = true;
                let len = cur.pep_seq.chars().count();
                if let Some(m) = tag.get("mod_nterm_mass") {
                    cur.pep_mods.insert(0, atof(m));
                }
                if let Some(m) = tag.get("mod_cterm_mass") {
                    cur.pep_mods.insert(len + 1, atof(m));
                }
            }
            "mod_aminoacid_mass" if in_mod && sh_go => {
                if let (Some(p), Some(m)) = (tag.get("position"), tag.get("mass")) {
                    cur.pep_mods.insert(atoi(p) as usize, atof(m));
                }
            }
            "modification_info" if tag.end && !tag.start => {
                in_mod = false;
            }
            "spectrum_query" if tag.end && !tag.start => {
                if in_query && has_mod && !cur.pep_seq.is_empty() {
                    f(&cur);
                }
                in_query = false;
            }
            _ => {}
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn psm_stream_real_f05() {
        let path = Path::new("I:/test_ptmprophet_rs/ccrcc_phospho/ptmprophet/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/interact-01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05.pep.xml");
        if !path.exists() {
            eprintln!("[skip] test data not present");
            return;
        }
        let content = std::fs::read_to_string(path).unwrap();
        let mut count = 0usize;
        let mut scan3185: Option<PsmRecord> = None;
        for_each_psm(&content, |psm| {
            count += 1;
            if psm.scan == 3185 {
                scan3185 = Some(psm.clone());
            }
        });
        assert!(count > 1000, "modified PSMs = {count}");
        let p = scan3185.expect("scan 3185");
        assert_eq!(p.pep_seq, "ALRASESGI");
        assert_eq!(p.charge, 2);
        assert!((p.pep_mods[&0] - 230.17076).abs() < 1e-4); // n-term TMT
        assert!((p.pep_mods[&5] - 166.99835).abs() < 1e-4); // S5 phospho
        assert!(p.data_file.is_empty(), "data_file is assigned from the file-list, not the pepXML: {}", p.data_file);
        eprintln!("[ok] {count} modified PSMs; scan3185 = {} z{}", p.pep_seq, p.charge);
    }
}
