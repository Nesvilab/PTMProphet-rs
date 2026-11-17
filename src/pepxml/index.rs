// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! PASS A — `parse_spec_probs`: one linear scan of the input pepXML.
//!
//! Builds (1) the search-summary static/variable modification tables and (2) the per-PSM probability
//! hash used for the `--min-prob` gate (keyed `experiment_label + spectrum`, **first-wins** on a
//! duplicate key).
//!
//! The input pepXML is streamed linearly in each pass (no byte-offset index into the
//! `spectrum_query`s); the probability is looked up by key. Parallelism happens later — over
//! spectra and PSMs — not over byte ranges of this file.

use std::collections::HashMap;

use crate::config::Config;
use crate::error::{PtmProphetRsError, Result};
use crate::pepxml::tag::Tag;

/// Static/variable modification tables harvested from `<search_summary>`.
#[derive(Debug, Default, Clone)]
pub struct ModTables {
    /// Residue letter (or `'n'`/`'c'` peptide terminus) → static mass shift.
    pub stat_mods: HashMap<char, f64>,
    /// Residue letter (or `'n'`/`'c'`) → variable mass shifts.
    pub var_mods: HashMap<char, Vec<f64>>,
    /// Protein-terminus static mods, keyed `'n'`/`'c'`.
    pub stat_prot_termods: HashMap<char, f64>,
    /// Protein-terminus variable mods, keyed `'n'`/`'c'`.
    pub var_prot_termods: HashMap<char, Vec<f64>>,
}

/// Everything PASS A produces.
#[derive(Debug, Default)]
pub struct InputIndex {
    /// key = `experiment_label + spectrum` (no separator); value = PSM prob (or -1).
    pub spec_prob: HashMap<String, f64>,
    pub mods: ModTables,
    pub num_spectrum_queries: usize,
}

/// C-locale `atof` equivalent: parse a leading numeric prefix, ignoring trailing junk.
pub fn atof(s: &str) -> f64 {
    let s = s.trim_start();
    let mut end = 0;
    let bytes = s.as_bytes();
    let mut seen_dot = false;
    let mut seen_e = false;
    while end < bytes.len() {
        let c = bytes[end] as char;
        let ok = match c {
            '0'..='9' => true,
            '+' | '-' => end == 0 || (seen_e && (bytes[end - 1] as char == 'e' || bytes[end - 1] as char == 'E')),
            '.' if !seen_dot && !seen_e => {
                seen_dot = true;
                true
            }
            'e' | 'E' if !seen_e && end > 0 => {
                seen_e = true;
                true
            }
            _ => false,
        };
        if !ok {
            break;
        }
        end += 1;
    }
    s[..end].parse().unwrap_or(0.0)
}

/// C-locale `atoi` equivalent.
pub fn atoi(s: &str) -> i64 {
    let s = s.trim_start();
    let mut end = 0;
    let bytes = s.as_bytes();
    while end < bytes.len() {
        let c = bytes[end] as char;
        if c.is_ascii_digit() || ((c == '+' || c == '-') && end == 0) {
            end += 1;
        } else {
            break;
        }
    }
    s[..end].parse().unwrap_or(0)
}

/// Iterate every `<...>` tag in `content`, calling `f` for each (line-by-line, multiple tags per
/// line). Shared by PASS A here and the read/write passes.
pub fn for_each_tag(content: &str, mut f: impl FnMut(&Tag)) {
    for line in content.lines() {
        let mut rest = line;
        while let Some(lt) = rest.find('<') {
            let from = &rest[lt..];
            if let Some(tag) = Tag::parse(from) {
                f(&tag);
            }
            // advance past this '<'
            rest = &from[1..];
        }
    }
}

/// Scan the configured input pepXML once and build the [`InputIndex`] (reads `config.in_file`).
/// Convenience wrapper over [`parse_index_from_content`] for the single-file path and the
/// test/validation helpers; the batch pipeline reads each file itself and merges with [`merge_index`].
pub fn parse_spec_probs(config: &Config) -> Result<InputIndex> {
    let content = std::fs::read_to_string(&config.in_file).map_err(|e| PtmProphetRsError::Io {
        path: config.in_file.clone(),
        source: e,
    })?;
    parse_index_from_content(&content)
}

/// Merge one file's [`InputIndex`] (`from`) into an accumulator (`into`) for batch runs.
///
/// Modification tables are unioned: static mods are conflict-checked (same residue, differing mass →
/// error, as within a single file), variable/terminal mods deduped by ~equality. The probability
/// hash is extended first-wins (keys are `experiment_label + spectrum`, unique across fractions), and
/// the run list and query count are accumulated.
pub fn merge_index(into: &mut InputIndex, from: InputIndex) -> Result<()> {
    for (aa, m) in from.mods.stat_mods {
        if let Some(&prev) = into.mods.stat_mods.get(&aa) {
            if (prev - m).abs() > 0.001 {
                return Err(PtmProphetRsError::StaticModConflict { residue: aa, m1: prev, m2: m });
            }
        }
        into.mods.stat_mods.insert(aa, m);
    }
    for (aa, m) in from.mods.stat_prot_termods {
        if let Some(&prev) = into.mods.stat_prot_termods.get(&aa) {
            if (prev - m).abs() > 0.001 {
                return Err(PtmProphetRsError::StaticModConflict { residue: aa, m1: prev, m2: m });
            }
        }
        into.mods.stat_prot_termods.insert(aa, m);
    }
    for (aa, shifts) in from.mods.var_mods {
        let v = into.mods.var_mods.entry(aa).or_default();
        for s in shifts {
            if !v.iter().any(|&x| (x - s).abs() < 1e-5) {
                v.push(s);
            }
        }
    }
    for (aa, shifts) in from.mods.var_prot_termods {
        let v = into.mods.var_prot_termods.entry(aa).or_default();
        for s in shifts {
            if !v.iter().any(|&x| (x - s).abs() < 1e-5) {
                v.push(s);
            }
        }
    }
    for (k, p) in from.spec_prob {
        into.spec_prob.entry(k).or_insert(p);
    }
    into.num_spectrum_queries += from.num_spectrum_queries;
    Ok(())
}

/// Scan an already-read pepXML once and build the [`InputIndex`].
pub fn parse_index_from_content(content: &str) -> Result<InputIndex> {
    let mut idx = InputIndex::default();

    // Per-spectrum_query transient state.
    let mut exp_lbl = String::new();
    let mut spectrum_name = String::new();
    let mut spec_prob = -1.0_f64;
    let mut in_query = false;

    let mut err: Option<PtmProphetRsError> = None;

    for_each_tag(content, |tag| {
        if err.is_some() {
            return;
        }
        let name = tag.name.as_str();
        match name {
            "aminoacid_modification" => {
                let aa = tag.get("aminoacid").and_then(|s| s.chars().next());
                let massdiff = tag.get("massdiff").map(atof);
                let variable = tag.get("variable").unwrap_or("N");
                if let (Some(aa), Some(md)) = (aa, massdiff) {
                    if variable == "N" {
                        if let Some(&prev) = idx.mods.stat_mods.get(&aa) {
                            if (prev - md).abs() > 0.001 {
                                err = Some(PtmProphetRsError::StaticModConflict { residue: aa, m1: prev, m2: md });
                                return;
                            }
                        }
                        idx.mods.stat_mods.insert(aa, md);
                    } else {
                        let v = idx.mods.var_mods.entry(aa).or_default();
                        if !v.iter().any(|&x| (x - md).abs() < 1e-5) {
                            v.push(md);
                        }
                    }
                }
            }
            "terminal_modification" => {
                let term = tag.get("terminus").and_then(|s| s.chars().next()).map(|c| c.to_ascii_lowercase());
                let massdiff = tag.get("massdiff").map(atof);
                let variable = tag.get("variable").unwrap_or("N");
                let prot = tag.get("protein_terminus").unwrap_or("N");
                if let (Some(t), Some(md)) = (term, massdiff) {
                    match (variable, prot) {
                        ("N", "N") => {
                            idx.mods.stat_mods.insert(t, md);
                        }
                        ("Y", "N") => {
                            let v = idx.mods.var_mods.entry(t).or_default();
                            if !v.iter().any(|&x| (x - md).abs() < 1e-5) {
                                v.push(md);
                            }
                        }
                        ("N", "Y") => {
                            idx.mods.stat_prot_termods.insert(t, md);
                        }
                        _ => {
                            let v = idx.mods.var_prot_termods.entry(t).or_default();
                            if !v.iter().any(|&x| (x - md).abs() < 1e-5) {
                                v.push(md);
                            }
                        }
                    }
                }
            }
            "spectrum_query" if tag.start => {
                in_query = true;
                exp_lbl = tag.get("experiment_label").unwrap_or("").to_string();
                spectrum_name = tag.get("spectrum").unwrap_or("").to_string();
                spec_prob = -1.0;
                idx.num_spectrum_queries += 1;
            }
            "peptideprophet_result" if tag.start => {
                if let Some(p) = tag.get("probability") {
                    spec_prob = atof(p);
                }
            }
            "interprophet_result" if tag.start => {
                if let Some(p) = tag.get("probability") {
                    spec_prob = atof(p);
                }
            }
            "spectrum_query" if tag.end && !tag.start => {
                if in_query {
                    let key = format!("{exp_lbl}{spectrum_name}");
                    // first-wins: insert only if absent.
                    idx.spec_prob.entry(key).or_insert(spec_prob);
                }
                in_query = false;
                exp_lbl.clear();
                spectrum_name.clear();
                spec_prob = -1.0;
            }
            _ => {}
        }
    });

    if let Some(e) = err {
        return Err(e);
    }
    Ok(idx)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::masses::IonKind;
    use std::path::PathBuf;

    #[test]
    fn atof_lenient() {
        assert!((atof("79.96633") - 79.96633).abs() < 1e-9);
        assert!((atof("79.96633x") - 79.96633).abs() < 1e-9);
        assert!((atof("  -18.010565") + 18.010565).abs() < 1e-9);
        assert!((atof("3.0315e-01") - 0.30315).abs() < 1e-9);
        assert_eq!(atof("abc"), 0.0);
    }

    #[test]
    fn merge_index_unions_mods_and_probs() {
        // Two "fractions" with identical static mods but per-fraction probabilities/runs.
        let a = parse_index_from_content(
            "<aminoacid_modification aminoacid=\"C\" massdiff=\"57.02146\" variable=\"N\"/>\
             <aminoacid_modification aminoacid=\"S\" massdiff=\"79.96633\" variable=\"Y\"/>\
             <msms_run_summary base_name=\"/d/f01\" raw_data=\"mzML\">\
             <spectrum_query experiment_label=\"f01\" spectrum=\"f01.1.1.2\">\
             <peptideprophet_result probability=\"0.9\"></spectrum_query>",
        )
        .unwrap();
        let b = parse_index_from_content(
            "<aminoacid_modification aminoacid=\"C\" massdiff=\"57.02146\" variable=\"N\"/>\
             <aminoacid_modification aminoacid=\"T\" massdiff=\"79.96633\" variable=\"Y\"/>\
             <msms_run_summary base_name=\"/d/f02\" raw_data=\"mzML\">\
             <spectrum_query experiment_label=\"f02\" spectrum=\"f02.1.1.2\">\
             <peptideprophet_result probability=\"0.8\"></spectrum_query>",
        )
        .unwrap();

        let mut merged = a;
        merge_index(&mut merged, b).unwrap();
        assert!((merged.mods.stat_mods[&'C'] - 57.02146).abs() < 1e-6);
        assert!(merged.mods.var_mods.contains_key(&'S') && merged.mods.var_mods.contains_key(&'T'));
        assert_eq!(merged.spec_prob.len(), 2);
        assert!((merged.spec_prob["f01f01.1.1.2"] - 0.9).abs() < 1e-6);
        assert!((merged.spec_prob["f02f02.1.1.2"] - 0.8).abs() < 1e-6);
    }

    #[test]
    fn merge_index_rejects_conflicting_static_mod() {
        let mut a = parse_index_from_content(
            "<aminoacid_modification aminoacid=\"C\" massdiff=\"57.02146\" variable=\"N\"/>",
        )
        .unwrap();
        let b = parse_index_from_content(
            "<aminoacid_modification aminoacid=\"C\" massdiff=\"58.00000\" variable=\"N\"/>",
        )
        .unwrap();
        assert!(merge_index(&mut a, b).is_err());
    }

    #[test]
    fn parse_real_f05_pepxml() {
        let path = PathBuf::from("I:/test_ptmprophet_rs/ccrcc_phospho/ptmprophet/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/interact-01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05.pep.xml");
        if !path.exists() {
            eprintln!("[skip] test data not present: {}", path.display());
            return;
        }
        let cfg = Config {
            in_file: path,
            out_file: PathBuf::new(),
            ..test_config()
        };
        let idx = parse_spec_probs(&cfg).expect("parse_spec_probs");

        // Static mods from the search summary (TMT run).
        assert!((idx.mods.stat_mods[&'C'] - 57.02146).abs() < 1e-4, "C static {:?}", idx.mods.stat_mods.get(&'C'));
        assert!((idx.mods.stat_mods[&'K'] - 229.16293).abs() < 1e-4);
        assert!((idx.mods.stat_mods[&'n'] - 229.16293).abs() < 1e-4, "n-term TMT {:?}", idx.mods.stat_mods.get(&'n'));
        // Variable phospho on S/T/Y and oxidation on M.
        assert!(idx.mods.var_mods[&'S'].iter().any(|&x| (x - 79.96633).abs() < 1e-4));
        assert!(idx.mods.var_mods[&'T'].iter().any(|&x| (x - 79.96633).abs() < 1e-4));
        assert!(idx.mods.var_mods[&'M'].iter().any(|&x| (x - 15.9949).abs() < 1e-4));

        // Probabilities: many queries; the known scan 3185 PSM has prob >= 0.5.
        assert!(idx.num_spectrum_queries > 1000, "queries={}", idx.num_spectrum_queries);
        assert_eq!(idx.spec_prob.len(), idx.num_spectrum_queries, "unique spectrum keys");
        let k = "01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05.03185.03185.2";
        let p = idx.spec_prob.get(k).copied().unwrap_or(-1.0);
        assert!(p >= 0.5, "scan 3185 prob = {p}");
        eprintln!("[ok] f05: {} queries, scan3185 prob={p}", idx.num_spectrum_queries);
    }

    fn test_config() -> Config {
        Config {
            mods: Vec::new(),
            mod_string: String::new(),
            opt_string: String::new(),
            frag_tol: 10.0,
            em: 1,
            min_prob: 0.5,
            mod_prec: 0,
            min_o: 0,
            max_frag_z: 0,
            threads: None,
            nostack: true,
            static_tol: true,
            nions: vec![IonKind::B],
            cions: vec![IonKind::Y],
            inputs: Vec::new(),
            in_file: PathBuf::new(),
            out_file: PathBuf::new(),
        }
    }
}
