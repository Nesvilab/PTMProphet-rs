// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! The localizing writer. Streams the input pepXML line-by-line, copies it
//! through, and inserts localization content:
//!   * after `<msms_pipeline_analysis ...>`: the `<analysis_summary analysis="ptmprophet">` header
//!     (params, the ROC / error-rate summary, and the EM `<mixturemodel>` curve);
//!   * in place of each rank-1 hit's `<modification_info>`: for a localized PSM, the rewritten
//!     `<modification_info>` with the variable mods on the best-localized sites, followed by the
//!     `<analysis_result analysis="ptmprophet">` block with `<ptmprophet_result>` +
//!     `<mod_aminoacid_probability>` rows; for any other hit, the original reformatted to
//!     PTMProphet's form (an empty one for an unmodified hit). See [`stream_output`].
//!
//! Output number formats: `prior` default precision; `ptm_peptide` probs fixed-4;
//! `probability` scientific; `oscore`/`mscore`/`*_term_score` fixed-4.
//!
//! The emitted file is valid pepXML carrying the per-site localization probabilities downstream
//! tools read.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write as _;
use std::path::Path;

use crate::config::Config;
use crate::em::kdmodel::KdModel;
use crate::error::{PtmProphetRsError, Result};
use crate::model::masses;
use crate::mpx::Mpx;
use crate::pepxml::index::InputIndex;
use crate::pepxml::psm::PsmRecord;
use crate::pepxml::roc::InfoLocs as RocInfoLocs;

/// The `version` attribute of `<ptmprophet_summary>`: PTMProphet-rs and its crate version. The
/// element names and `analysis="ptmprophet"` stay PTMProphet's, because Philosopher and other
/// pepXML readers look for them.
const SUMMARY_VERSION: &str =
    concat!("PTMProphet-rs ", env!("CARGO_PKG_VERSION"), " (Rust port of TPP PTMProphet)");

/// Build one PSM's output contribution from its **already processed and prob-scored** `Mpx`:
/// the rewritten `<modification_info>` (mods moved to localized sites) immediately followed by
/// the `<analysis_result analysis="ptmprophet">` block, plus the per-type ROC `InfoLocs` to
/// accumulate. Returns `None` if no type localizes. Pure — does no I/O — so it can be computed
/// for every PSM in parallel.
pub fn build_psm_block(
    config: &Config,
    index: &InputIndex,
    psm: &PsmRecord,
    mpx: &Mpx,
) -> Option<(String, Vec<(usize, RocInfoLocs)>)> {
    let ar = build_analysis_result(config, mpx)?;
    let modinfo = build_modification_info(config, index, psm, mpx);
    let block = format!("{modinfo}{ar}");

    let spec_prob = index.spec_prob.get(&psm.prob_key()).copied().unwrap_or(-1.0);
    let mut rocs = Vec::new();
    for (t, tr) in mpx.types.iter().enumerate() {
        if tr.terminus.is_some() {
            // Terminal mods are unambiguous → fully localized (info_gain = localized_mods = 1).
            rocs.push((t, RocInfoLocs::new(1.0, 1.0, 1.0, spec_prob as f32, 1.0)));
            continue;
        }
        if tr.sites.is_empty() || tr.nmods == 0 {
            continue;
        }
        let (nig, localized) = mpx.info_gain(t);
        if nig >= 0.0 {
            let mbp = mean_best_prob(tr);
            rocs.push((
                t,
                RocInfoLocs::new(localized as f32, nig as f32, mbp as f32, spec_prob as f32, tr.nmods as f32),
            ));
        }
    }
    Some((block, rocs))
}

/// Stream the localized output pepXML: copy the input through, splicing the
/// `<analysis_summary>` header after `<msms_pipeline_analysis ...>` and rewriting **every** rank-1
/// hit's `<modification_info>` to PTMProphet's canonical form, exactly as PTMProphet does — every
/// `search_hit` in its output carries a `<modification_info>`. Three cases:
///   * **localized** (scan in `blocks`): replace the original with the precomputed block (rewritten
///     `<modification_info>` + `<analysis_result>`), full `%.8f` precision.
///   * **modified but not localized** (static mods only, low-prob, or otherwise unscored): reformat
///     the original `<modification_info>` into the canonical form with `static`/`variable`/`source`
///     attributes and bracketed `modified_peptide`, using PTMProphet's `%.8g` "pass-through" stream
///     precision (see [`build_passthrough_modinfo`]).
///   * **unmodified** (no input `<modification_info>`): insert an empty
///     `<modification_info modified_peptide="SEQ"></modification_info>` before the first downstream
///     element, matching PTMProphet which emits one for every hit.
///
/// `blocks` is keyed by `start_scan`. `stat_mods` is the merged search-summary static-mod table (used
/// to decompose each pass-through site into residue + static + variable).
///
/// `in_path`/`out_path` are this input's own paths (each batch file has its own), so the header's
/// `<inputfile>` names the right file and the output lands beside it.
pub fn stream_output(
    config: &Config,
    stat_mods: &HashMap<char, f64>,
    in_path: &Path,
    out_path: &Path,
    content: &str,
    intens_model: &KdModel,
    roc_block: &str,
    blocks: &HashMap<i64, String>,
) -> Result<()> {
    let out = std::fs::File::create(out_path).map_err(|e| PtmProphetRsError::Io {
        path: out_path.to_path_buf(),
        source: e,
    })?;
    let mut w = std::io::BufWriter::new(out);

    let header = build_header(config, in_path, intens_model, roc_block);
    let mut cur_scan: i64 = 0;
    let mut hit_rank = 1i64;
    let mut pep_seq = String::new();
    let mut seen_modinfo = false; // a <modification_info> was emitted for the current rank-1 hit

    let mut skipping = false; // inside a localized input <modification_info> we are replacing
    // While set, we are buffering a non-localized input <modification_info> to reformat it.
    let mut pass_active = false;
    let mut pass_masses: Vec<(usize, f64)> = Vec::new();
    let mut pass_nterm: Option<f64> = None;
    let mut pass_cterm: Option<f64> = None;

    for line in content.lines() {
        let trimmed = line.trim_start();
        // Update PSM state from structural tags (these are never skipped).
        if trimmed.starts_with("<spectrum_query") {
            cur_scan = attr(line, "start_scan").map(|s| s.parse().unwrap_or(0)).unwrap_or(0);
            hit_rank = 1;
        } else if trimmed.starts_with("<search_hit") {
            hit_rank = attr(line, "hit_rank").map(|s| s.parse().unwrap_or(1)).unwrap_or(1);
            pep_seq = attr(line, "peptide").unwrap_or_default();
            seen_modinfo = false;
        }

        // Finishing a localized replacement: drop the original lines, emit the precomputed block.
        if skipping {
            if trimmed.starts_with("</modification_info>") {
                skipping = false;
                if let Some(block) = blocks.get(&cur_scan) {
                    // rewritten <modification_info> + <analysis_result>
                    w.write_all(block.as_bytes()).map_err(io_err(out_path))?;
                }
            }
            continue; // drop the original modification_info lines
        }

        // Buffering a non-localized <modification_info>: collect its sites, then reformat on close.
        if pass_active {
            if trimmed.starts_with("<mod_aminoacid_mass") {
                if let (Some(p), Some(m)) = (attr(line, "position"), attr(line, "mass")) {
                    pass_masses.push((p.parse().unwrap_or(0), m.parse().unwrap_or(0.0)));
                }
            } else if trimmed.starts_with("</modification_info>") {
                pass_active = false;
                let block = build_passthrough_modinfo(&pep_seq, &pass_masses, pass_nterm, pass_cterm, stat_mods);
                w.write_all(block.as_bytes()).map_err(io_err(out_path))?;
            }
            continue; // drop the original modification_info lines
        }

        // Start of a rank-1 hit's <modification_info>: either replace (localized) or reformat.
        if trimmed.starts_with("<modification_info") && hit_rank == 1 {
            seen_modinfo = true;
            if blocks.contains_key(&cur_scan) {
                skipping = true;
                continue; // drop the original <modification_info ...> start line
            }
            pass_active = true;
            pass_masses.clear();
            pass_nterm = attr(line, "mod_nterm_mass").and_then(|s| s.parse().ok());
            pass_cterm = attr(line, "mod_cterm_mass").and_then(|s| s.parse().ok());
            continue; // drop the original <modification_info ...> start line
        }

        // Unmodified rank-1 hit: PTMProphet still emits a (now empty) <modification_info>. Insert it
        // before the first element that follows where one would have been.
        if hit_rank == 1
            && !seen_modinfo
            && (trimmed.starts_with("<search_score")
                || trimmed.starts_with("<analysis_result")
                || trimmed.starts_with("</search_hit>"))
        {
            seen_modinfo = true;
            writeln!(w, "<modification_info modified_peptide=\"{pep_seq}\">").map_err(io_err(out_path))?;
            writeln!(w, "</modification_info>").map_err(io_err(out_path))?;
        }

        writeln!(w, "{line}").map_err(io_err(out_path))?;
        if trimmed.starts_with("<msms_pipeline_analysis") {
            w.write_all(header.as_bytes()).map_err(io_err(out_path))?;
        }
    }
    w.flush().map_err(io_err(out_path))?;
    Ok(())
}

fn io_err(out_path: &Path) -> impl Fn(std::io::Error) -> PtmProphetRsError + '_ {
    move |e| PtmProphetRsError::Io { path: out_path.to_path_buf(), source: e }
}

/// Crude single-attribute extractor for the streaming copy (the values we read are simple).
fn attr(line: &str, name: &str) -> Option<String> {
    let pat = format!("{name}=\"");
    let i = line.find(&pat)? + pat.len();
    let rest = &line[i..];
    let j = rest.find('"')?;
    Some(rest[..j].to_string())
}

/// Build the `<analysis_summary analysis="ptmprophet">` header block. The ROC block is
/// spliced right after `<inputfile>`, before the `<mixturemodel>` curves. `in_path` is this output's
/// own input pepXML (so `<inputfile name>` is correct per batch file).
fn build_header(config: &Config, in_path: &Path, model: &KdModel, roc_block: &str) -> String {
    let mut s = String::new();
    let infile = in_path
        .file_name()
        .map(|f| f.to_string_lossy().to_string())
        .unwrap_or_default();
    let _ = writeln!(s, "<analysis_summary analysis=\"ptmprophet\" time=\"\">");
    let _ = writeln!(
        s,
        "<ptmprophet_summary version=\"{}\" options=\"{}\" mod_string=\"{}\" min_o=\"{}\" min_o_factors=\"\" frag_ppm_tol=\"{}\">",
        SUMMARY_VERSION, config.opt_string, config.mod_string, config.min_o, config.frag_tol
    );
    let _ = writeln!(s, "<inputfile name=\"{infile}\"/>");
    s.push_str(roc_block);
    if model.ready {
        let _ = model.report(&mut s);
    }
    let _ = writeln!(s, "</ptmprophet_summary>");
    let _ = writeln!(s, "</analysis_summary>");
    s
}

/// Build the per-PSM `<analysis_result>` block. Returns `None` if no type localizes.
fn build_analysis_result(config: &Config, mpx: &Mpx) -> Option<String> {
    // Only emit if at least one type localizes: candidate residue sites, or a terminal mod.
    let any = (0..config.mods.len())
        .any(|t| !mpx.types[t].sites.is_empty() || mpx.types[t].terminus.is_some());
    if !any {
        return None;
    }
    let mut s = String::new();
    let _ = writeln!(s, "<analysis_result analysis=\"ptmprophet\">");
    let _ = writeln!(s, "<parameter name=\"mass_ppm_offset\" value=\"0\"/>");
    let _ = writeln!(s, "<parameter name=\"mass_ppm_tolerance\" value=\"{}\"/>", config.frag_tol);

    // PTMProphet builds each PSM's block in a fresh stream: the FIRST `<ptmprophet_result>` prints
    // `prior` with the default float format (`%.6g`); printing its per-site scores then switches the
    // stream to `fixed`+`setprecision(4)`, so every SUBSEQUENT result's `prior` comes out as `%.4f`
    // (e.g. first `prior="0.5"`, second `prior="0.5000"`). We replicate that sticky-stream quirk.
    let mut first_result = true;
    for (t, mt) in config.mods.iter().enumerate() {
        let tr = &mpx.types[t];
        // Terminal mod (n/c-term): unambiguous, fully localized at probability 1. Reported via
        // `<mod_terminal_probability>` with a `n(1.0000)…`/`…c(1.0000)` ptm_peptide.
        if let Some(term) = tr.terminus {
            let prior_str = if first_result { fmt_default(1.0) } else { format!("{:.4}", 1.0) };
            first_result = false;
            let seq: String = mpx_residues(mpx).iter().collect();
            let ptm_pep = if term == 'n' {
                format!("n({:.4}){}", 1.0, seq)
            } else {
                format!("{}c({:.4})", seq, 1.0)
            };
            let _ = writeln!(
                s,
                "<ptmprophet_result ptm=\"{}\" prior=\"{}\" ptm_peptide=\"{}\">",
                mt.ptm_label, prior_str, ptm_pep
            );
            let _ = writeln!(s, "<parameter name=\"mean_best_prob\" value=\"{:.4}\"/>", 1.0);
            let _ = writeln!(s, "<parameter name=\"norm_info_gain\" value=\"{:.4}\"/>", 1.0);
            let _ = writeln!(s, "<parameter name=\"localized_mods\" value=\"{:.4}\"/>", 1.0);
            let _ = writeln!(
                s,
                "<mod_terminal_probability terminus=\"{}\" probability=\"{}\" oscore=\"{:.4}\" mscore=\"{:.4}\" direct_oscore=\"{:.4}\" direct_mscore=\"{:.4}\" cterm_score=\"{:.4}\" nterm_score=\"{:.4}\" shift=\"c\"/>",
                term, fmt_sci(1.0), 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            );
            let _ = writeln!(s, "</ptmprophet_result>");
            continue;
        }
        if tr.sites.is_empty() || tr.nmods == 0 {
            continue;
        }
        let prior = tr.nmods as f64 / tr.sites.len() as f64;
        let prior_str = if first_result { fmt_default(prior) } else { format!("{:.4}", prior) };
        first_result = false;
        // ptm_peptide string: residues with (prob) at candidate sites, fixed precision 4.
        let mut ptm_pep = String::new();
        for (i, &aa) in mpx_residues(mpx).iter().enumerate() {
            ptm_pep.push(aa);
            if let Some(sc) = tr.scores.get(&i) {
                let _ = write!(ptm_pep, "({:.4})", sc.prob);
                let _ = sc; // (prob already used)
            }
        }
        let _ = writeln!(
            s,
            "<ptmprophet_result ptm=\"{}\" prior=\"{}\" ptm_peptide=\"{}\">",
            mt.ptm_label, prior_str, ptm_pep
        );
        // mean_best_prob (mean of the top-nmods probs); norm_info_gain + localized_mods.
        let mbp = mean_best_prob(tr);
        let (nig, localized) = mpx.info_gain(t);
        let _ = writeln!(s, "<parameter name=\"mean_best_prob\" value=\"{:.4}\"/>", mbp);
        let _ = writeln!(s, "<parameter name=\"norm_info_gain\" value=\"{:.4}\"/>", nig);
        let _ = writeln!(s, "<parameter name=\"localized_mods\" value=\"{:.4}\"/>", localized);
        for (&site, sc) in &tr.scores {
            let _ = writeln!(
                s,
                "<mod_aminoacid_probability position=\"{}\" probability=\"{}\" oscore=\"{:.4}\" mscore=\"{:.4}\" direct_oscore=\"{:.4}\" direct_mscore=\"{:.4}\" cterm_score=\"{:.4}\" nterm_score=\"{:.4}\" shift=\"{}\"/>",
                site + 1,
                fmt_sci(sc.prob),
                sc.oscore, sc.mscore, sc.direct_oscore, sc.direct_mscore, sc.cterm_score, sc.nterm_score, sc.shift
            );
        }
        let _ = writeln!(s, "</ptmprophet_result>");
    }
    let _ = writeln!(s, "</analysis_result>");
    Some(s)
}

/// Build the rewritten `<modification_info>` block: variable mods placed at the
/// localized sites, static mods (C/K/termini) at their fixed positions, masses recomputed from the
/// search-summary tables (full precision). `modified_peptide` uses integer brackets (mod_prec=0);
/// `mod_aminoacid_mass`/`mod_*term_mass` use fixed precision 8/7.
fn build_modification_info(config: &Config, index: &InputIndex, psm: &PsmRecord, mpx: &Mpx) -> String {
    let residues: Vec<char> = psm.pep_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let len = residues.len();
    let stat = &index.mods.stat_mods;

    // Localized variable-mod positions → shift.
    let localized = mpx.localized_positions();
    let mut var_at: HashMap<usize, f64> = HashMap::new();
    for (t, positions) in localized.iter().enumerate() {
        for &p in positions {
            var_at.insert(p, config.mods[t].shift);
        }
    }

    let nterm_mass = if psm.pep_mods.contains_key(&0) {
        Some(match stat.get(&'n') {
            Some(&s) => masses::GROUP_N_TERM + s,
            None => psm.pep_mods[&0],
        })
    } else {
        None
    };
    let cterm_mass = if psm.pep_mods.contains_key(&(len + 1)) {
        Some(match stat.get(&'c') {
            Some(&s) => masses::GROUP_C_TERM + s,
            None => psm.pep_mods[&(len + 1)],
        })
    } else {
        None
    };

    // modified_peptide string (integer brackets).
    let mut mp = String::new();
    if let Some(nm) = nterm_mass {
        let _ = write!(mp, "n[{}]", nm.round() as i64);
    }
    for (i, &aa) in residues.iter().enumerate() {
        mp.push(aa);
        let base = masses::residue_mono_mass(aa).unwrap_or(0.0);
        if let Some(&shift) = var_at.get(&i) {
            let _ = write!(mp, "[{}]", (base + shift).round() as i64);
        } else if let Some(&ss) = stat.get(&aa) {
            let _ = write!(mp, "[{}]", (base + ss).round() as i64);
        }
    }
    if let Some(cm) = cterm_mass {
        let _ = write!(mp, "c[{}]", cm.round() as i64);
    }

    let mut s = String::new();
    let _ = write!(s, "<modification_info modified_peptide=\"{mp}\"");
    if let Some(nm) = nterm_mass {
        let _ = write!(s, " mod_nterm_mass=\"{:.7}\"", nm);
    }
    if let Some(cm) = cterm_mass {
        let _ = write!(s, " mod_cterm_mass=\"{:.7}\"", cm);
    }
    let _ = writeln!(s, ">");
    for (i, &aa) in residues.iter().enumerate() {
        let base = masses::residue_mono_mass(aa).unwrap_or(0.0);
        if let Some(&shift) = var_at.get(&i) {
            let _ = writeln!(
                s,
                "<mod_aminoacid_mass position=\"{}\" mass=\"{:.8}\" variable=\"{:.8}\" source=\"param\"/>",
                i + 1, base + shift, shift
            );
        } else if let Some(&ss) = stat.get(&aa) {
            let _ = writeln!(
                s,
                "<mod_aminoacid_mass position=\"{}\" mass=\"{:.8}\" static=\"{:.8}\" variable=\"0.00000000\" source=\"param\"/>",
                i + 1, base + ss, ss
            );
        }
    }
    let _ = writeln!(s, "</modification_info>");
    s
}

/// Reformat a **non-localized** PSM's original `<modification_info>` into PTMProphet's canonical
/// "pass-through" form. Unlike [`build_modification_info`] (which rebuilds a *localized* PSM from its
/// `Mpx` at full `%.8f` precision), this works straight from the input sites and reproduces the
/// distinct precision PTMProphet uses on the peptides it does **not** score: its output stream is in
/// its default float mode, so every mass is `%.8g` (8 significant figures), not fixed `%.8f`.
///
/// Per static site PTMProphet emits `mass`, `static` and `variable` where, quirkily, `mass` is the
/// recomputed residue (full precision) + static while `variable` is `input_mass − static − residue`
/// using the *parameter-rounded* residue mass — a near-zero float residue (e.g. `1.4210855e-14` for
/// C+57). We reproduce that residue exactly: the parameter residue is `round5(input_mass − static)`
/// (the input mass already encodes `residue₅ + static`), so `variable = (input_mass − static) −
/// round5(input_mass − static)`. `modified_peptide` carries integer mod brackets (`mod_prec = 0`).
fn build_passthrough_modinfo(
    pep_seq: &str,
    masses: &[(usize, f64)],
    nterm_mass: Option<f64>,
    cterm_mass: Option<f64>,
    stat_mods: &HashMap<char, f64>,
) -> String {
    let residues: Vec<char> = pep_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let mass_at: HashMap<usize, f64> = masses.iter().copied().collect();

    // modified_peptide string (integer brackets at every modified site).
    let mut mp = String::new();
    if let Some(nm) = nterm_mass {
        let _ = write!(mp, "n[{}]", nm.round() as i64);
    }
    for (i, &aa) in residues.iter().enumerate() {
        mp.push(aa);
        if let Some(&m) = mass_at.get(&(i + 1)) {
            let _ = write!(mp, "[{}]", m.round() as i64);
        }
    }
    if let Some(cm) = cterm_mass {
        let _ = write!(mp, "c[{}]", cm.round() as i64);
    }

    let mut s = String::new();
    let _ = write!(s, "<modification_info modified_peptide=\"{mp}\"");
    if let Some(nm) = nterm_mass {
        let _ = write!(s, " mod_nterm_mass=\"{:.7}\"", nm);
    }
    if let Some(cm) = cterm_mass {
        let _ = write!(s, " mod_cterm_mass=\"{:.7}\"", cm);
    }
    let _ = writeln!(s, ">");
    for &(pos, m_in) in masses {
        let aa = residues.get(pos.wrapping_sub(1)).copied().unwrap_or('X');
        let base = masses::residue_mono_mass(aa).unwrap_or(0.0);
        if let Some(&stat) = stat_mods.get(&aa) {
            // Static site: residue + static, with the near-zero float residue PTMProphet leaves.
            let diff = m_in - stat;
            let std_aa = (diff * 1e5).round() / 1e5;
            let variable = diff - std_aa;
            let _ = writeln!(
                s,
                "<mod_aminoacid_mass position=\"{}\" mass=\"{}\" static=\"{}\" variable=\"{}\" source=\"param\"/>",
                pos,
                fmt_g8(base + stat),
                fmt_g8(stat),
                fmt_g8(variable)
            );
        } else {
            // Variable-only site (not seen in the supported configs' pass-through, kept for safety):
            // mass is the input token mass; variable is mass − bare residue.
            let _ = writeln!(
                s,
                "<mod_aminoacid_mass position=\"{}\" mass=\"{}\" variable=\"{}\" source=\"param\"/>",
                pos,
                fmt_g8(m_in),
                fmt_g8(m_in - base)
            );
        }
    }
    let _ = writeln!(s, "</modification_info>");
    s
}

/// C++ `ostream <<` with `precision(8)` and no `fixed` flag — i.e. `%.8g`. PTMProphet's pass-through
/// `<mod_aminoacid_mass>` masses are written in this default stream mode.
fn fmt_g8(x: f64) -> String {
    crate::fmt::fmt_g(x, 8)
}

fn mean_best_prob(tr: &crate::mpx::TypeResult) -> f64 {
    let mut probs: Vec<f64> = tr.scores.values().map(|s| s.prob).collect();
    probs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
    let k = tr.nmods.min(probs.len()).max(1);
    probs.iter().take(k).sum::<f64>() / k as f64
}

fn mpx_residues<'a>(mpx: &'a Mpx) -> &'a [char] {
    mpx.residues_view()
}

/// Default-ostream formatting (C++ `<<`, i.e. `%.6g`) — used for the first `prior` in each PSM block.
/// See [`crate::fmt::fmt_g`] for the shared C-printf reproduction.
fn fmt_default(x: f64) -> String {
    crate::fmt::fmt_g6(x)
}

/// Scientific with 4 decimals + 2-digit exponent (e.g. `3.0315e-01`).
fn fmt_sci(x: f64) -> String {
    let s = format!("{:.4e}", x); // e.g. "3.0315e-1"
    // pad exponent to 2 digits with sign
    if let Some(pos) = s.find('e') {
        let (mant, exp) = s.split_at(pos);
        let exp = &exp[1..]; // skip 'e'
        let (sign, digits) = if let Some(d) = exp.strip_prefix('-') {
            ("-", d)
        } else {
            ("+", exp.trim_start_matches('+'))
        };
        format!("{mant}e{sign}{:0>2}", digits)
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::Config;
    use crate::em::kdmodel::KdModel;
    use crate::model::masses::IonKind;
    use std::path::PathBuf;

    fn stat_c() -> HashMap<char, f64> {
        let mut m = HashMap::new();
        m.insert('C', 57.02146);
        m
    }

    /// `%.8g` (C++ `precision(8)`, no `fixed`) for the three masses on a pass-through static site.
    #[test]
    fn fmt_g8_matches_cpp_precision8() {
        assert_eq!(fmt_g8(160.0306445), "160.03064");
        assert_eq!(fmt_g8(57.02146), "57.02146");
        // The near-zero float residue PTMProphet leaves on a static C+57 site.
        assert_eq!(fmt_g8(1.4210854715202004e-14), "1.4210855e-14");
    }

    /// An n-terminal mod (acetyl) is localized: the `<analysis_result>` carries a
    /// `<ptmprophet_result ptm="n:…">` with `n(1.0000)…` and a `<mod_terminal_probability>` — exactly
    /// like PTMProphet, which reports terminal mods at probability 1.
    #[test]
    fn nterm_mod_localizes_and_emits_terminal_probability() {
        use crate::model::modspec::ModType;
        use crate::pepxml::index::{InputIndex, ModTables};

        let mods = vec![
            ModType { aas: vec!['M'], shift: 15.9949, ptm_label: "M:15.9949".into() },
            ModType { aas: vec!['n'], shift: 42.0106, ptm_label: "n:42.0106".into() },
        ];
        let mut tables = ModTables::default();
        tables.stat_mods.insert('C', 57.02146);

        // n-term acetyl (43.0184250) + static C at residue 12; no M present.
        let mut pep_mods = std::collections::BTreeMap::new();
        pep_mods.insert(0usize, 43.0184250);
        pep_mods.insert(12usize, 160.03064450);
        let mpx = Mpx::build(&min_config(), &mods, &tables, "AADISESSGADCK", &pep_mods, 2, 0.0);

        // The n-term type is flagged terminal (one unambiguous site) → there is work to do.
        assert_eq!(mpx.types[1].terminus, Some('n'));
        assert!(mpx.has_work());

        let mut config = min_config();
        config.mods = mods.clone();
        config.frag_tol = 10.0;
        let mut index = InputIndex::default();
        index.mods = tables;
        let psm = PsmRecord {
            pep_seq: "AADISESSGADCK".into(),
            pep_mods,
            ..Default::default()
        };
        let (block, _rocs) = build_psm_block(&config, &index, &psm, &mpx).expect("localizes");
        assert!(block.contains("<ptmprophet_result ptm=\"n:42.0106\" prior=\"1\" ptm_peptide=\"n(1.0000)AADISESSGADCK\">"), "block:\n{block}");
        assert!(block.contains("<mod_terminal_probability terminus=\"n\" probability=\"1.0000e+00\" oscore=\"1.0000\" mscore=\"1.0000\" direct_oscore=\"1.0000\" direct_mscore=\"1.0000\" cterm_score=\"1.0000\" nterm_score=\"1.0000\" shift=\"c\"/>"), "block:\n{block}");
        // Static C still rendered in the localized %.8f form.
        assert!(block.contains("<mod_aminoacid_mass position=\"12\" mass=\"160.03064450\" static=\"57.02146000\" variable=\"0.00000000\" source=\"param\"/>"), "block:\n{block}");
    }

    /// A static-only PSM is reformatted byte-for-byte like PTMProphet's pass-through output
    /// (`C[160]C[160]TESLVNR`, position-first sites with `static`/`variable`/`source`).
    #[test]
    fn passthrough_static_c_matches_ptmprophet() {
        let got = build_passthrough_modinfo(
            "CCTESLVNR",
            &[(1, 160.03065), (2, 160.03065)],
            None,
            None,
            &stat_c(),
        );
        let want = "<modification_info modified_peptide=\"C[160]C[160]TESLVNR\">\n\
<mod_aminoacid_mass position=\"1\" mass=\"160.03064\" static=\"57.02146\" variable=\"1.4210855e-14\" source=\"param\"/>\n\
<mod_aminoacid_mass position=\"2\" mass=\"160.03064\" static=\"57.02146\" variable=\"1.4210855e-14\" source=\"param\"/>\n\
</modification_info>\n";
        assert_eq!(got, want);
    }

    fn min_config() -> Config {
        Config {
            mods: Vec::new(),
            mod_string: String::new(),
            opt_string: String::new(),
            frag_tol: 15.0,
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

    /// End-to-end stream: an unmodified hit gets an empty `<modification_info>`, a static-only hit is
    /// reformatted, and the `<analysis_summary>` header is spliced in — matching PTMProphet's
    /// "every hit carries a modification_info" invariant.
    #[test]
    fn stream_output_rewrites_every_hit() {
        let content = "<msms_pipeline_analysis>\n\
<spectrum_query start_scan=\"1\">\n\
<search_result>\n\
<search_hit peptide=\"QLVQDENVR\" hit_rank=\"1\">\n\
<search_score name=\"hyperscore\" value=\"1\"/>\n\
</search_hit>\n\
</search_result>\n\
</spectrum_query>\n\
<spectrum_query start_scan=\"2\">\n\
<search_result>\n\
<search_hit peptide=\"CCTESLVNR\" hit_rank=\"1\">\n\
<modification_info modified_peptide=\"CCTESLVNR\">\n\
<mod_aminoacid_mass mass=\"160.03065\" position=\"1\"/>\n\
<mod_aminoacid_mass mass=\"160.03065\" position=\"2\"/>\n\
</modification_info>\n\
<search_score name=\"hyperscore\" value=\"2\"/>\n\
</search_hit>\n\
</search_result>\n\
</spectrum_query>\n\
</msms_pipeline_analysis>\n";

        let dir = std::env::temp_dir();
        let in_path = dir.join("ptmprophet_rs_stream_in.pep.xml");
        let out_path = dir.join("ptmprophet_rs_stream_out.mod.pep.xml");
        let model = KdModel::new("PTMProphet_Oscore");
        let blocks: HashMap<i64, String> = HashMap::new();

        stream_output(&min_config(), &stat_c(), &in_path, &out_path, content, &model, "", &blocks)
            .unwrap();
        let got = std::fs::read_to_string(&out_path).unwrap();
        std::fs::remove_file(&out_path).ok();

        // Unmodified peptide: an empty modification_info was inserted before its search_score.
        assert!(
            got.contains("<modification_info modified_peptide=\"QLVQDENVR\">\n</modification_info>\n<search_score name=\"hyperscore\" value=\"1\"/>"),
            "missing empty modinfo for unmodified peptide:\n{got}"
        );
        // Static-only peptide: reformatted to the canonical pass-through form.
        assert!(
            got.contains("<modification_info modified_peptide=\"C[160]C[160]TESLVNR\">\n<mod_aminoacid_mass position=\"1\" mass=\"160.03064\" static=\"57.02146\" variable=\"1.4210855e-14\" source=\"param\"/>"),
            "static-only modinfo not reformatted:\n{got}"
        );
        // Original mass-first lines are gone; the header was spliced in.
        assert!(!got.contains("mass=\"160.03065\""), "original input mass leaked through:\n{got}");
        assert!(got.contains("<analysis_summary analysis=\"ptmprophet\""), "header not spliced:\n{got}");
        // The summary names PTMProphet-rs and its version; the element names stay PTMProphet's.
        let summary = format!(
            "<ptmprophet_summary version=\"PTMProphet-rs {} (Rust port of TPP PTMProphet)\" options=",
            env!("CARGO_PKG_VERSION")
        );
        assert!(got.contains(&summary), "summary version not PTMProphet-rs's:\n{got}");
    }
}
