// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp and
// src/Validation/PTMProphetParser/PTMProphetParser.h (LGPL-2.1-or-later); modified by the
// University of Michigan, 2026.

//! ROC / error-rate summary tables.
//!
//! After the localizing write pass, the ROC tables are spliced into the header
//! `<ptmprophet_summary>` right after `<inputfile>`, emitted inline. All ROC quantities are
//! **PSM-probability-weighted** via [`InfoLocs`].

/// One localization record fed into the ROC tables. Stored values are weighted by the PSM
/// probability: `locs = locs·prob`, `num_mods = mods·prob`, `norm_info = info·prob`
/// (when `prob >= 0`). `mean_best_prob` and `prob` are stored unweighted.
#[derive(Debug, Clone, Copy)]
pub struct InfoLocs {
    pub norm_info: f32,
    pub mean_best_prob: f32,
    pub locs: f32,
    pub prob: f32,
    pub num_mods: f32,
}

impl InfoLocs {
    /// Build a record, applying the PSM probability weighting.
    pub fn new(locs: f32, info: f32, mean_best_prob: f32, prob: f32, mods: f32) -> InfoLocs {
        let (locs, num_mods, norm_info) = if prob >= 0.0 {
            (locs * prob, mods * prob, info * prob)
        } else {
            (locs, mods, info)
        };
        InfoLocs { norm_info, mean_best_prob, locs, prob, num_mods }
    }
}

use std::fmt::Write as _;

const THRESH: [f64; 24] = [
    0.9999, 0.999, 0.99, 0.98, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4,
    0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.0,
];
const ERROR_RATES: [f64; 31] = [
    0.0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009, 0.001, 0.0015,
    0.002, 0.0025, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03,
    0.04, 0.05, 0.075, 0.1, 0.15,
];

/// Build the full `<roc_error_data>` block for all mod types.
/// One `<roc_error_data charge="all" ptm="...">` per type, containing `<roc_data_point>` rows
/// (sensitivity/error walk down the probability thresholds) then `<error_point>` rows.
pub fn build_roc_block(by_type: &[Vec<InfoLocs>], ptm_labels: &[String]) -> String {
    let mut out = String::new();
    for (ty, entries) in by_type.iter().enumerate() {
        // Sort: descending mean_best_prob, ties descending locs.
        let mut e = entries.clone();
        e.sort_by(|a, b| {
            b.mean_best_prob
                .partial_cmp(&a.mean_best_prob)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(b.locs.partial_cmp(&a.locs).unwrap_or(std::cmp::Ordering::Equal))
        });
        let totcorrect: f64 = e.iter().map(|x| x.locs as f64).sum();
        let label = ptm_labels.get(ty).map(|s| s.as_str()).unwrap_or("");
        let _ = writeln!(out, "<roc_error_data charge=\"all\" ptm=\"{label}\">");

        // ---- roc_data_point rows ----
        if totcorrect > 0.0 {
            let mut correct = 0.0;
            let mut incorrect = 0.0;
            let mut ti = 0usize;
            let mut k = 0usize;
            while k < e.len() && ti < THRESH.len() {
                let rnd: f64 = format!("{:.4}", e[k].mean_best_prob as f64).parse().unwrap_or(0.0);
                if rnd >= THRESH[ti] {
                    correct += e[k].locs as f64;
                    incorrect += e[k].num_mods as f64 - e[k].locs as f64;
                    k += 1;
                } else {
                    let sens = (correct / totcorrect).max(0.0);
                    let err = if correct + incorrect > 0.0 {
                        format!("{:.4}", (incorrect / (correct + incorrect)).max(0.0))
                    } else {
                        format!("{:.0}", 0.0)
                    };
                    let _ = writeln!(
                        out,
                        "<roc_data_point min_prob=\"{:.4}\" sensitivity=\"{:.4}\" error=\"{}\" num_corr=\"{:.0}\" num_incorr=\"{:.0}\"/>",
                        THRESH[ti], sens, err, correct, incorrect
                    );
                    ti += 1; // re-test the same entry against the next threshold
                }
            }
            // final point for the threshold not yet crossed
            if ti < THRESH.len() {
                let _ = writeln!(
                    out,
                    "<roc_data_point min_prob=\"{:.4}\" sensitivity=\"{:.4}\" error=\"{:.4}\" num_corr=\"{:.0}\" num_incorr=\"{:.0}\"/>",
                    THRESH[ti], correct / totcorrect, incorrect / (correct + incorrect), correct, incorrect
                );
            }
        } else {
            let _ = writeln!(
                out,
                "<roc_data_point min_prob=\"0.0000\" sensitivity=\"0.0000\" error=\"0.0000\" num_corr=\"0.0\" num_incorr=\"0.0\"/>"
            );
        }

        // ---- error_point rows ----
        if totcorrect > 0.0 {
            let mut correct = 0.0;
            let mut incorrect = 0.0;
            let mut ti = 0usize;
            let mut k = 0usize;
            while k < e.len() {
                let ratio = if correct + incorrect > 0.0 {
                    incorrect / (correct + incorrect)
                } else {
                    f64::NAN // 0/0 → NaN; NaN > x is false → accumulate
                };
                if ratio > ERROR_RATES[ti] {
                    let idx = if k == 0 { 0 } else { k - 1 };
                    let _ = writeln!(
                        out,
                        "<error_point error=\"{:.4}\" min_prob=\"{:.4}\" num_corr=\"{:.0}\" num_incorr=\"{:.0}\"/>",
                        ERROR_RATES[ti],
                        e[idx].mean_best_prob as f64,
                        correct + e[idx].locs as f64,
                        incorrect + e[idx].num_mods as f64 - e[idx].locs as f64
                    );
                    ti += 1; // re-test same entry
                } else {
                    correct += e[k].locs as f64;
                    incorrect += e[k].num_mods as f64 - e[k].locs as f64;
                    k += 1;
                }
                if ti == ERROR_RATES.len() - 1 {
                    break; // done
                }
            }
        }

        let _ = writeln!(out, "</roc_error_data>");
    }
    out
}
