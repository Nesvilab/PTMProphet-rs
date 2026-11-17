// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// SPDX-FileCopyrightText: 2007 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetParser.cpp and
// src/Validation/InterProphet/InterProphetParser/KDModel.cpp (LGPL-2.1-or-later); modified by the
// University of Michigan, 2026.

//! EM driver — fits the intensity model and iterates the E-step/M-step.
//!
//! Inputs are the per-site arrays harvested in the read pass: `priors`
//! (= `0.999999·Oscore`) and `ovals` (= Oscore). For EM=1 the inserted weight and the per-site
//! prior are both `0.999999·Oscore`.
//!
//! Loop, ≤ 99 iterations:
//!   E-step: `prob_i = π_i·f_pos(o_i) / (π_i·f_pos(o_i) + (1−π_i)·f_neg(o_i))`, `π_i = priors[i]`
//!           (the per-site prior, every iteration — the global-prior path is unused for EM=1).
//!   M-step: `model.update(prob_i, i)` for all i; `update_prob_sums(Σprob)`; `make_ready()` (re-fit).
//!   Converge when `|Σprob − lastΣprob| < 0.001`.

use crate::config::Config;
use crate::em::kdmodel::KdModel;
use crate::error::Result;

const MAX_EM_ITERS: usize = 99; // while iter_ < 100
const CONVERGENCE: f64 = 0.001; // |Δprobsum| < 0.001
/// EM prior/weight cap applied to every Oscore before insertion.
pub const PRIOR_CAP: f64 = 0.999999;

/// Per-site evidence harvested in the read pass and consumed by the EM.
/// Parallel vectors, one entry per ambiguous `(PSM, type, site)` triple, in parse order.
#[derive(Debug, Default, Clone)]
pub struct EmEvidence {
    /// `0.999999·Oscore` — both the inserted weight and the per-site EM prior π_i.
    pub priors: Vec<f64>,
    /// The Oscore value binned into the model (the density abscissa).
    pub ovals: Vec<f64>,
    pub total_mods: usize,
    pub total_sites: usize,
}

impl EmEvidence {
    /// Harvest one ambiguous site: `prior = 0.999999·oscore`, value = `oscore`.
    pub fn push(&mut self, oscore: f64) {
        self.priors.push(PRIOR_CAP * oscore);
        self.ovals.push(oscore);
    }
    pub fn len(&self) -> usize {
        self.ovals.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ovals.is_empty()
    }
}

/// Fit the intensity KDModel and EM-iterate the priors to convergence.
pub fn compute_ptm_model(_config: &Config, model: &mut KdModel, evidence: &EmEvidence) -> Result<()> {
    if evidence.is_empty() {
        return Ok(());
    }
    // Initial fit: insert (prior, oscore) per site, then build densities.
    for i in 0..evidence.len() {
        model.insert(evidence.priors[i], evidence.ovals[i]);
    }
    model.make_ready();

    let n = model.vals.len();
    let mut probarr = vec![0.0_f64; n];
    let mut last_probsum = 0.0_f64;

    for iter in 1..=MAX_EM_ITERS {
        // E-step.
        let mut probsum = 0.0;
        for (i, &o) in model.vals.iter().enumerate() {
            let pos = model.get_pos_prob(o);
            let neg = model.get_neg_prob(o);
            let pi = evidence.priors[i]; // per-site prior, every iteration
            let denom = pi * pos + (1.0 - pi) * neg;
            let prob = if denom > 0.0 { pi * pos / denom } else { pi };
            probarr[i] = prob;
            probsum += prob;
        }

        // Convergence: stop BEFORE the M-step of the converging iteration.
        if iter > 1 && (probsum - last_probsum).abs() < CONVERGENCE {
            break;
        }

        // M-step: re-weight the model with the new posteriors and re-fit.
        for (i, &prob) in probarr.iter().enumerate() {
            model.update(prob, i);
        }
        model.update_prob_sums(probsum);
        model.make_ready();

        last_probsum = probsum;
    }
    // Fit the observed densities from the converged posteriors (the `<mixturemodel>` obs columns).
    model.compute_obs_fits();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::masses::IonKind;
    use std::path::PathBuf;

    fn cfg() -> Config {
        Config {
            mods: Vec::new(), mod_string: String::new(), opt_string: String::new(),
            frag_tol: 10.0, em: 1, min_prob: 0.5,
            mod_prec: 0, min_o: 0, max_frag_z: 0, threads: None, nostack: true,
            static_tol: true,
            nions: vec![IonKind::B], cions: vec![IonKind::Y],
            inputs: Vec::new(), in_file: PathBuf::new(), out_file: PathBuf::new(),
        }
    }

    #[test]
    fn em_converges_and_separates() {
        // Bimodal oscores: many near 0.8 (true sites) and many near 0.2 (decoy sites).
        let mut ev = EmEvidence::default();
        for k in 0..200 {
            let hi = 0.75 + (k % 10) as f64 * 0.005;
            let lo = 0.15 + (k % 10) as f64 * 0.005;
            ev.push(hi);
            ev.push(lo);
        }
        let mut model = KdModel::new("PTMProphet_Oscore");
        compute_ptm_model(&cfg(), &mut model, &ev).unwrap();
        assert!(model.ready);
        // A high oscore should map to a high posterior, a low oscore to a low posterior.
        let post = |o: f64| {
            let pos = model.get_pos_prob(o);
            let neg = model.get_neg_prob(o);
            let pi = PRIOR_CAP * o;
            pi * pos / (pi * pos + (1.0 - pi) * neg)
        };
        assert!(post(0.8) > 0.6, "post(0.8)={}", post(0.8));
        assert!(post(0.2) < 0.4, "post(0.2)={}", post(0.2));
    }
}
