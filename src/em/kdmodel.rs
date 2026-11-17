// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2007 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/InterProphet/InterProphetParser/KDModel.cpp
// and src/Validation/InterProphet/InterProphetParser/KDModel.h (LGPL-2.1-or-later); modified by the
// University of Michigan, 2026.

//! `KDModel` — the 1-D mixture-model container. Holds the inserted `(weight, value)`
//! pairs and, after [`KdModel::make_ready`], the fitted positive/negative densities over a 500-point
//! uniform grid spanning `[min−0.01, max+0.01]` with fixed Gaussian bandwidth `(range+0.02)/5`
//! (breaks=5). Density queries return **1** when the model is not ready (degenerate data) so the EM
//! posterior passes the prior through unchanged.
//!
//! Insertion stores `weight` and `1−weight` per point (parallel arrays addressed by index during
//! the M-step `update`).

use crate::em::density;

const GRIDS: usize = 500; // number of uniform grid points
const BREAKS: f64 = 5.0; // bandwidth divisor: bw = range / breaks
const PAD: f64 = 0.01; // grid padding: minval -= 0.01; maxval += 0.01
const OBS_HIGH_PR: f64 = 0.80; // posterior ≥ this → "observed positive" value
const OBS_LOW_PR: f64 = 0.20; // posterior ≤ this → "observed negative" value

/// 1-D non-parametric mixture model.
#[derive(Debug, Clone, Default)]
pub struct KdModel {
    pub name: String,
    /// Parallel insertion-ordered arrays.
    pub vals: Vec<f64>,
    pub posprobs: Vec<f64>,
    pub negprobs: Vec<f64>,
    /// Fitted state (set by `make_ready`).
    pub grid: Vec<f64>,
    pub pos_dens: Vec<f64>,
    pub neg_dens: Vec<f64>,
    /// Observed densities — unweighted KDE of the high- (`pos_obs_dens`) and low-posterior
    /// (`neg_obs_dens`) values. Computed once at the first `make_ready` and frozen thereafter.
    pub pos_obs_dens: Vec<f64>,
    pub neg_obs_dens: Vec<f64>,
    pub pos_bw: f64,
    pub neg_bw: f64,
    pub postot: f64,
    pub negtot: f64,
    pub min_val: f64,
    pub max_val: f64,
    pub ready: bool,
}

impl KdModel {
    pub fn new(name: &str) -> KdModel {
        KdModel { name: name.to_string(), min_val: 0.0, max_val: 0.0, ..Default::default() }
    }

    /// `insert(prob, val)`: store `prob` and `1−prob` and `val` in parallel; track
    /// min/max and running totals; silently drop `prob ∉ [0,1]`.
    pub fn insert(&mut self, prob: f64, val: f64) {
        if !(0.0..=1.0).contains(&prob) {
            return;
        }
        if val > self.max_val {
            self.max_val = val;
        }
        if val < self.min_val {
            self.min_val = val;
        }
        self.postot += prob;
        self.negtot += 1.0 - prob;
        self.posprobs.push(prob);
        self.negprobs.push(1.0 - prob);
        self.vals.push(val);
    }

    /// Build the grid (once) and fit the pos/neg densities from current weights.
    /// Sets `ready=false` on degenerate data (min==max), so queries
    /// then return 1.
    pub fn make_ready(&mut self) {
        if self.grid.is_empty() {
            // Build grid from data min/max: pad, derive bandwidth, lay down 500 uniform points.
            if self.vals.is_empty() {
                self.ready = false;
                return;
            }
            let mut lo = self.vals[0];
            let mut hi = self.vals[0];
            for &v in &self.vals {
                if v < lo {
                    lo = v;
                }
                if v > hi {
                    hi = v;
                }
            }
            if lo == hi {
                // degenerate: cannot classify → not ready.
                self.ready = false;
                return;
            }
            lo -= PAD;
            hi += PAD;
            self.pos_bw = (hi - lo) / BREAKS;
            self.neg_bw = (hi - lo) / BREAKS;
            self.grid = (0..GRIDS)
                .map(|i| lo + (i as f64) * (hi - lo) / (GRIDS as f64 - 1.0))
                .collect();
        }
        self.pos_dens = density::wt_density_fit(&self.grid, &self.vals, self.pos_bw, &self.posprobs);
        self.neg_dens = density::wt_density_fit(&self.grid, &self.vals, self.neg_bw, &self.negprobs);
        self.ready = true;
    }

    /// Fit the observed densities from the CURRENT posteriors: an **unweighted** KDE of the high-
    /// (`posprob ≥ 0.80`, → `pos_obs_dens`) and low- (`posprob ≤ 0.20`, → `neg_obs_dens`) posterior
    /// values, both using `pos_bw`. PTMProphet emits these as `pos_obs_dens`/`neg_obs_dens` in the
    /// `<mixturemodel>` curve. Call once after the EM converges so the split reflects the final
    /// posteriors. Unweighted fit = weighted fit with all weights 1 (`ΣwK/Σw = ΣK/N`).
    pub fn compute_obs_fits(&mut self) {
        if self.grid.is_empty() {
            return;
        }
        let mut highvals = Vec::new();
        let mut lowvals = Vec::new();
        for (i, &v) in self.vals.iter().enumerate() {
            if self.posprobs[i] >= OBS_HIGH_PR {
                highvals.push(v);
            }
            if self.posprobs[i] <= OBS_LOW_PR {
                lowvals.push(v);
            }
        }
        let hi_w = vec![1.0; highvals.len()];
        let lo_w = vec![1.0; lowvals.len()];
        self.pos_obs_dens = density::wt_density_fit(&self.grid, &highvals, self.pos_bw, &hi_w);
        self.neg_obs_dens = density::wt_density_fit(&self.grid, &lowvals, self.pos_bw, &lo_w);
    }

    /// M-step weight update for point `i`: adjust totals by the
    /// delta and set `posprobs[i]=prob`, `negprobs[i]=1−prob`. The `val` is ignored (fit points fixed).
    pub fn update(&mut self, prob: f64, i: usize) {
        if i >= self.posprobs.len() {
            return;
        }
        let old = self.posprobs[i];
        self.postot += prob - old;
        self.negtot += (1.0 - prob) - (1.0 - old);
        self.posprobs[i] = prob;
        self.negprobs[i] = 1.0 - prob;
    }

    /// Hard-reset the running totals.
    pub fn update_prob_sums(&mut self, probsum: f64) {
        self.postot = probsum;
        self.negtot = self.vals.len() as f64 - probsum;
    }

    /// Positive density at `val`. Returns 1 if `!ready`.
    pub fn get_pos_prob(&self, val: f64) -> f64 {
        if !self.ready {
            return 1.0;
        }
        density::predict_density_fit(&self.grid, &self.pos_dens, val)
    }

    /// Negative density at `val`. Returns 1 if `!ready`.
    pub fn get_neg_prob(&self, val: f64) -> f64 {
        if !self.ready {
            return 1.0;
        }
        density::predict_density_fit(&self.grid, &self.neg_dens, val)
    }

    /// Serialize the `<mixturemodel>` report block for the output header: the full grid, every number
    /// formatted with C++ `ostream <<` default precision (`%.6g`). Each `<point>` also carries the
    /// observed-density columns (`neg_obs_dens`/`pos_obs_dens`) when they were fitted, matching
    /// PTMProphet's `KDModel::report`. Returns nothing extra if not ready.
    pub fn report(&self, out: &mut dyn std::fmt::Write) -> std::fmt::Result {
        use crate::fmt::fmt_g6;
        writeln!(
            out,
            "<mixturemodel name=\"{}\" pos_bandwidth=\"{}\" neg_bandwidth=\"{}\">",
            self.name,
            fmt_g6(self.pos_bw),
            fmt_g6(self.neg_bw)
        )?;
        if self.ready {
            let has_obs = self.pos_obs_dens.len() == self.grid.len()
                && self.neg_obs_dens.len() == self.grid.len();
            for i in 0..self.grid.len() {
                write!(
                    out,
                    "<point value=\"{}\" pos_dens=\"{}\" neg_dens=\"{}\"",
                    fmt_g6(self.grid[i]),
                    fmt_g6(self.pos_dens[i]),
                    fmt_g6(self.neg_dens[i])
                )?;
                if has_obs {
                    write!(
                        out,
                        " neg_obs_dens=\"{}\" pos_obs_dens=\"{}\"",
                        fmt_g6(self.neg_obs_dens[i]),
                        fmt_g6(self.pos_obs_dens[i])
                    )?;
                }
                writeln!(out, "/>")?;
            }
        }
        writeln!(out, "</mixturemodel>")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ready_and_separation() {
        // Build a separable model: high oscores are "positive", low are "negative".
        let mut m = KdModel::new("test");
        for _ in 0..50 {
            m.insert(0.9, 0.9); // weight 0.9 at value 0.9
            m.insert(0.1, 0.1); // weight 0.1 at value 0.1
        }
        m.make_ready();
        assert!(m.ready);
        assert_eq!(m.grid.len(), 500);
        // pos density should favor high values; neg density should favor low values.
        assert!(m.get_pos_prob(0.9) > m.get_neg_prob(0.9));
        assert!(m.get_neg_prob(0.1) > m.get_pos_prob(0.1));
    }

    #[test]
    fn degenerate_not_ready_returns_one() {
        let mut m = KdModel::new("deg");
        for _ in 0..10 {
            m.insert(0.5, 0.5); // all identical values → min==max
        }
        m.make_ready();
        assert!(!m.ready);
        assert_eq!(m.get_pos_prob(0.5), 1.0);
        assert_eq!(m.get_neg_prob(0.5), 1.0);
    }
}
