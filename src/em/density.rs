// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2007 David Shteynberg
// SPDX-FileCopyrightText: 2003 Andrew Keller
// SPDX-FileCopyrightText: 2007 Henry Lam
// Translated into Rust from TPP 6.3.2 src/Validation/Distribution/NonParametricDistribution.cpp,
// src/Validation/Distribution/NonParametricDistribution.h,
// src/Validation/Distribution/GaussianDistribution.cpp and
// src/Validation/InterProphet/InterProphetParser/KDModel.cpp (LGPL-2.1-or-later); modified by the
// University of Michigan, 2026.

//! Weighted Gaussian kernel-density estimate + grid interpolation, providing
//! the non-parametric density model used by the EM fit.
//!
//! Density at grid point g:
//! ```text
//!            Σ_j  w_j · (1/bw)·φ((x_j − g)/bw)
//! f(g)  =  ──────────────────────────────────────      φ(x) = exp(−0.5 x²) / sqrt(6.28318)
//!                       Σ_j  w_j
//! ```
//! The normalizer uses the literal truncated 2π constant `6.28318` — kept exact for reproducibility.
//! Summation order: outer over the grid, inner over data points in insertion order;
//! one `f64` accumulation per grid point, then one divide by `Σw`.

/// Standard-normal pdf using the literal `sqrt(6.28318)` denominator.
///
/// The `6.28318` is a deliberately truncated 2π — it reproduces PTMProphet's C++ `getGaussianProb`
/// bit-for-bit, so it must NOT be replaced with `std::f64::consts::TAU`.
#[inline]
#[allow(clippy::approx_constant)]
pub fn gaussian_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (6.28318f64).sqrt()
}

/// Weighted KDE over `grid`, evaluated at every grid point. `weights[j]` is the
/// effective weight of data point `vals[j]` (pass `posprobs` for the positive density, `negprobs`
/// for the negative). Returns a length-`grid.len()` density vector.
///
/// This is the EM hot loop: `make_ready` calls it twice per iteration over the full evidence set
/// (O(|grid|·|vals|) Gaussian evaluations). The outer loop over grid points is independent, and each
/// grid point's inner accumulation keeps the fixed data (insertion) order, so fanning grid points out
/// across threads is **bit-identical** to the serial fit while using all cores.
pub fn wt_density_fit(grid: &[f64], vals: &[f64], bw: f64, weights: &[f64]) -> Vec<f64> {
    let wtsum: f64 = weights.iter().sum();
    let inv_bw = 1.0 / bw;
    // Density at one grid point: Σ_j w_j·(1/bw)·φ((x_j − g)/bw) / Σw, summed in data order.
    let at = |g: f64| -> f64 {
        let mut d = 0.0;
        for (j, &v) in vals.iter().enumerate() {
            d += weights[j] * inv_bw * gaussian_pdf((v - g) * inv_bw);
        }
        if wtsum > 0.0 {
            d / wtsum
        } else {
            d
        }
    };
    // Parallelize across grid points only when the work is large enough to amortize thread spawn.
    if grid.len().saturating_mul(vals.len()) >= 50_000 {
        crate::threads::par_map(grid, |&g| at(g))
    } else {
        grid.iter().map(|&g| at(g)).collect()
    }
}

/// Linear interpolation of `dens` (aligned to the uniform `grid`) at `val`.
/// Returns 0 outside `[grid[0], grid[last]]`; `val == grid node` returns the
/// node value exactly.
pub fn predict_density_fit(grid: &[f64], dens: &[f64], val: f64) -> f64 {
    let n = grid.len();
    if n == 0 || dens.len() != n {
        return 0.0;
    }
    if val < grid[0] || val > grid[n - 1] {
        return 0.0;
    }
    if n == 1 {
        return dens[0];
    }
    let gridwid = (grid[n - 1] - grid[0]) / (n as f64 - 1.0);
    // integer index guess by uniform spacing
    let mut i = ((val - grid[0]) / gridwid) as usize;
    if i >= n - 1 {
        i = n - 2;
    }
    let a = grid[i];
    let b = grid[i + 1];
    let aa = dens[i];
    let bb = dens[i + 1];
    if val == a {
        return aa;
    }
    if b == a {
        return aa;
    }
    aa + (bb - aa) / (b - a) * (val - a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_at_zero() {
        // φ(0) = 1/sqrt(6.28318) ≈ 0.39894
        assert!((gaussian_pdf(0.0) - 1.0 / 6.28318f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn density_normalization_and_interp() {
        // uniform grid 0..1, 11 points; data all at 0.5 with weight 1 → symmetric peak at 0.5.
        let grid: Vec<f64> = (0..11).map(|i| i as f64 / 10.0).collect();
        let vals = vec![0.5_f64];
        let w = vec![1.0_f64];
        let d = wt_density_fit(&grid, &vals, 0.2, &w);
        // peak at the 0.5 node (index 5)
        let maxi = d.iter().cloned().fold(f64::MIN, f64::max);
        assert!((d[5] - maxi).abs() < 1e-9);
        // interpolation at a node returns the node value
        assert!((predict_density_fit(&grid, &d, 0.5) - d[5]).abs() < 1e-12);
        // midway between nodes is between the two
        let p = predict_density_fit(&grid, &d, 0.55);
        assert!(p <= d[5] && p >= d[6]);
        // outside the grid → 0
        assert_eq!(predict_density_fit(&grid, &d, 2.0), 0.0);
    }
}
