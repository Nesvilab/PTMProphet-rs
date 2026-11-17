// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Number formatting that reproduces C/C++ stream output byte-for-bit.
//!
//! [`fmt_g`] mirrors C `printf("%.*g", prec, x)` — equivalently C++ `ostream <<` with the given
//! `precision()` and no `fixed`/`scientific` flag. PTMProphet relies on the stream default
//! (`prec = 6`) for the `prior` attribute *and* the whole `<mixturemodel>` curve, so we reproduce it
//! exactly. [`fmt_g6`] is the `prec = 6` default used in both of those places.

/// Format like C `printf("%.*g", prec, x)`: `prec` significant digits, trailing zeros trimmed,
/// switching to scientific notation when the decimal exponent is `< -4` or `>= prec`. This is exactly
/// what C++ `ostream <<` does by default (prec = 6), which PTMProphet relies on — so we reproduce it
/// bit-for-bit.
pub fn fmt_g(x: f64, prec: usize) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    if !x.is_finite() {
        return format!("{x}");
    }
    let p = prec.max(1);
    let e = x.abs().log10().floor() as i32;
    if e < -4 || e >= p as i32 {
        // Scientific: p-1 fractional digits; trim mantissa zeros; C-style exponent (e±NN, ≥2 digits).
        let s = format!("{:.*e}", p - 1, x);
        let (mant, exp) = s.split_once('e').unwrap_or((s.as_str(), "0"));
        let mant = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        let exp_i: i32 = exp.parse().unwrap_or(0);
        let sign = if exp_i < 0 { '-' } else { '+' };
        format!("{mant}e{sign}{:02}", exp_i.abs())
    } else {
        // Fixed: (p-1-e) fractional digits; trim trailing zeros and a dangling dot.
        let decimals = (p as i32 - 1 - e).max(0) as usize;
        let s = format!("{:.*}", decimals, x);
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    }
}

/// C++ stream default formatting (`%.6g`) — used for the first `prior` and the `<mixturemodel>` curve.
pub fn fmt_g6(x: f64) -> String {
    fmt_g(x, 6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_cpp_stream_defaults() {
        // 6 significant figures, trailing zeros trimmed.
        assert_eq!(fmt_g6(0.203933072808358), "0.203933");
        assert_eq!(fmt_g6(-0.00795658243679), "-0.00795658");
        assert_eq!(fmt_g6(0.11834683990789886), "0.118347");
        assert_eq!(fmt_g6(-0.01), "-0.01");
        assert_eq!(fmt_g6(0.0), "0");
        // Switches to scientific below 1e-4.
        assert_eq!(fmt_g6(0.0000123456789), "1.23457e-05");
    }
}
