// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2006 Henry Lam
// Translated into Rust from TPP 6.3.2 src/Search/SpectraST/Peptide.cpp (LGPL-2.1-or-later);
// modified by the University of Michigan, 2026.

//! Peptide model. Parses a modified-peptide string into a stripped sequence + per-position mod
//! Δmasses, computes the neutral monoisotopic mass, and emits theoretical fragment-ion m/z.
//!
//! Conventions:
//! * Residue masses are *residue* masses; the whole peptide gains H ('n') + OH ('c') = one water.
//! * `mods` iterates ascending by position (BTreeMap) for reproducible summation.
//! * b = Σ(prefix residues+mods) + nterm_Δ + z·PROTON_FRAG;  a = b − CO;  c = b + NH3;
//!   y = Σ(suffix residues+mods) + cterm_Δ + H2O + z·PROTON_FRAG;  z = y − (NH3 − proton);
//!   x = y + CO − 2H. Which series are emitted is set by `--n-ions`/`--c-ions`.
//! * Stored mod values are **Δmass** (so `residue + Δ` reconstructs the bracketed absolute mass).
//!   Placement masses come from the command-line `shift` (set via [`Peptide::set_mod`]); parsing a
//!   modified-peptide string yields the (possibly integer-rounded) bracket-derived Δ used only to
//!   *identify* modified positions.

use std::collections::BTreeMap;

use crate::error::{PtmProphetRsError, Result};
use crate::model::masses::{self, IonKind};

/// A parsed (possibly modified) peptide.
#[derive(Debug, Clone, Default)]
pub struct Peptide {
    /// Bare residue letters (no n/c, no brackets), uppercase.
    pub stripped: String,
    /// 0-based residue index → modification Δmass (ascending iteration). Terminal mods are kept
    /// separately, never as negative keys.
    pub mods: BTreeMap<usize, f64>,
    /// N-terminal mod Δmass (on top of the H group), if any.
    pub nterm_mod: Option<f64>,
    /// C-terminal mod Δmass (on top of the OH group), if any.
    pub cterm_mod: Option<f64>,
    pub charge: i32,
    /// Bracket-rendering precision (0 in our config).
    pub mod_prec: u32,
    /// Cached residue chars for O(1) indexed mass lookups.
    residues: Vec<char>,
}

impl Peptide {
    /// Parse a modified-peptide string: an optional leading `n[mass]`, residues each optionally
    /// followed by `[mass]`, an optional trailing `c[mass]`, and an optional trailing `/charge`.
    /// Bracketed values are absolute residue+mod (or group+mod for n/c) masses; stored as Δ from
    /// the base residue/group.
    ///
    /// (ProForma `[mod]-…`/`…-[mod]` and `P.SEQ.N` flanks are not handled by this parser; a TODO
    /// marks adding them if an input ever needs them.)
    pub fn parse(modified_peptide: &str, charge: i32, mod_prec: u32) -> Result<Peptide> {
        let mut p = Peptide { charge, mod_prec, ..Default::default() };
        let bytes: Vec<char> = modified_peptide.chars().collect();
        let n = bytes.len();
        let mut i = 0usize;

        // Optional trailing "/charge".
        let mut end = n;
        if let Some(slash) = modified_peptide.rfind('/') {
            if charge == 0 {
                if let Ok(z) = modified_peptide[slash + 1..].trim().parse::<i32>() {
                    p.charge = z;
                }
            }
            end = slash;
        }

        // Helper: read a bracketed mass starting at '[' (bytes[j]=='['); returns (value, index after ']').
        let read_bracket = |bytes: &[char], mut j: usize| -> Result<(f64, usize)> {
            debug_assert_eq!(bytes[j], '[');
            j += 1;
            let start = j;
            while j < bytes.len() && bytes[j] != ']' {
                j += 1;
            }
            if j >= bytes.len() {
                return Err(PtmProphetRsError::PepXml("unterminated '[' in peptide".to_string()));
            }
            let s: String = bytes[start..j].iter().collect();
            let v: f64 = s
                .trim()
                .parse()
                .map_err(|_| PtmProphetRsError::PepXml(format!("bad bracket mass '{s}'")))?;
            Ok((v, j + 1)) // skip ']'
        };

        while i < end {
            let c = bytes[i];
            match c {
                'n' if i + 1 < end && bytes[i + 1] == '[' => {
                    let (mass, ni) = read_bracket(&bytes, i + 1)?;
                    p.nterm_mod = Some(mass - masses::GROUP_N_TERM);
                    i = ni;
                }
                'c' if i + 1 < end && bytes[i + 1] == '[' => {
                    let (mass, ni) = read_bracket(&bytes, i + 1)?;
                    p.cterm_mod = Some(mass - masses::GROUP_C_TERM);
                    i = ni;
                }
                _ if c.is_ascii_alphabetic() => {
                    let aa = c.to_ascii_uppercase();
                    let res = masses::residue_mono_mass(aa)
                        .ok_or_else(|| PtmProphetRsError::PepXml(format!("unknown residue '{aa}'")))?;
                    let idx = p.residues.len();
                    p.residues.push(aa);
                    p.stripped.push(aa);
                    i += 1;
                    if i < end && bytes[i] == '[' {
                        let (mass, ni) = read_bracket(&bytes, i)?;
                        p.mods.insert(idx, mass - res);
                        i = ni;
                    }
                }
                // Skip any other separators defensively.
                _ => {
                    i += 1;
                }
            }
        }
        Ok(p)
    }

    /// Number of residues (no n/c).
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// The residue characters (read-only).
    pub fn residues(&self) -> &[char] {
        &self.residues
    }

    /// Set (or replace) the Δmass at 0-based residue position `pos`. Used to build a concrete
    /// placement with the command-line `shift`.
    pub fn set_mod(&mut self, pos: usize, delta: f64) {
        self.mods.insert(pos, delta);
    }

    pub fn clear_mods(&mut self) {
        self.mods.clear();
    }

    /// Neutral monoisotopic peptide mass = Σ residues (ascending) + H('n') + OH('c') + terminal Δ
    /// + Σ per-residue Δ (ascending).
    pub fn monoisotopic_neutral_mass(&self) -> f64 {
        let mut sum = 0.0;
        for (i, &aa) in self.residues.iter().enumerate() {
            sum += masses::residue_mono_mass(aa).unwrap_or(0.0);
            if let Some(&d) = self.mods.get(&i) {
                sum += d;
            }
        }
        sum += masses::GROUP_N_TERM; // H
        sum += masses::GROUP_C_TERM; // OH
        if let Some(d) = self.nterm_mod {
            sum += d;
        }
        if let Some(d) = self.cterm_mod {
            sum += d;
        }
        sum
    }

    /// Fragment m/z for `kind`, taking `num_aa` residues (prefix for N-term ions a/b/c, suffix for
    /// C-term ions x/y/z) at `charge`. Returns 0.0 unless `1 <= num_aa < len` and
    /// `charge >= 1`.
    pub fn monoisotopic_mz_fragment(&self, kind: IonKind, num_aa: usize, charge: u32) -> f64 {
        let len = self.len();
        if num_aa < 1 || num_aa >= len || charge < 1 {
            return 0.0;
        }
        let z = charge as f64;
        match kind {
            IonKind::B | IonKind::A | IonKind::C => {
                // N-terminal: first `num_aa` residues.
                let mut sum = self.nterm_mod.unwrap_or(0.0);
                for i in 0..num_aa {
                    sum += masses::residue_mono_mass(self.residues[i]).unwrap_or(0.0);
                    if let Some(&d) = self.mods.get(&i) {
                        sum += d;
                    }
                }
                sum += z * masses::PROTON_FRAG;
                match kind {
                    IonKind::A => sum -= masses::MASS_CO,
                    IonKind::C => sum += masses::MASS_NH3,
                    _ => {}
                }
                sum / z
            }
            IonKind::Y | IonKind::Z | IonKind::X => {
                // C-terminal: last `num_aa` residues.
                let mut sum = self.cterm_mod.unwrap_or(0.0);
                for i in (len - num_aa)..len {
                    sum += masses::residue_mono_mass(self.residues[i]).unwrap_or(0.0);
                    if let Some(&d) = self.mods.get(&i) {
                        sum += d;
                    }
                }
                sum += masses::MASS_H2O;
                sum += z * masses::PROTON_FRAG;
                match kind {
                    // z• = y − (NH3 − proton)
                    IonKind::Z => sum -= masses::MASS_NH3 - masses::PROTON_FRAG,
                    // x = y + CO − 2H (chemically correct ≈ y + 25.979). NOTE: PTMProphet's x uses
                    // a hydroxyl (OH) where carbon (C) belongs, making its x-ions ~5 Da heavy; we
                    // deliberately do not replicate that.
                    IonKind::X => sum += masses::MASS_CO - 2.0 * masses::MASS_H,
                    _ => {}
                }
                sum / z
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-4;

    #[test]
    fn parse_tmt_phospho_peptide() {
        // scan 3185 input modified_peptide, full-precision brackets.
        let p = Peptide::parse("n[230.17076]ALRAS[166.99835]ESGI", 2, 0).unwrap();
        assert_eq!(p.stripped, "ALRASESGI");
        assert_eq!(p.len(), 9);
        // n-term TMT Δ = 230.17076 - H(1.007825035) ≈ 229.16293
        assert!((p.nterm_mod.unwrap() - 229.16293).abs() < EPS, "{:?}", p.nterm_mod);
        // S at 0-based index 4 (residue 5) carries phospho Δ ≈ 79.96633.
        assert!((p.mods[&4] - 79.96633).abs() < 1e-3, "{:?}", p.mods.get(&4));
        assert!(!p.mods.contains_key(&6)); // the other S (index 6) is unmodified in the input
    }

    #[test]
    fn residue_plus_phospho_matches_oracle_brackets() {
        // Oracle absolute masses (mass_aminoacid_mass / search_summary).
        assert!((masses::residue_mono_mass('S').unwrap() + 79.96633 - 166.99835).abs() < 1e-4);
        assert!((masses::residue_mono_mass('T').unwrap() + 79.96633 - 181.014).abs() < 1e-4);
        assert!((masses::residue_mono_mass('M').unwrap() + 15.9949 - 147.0354).abs() < 1e-4);
        // n-term TMT group + Δ.
        assert!((masses::GROUP_N_TERM + 229.16293 - 230.17076).abs() < 1e-4);
    }

    #[test]
    fn neutral_mass_unmodified() {
        // PEPTIDE neutral monoisotopic mass = Σresidues + H2O.
        let p = Peptide::parse("PEPTIDE", 0, 0).unwrap();
        let expected: f64 = ['P', 'E', 'P', 'T', 'I', 'D', 'E']
            .iter()
            .map(|&c| masses::residue_mono_mass(c).unwrap())
            .sum::<f64>()
            + masses::MASS_H2O;
        assert!((p.monoisotopic_neutral_mass() - expected).abs() < EPS);
    }

    #[test]
    fn b_and_y_ions_unmodified() {
        let p = Peptide::parse("PEPTIDE", 1, 0).unwrap();
        // b2 (PE, +1) = P + E + proton.
        let b2 = masses::residue_mono_mass('P').unwrap()
            + masses::residue_mono_mass('E').unwrap()
            + masses::PROTON_FRAG;
        assert!((p.monoisotopic_mz_fragment(IonKind::B, 2, 1) - b2).abs() < EPS);
        // y2 (DE, +1) = D + E + H2O + proton.
        let y2 = masses::residue_mono_mass('D').unwrap()
            + masses::residue_mono_mass('E').unwrap()
            + masses::MASS_H2O
            + masses::PROTON_FRAG;
        assert!((p.monoisotopic_mz_fragment(IonKind::Y, 2, 1) - y2).abs() < EPS);
        // a2 = b2 - CO.
        assert!(
            (p.monoisotopic_mz_fragment(IonKind::A, 2, 1) - (b2 - masses::MASS_CO)).abs() < EPS
        );
        // out-of-range guards.
        assert_eq!(p.monoisotopic_mz_fragment(IonKind::B, 0, 1), 0.0);
        assert_eq!(p.monoisotopic_mz_fragment(IonKind::B, 7, 1), 0.0); // num_aa == len
    }

    #[test]
    fn a_c_x_z_ions_unmodified() {
        let p = Peptide::parse("PEPTIDE", 1, 0).unwrap();
        // b2 (PE, +1) and y2 (DE, +1) as references.
        let b2 = masses::residue_mono_mass('P').unwrap()
            + masses::residue_mono_mass('E').unwrap()
            + masses::PROTON_FRAG;
        let y2 = masses::residue_mono_mass('D').unwrap()
            + masses::residue_mono_mass('E').unwrap()
            + masses::MASS_H2O
            + masses::PROTON_FRAG;
        // a = b − CO ; c = b + NH3 (N-terminal series).
        assert!((p.monoisotopic_mz_fragment(IonKind::A, 2, 1) - (b2 - masses::MASS_CO)).abs() < EPS);
        assert!((p.monoisotopic_mz_fragment(IonKind::C, 2, 1) - (b2 + masses::MASS_NH3)).abs() < EPS);
        // z• = y − (NH3 − proton) ; x = y + CO − 2H (C-terminal series).
        let z2 = y2 - (masses::MASS_NH3 - masses::PROTON_FRAG);
        assert!((p.monoisotopic_mz_fragment(IonKind::Z, 2, 1) - z2).abs() < EPS);
        let x_off = masses::MASS_CO - 2.0 * masses::MASS_H;
        assert!((p.monoisotopic_mz_fragment(IonKind::X, 2, 1) - (y2 + x_off)).abs() < EPS);
        // The x offset is ~+25.979 (CO − H2), NOT PTMProphet's OH-based ~+30.98.
        assert!((x_off - 25.9793).abs() < 1e-3, "x offset {x_off}");
    }

    #[test]
    fn doubly_charged_fragment_halves() {
        let p = Peptide::parse("PEPTIDE", 2, 0).unwrap();
        let b3_1 = p.monoisotopic_mz_fragment(IonKind::B, 3, 1);
        let b3_2 = p.monoisotopic_mz_fragment(IonKind::B, 3, 2);
        // (neutral + 2H)/2 vs (neutral + 1H)/1 → check the +1 neutral relation.
        let neutral_plus_1 = b3_1; // = neutral_b + proton
        let expected_2 = (neutral_plus_1 - masses::PROTON_FRAG + 2.0 * masses::PROTON_FRAG) / 2.0;
        assert!((b3_2 - expected_2).abs() < EPS);
    }

    #[test]
    fn set_mod_uses_exact_shift() {
        // A placement uses the command-line shift, not the integer bracket.
        let mut p = Peptide::parse("ALRASESGI", 2, 0).unwrap();
        p.set_mod(4, 79.96633); // phospho on S5 (index 4)
        let m = p.monoisotopic_neutral_mass();
        let base = Peptide::parse("ALRASESGI", 2, 0).unwrap().monoisotopic_neutral_mass();
        assert!((m - base - 79.96633).abs() < EPS);
    }
}
