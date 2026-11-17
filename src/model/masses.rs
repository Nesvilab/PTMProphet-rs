// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2006 Henry Lam
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Search/SpectraST/Analyte.cpp,
// src/Search/SpectraST/Peptide.cpp and src/Validation/PTMProphetParser/PTMProphetMpx.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! Monoisotopic mass constants. These are *reference data* (not algorithm).
//!
//! CRITICAL: there are **two** proton masses:
//!   * [`PROTON_FRAG`] = 1.00728        — used by ALL a/b/c/x/y/z fragment m/z.
//!   * [`PROTON_PRECURSOR`] = 1.007276466 — used ONLY for the precursor/isotope envelope.
//!
//! They must not be unified.

/// Proton mass used for fragment-ion m/z.
pub const PROTON_FRAG: f64 = 1.00728;
/// Proton mass used only for the precursor / isotope envelope.
pub const PROTON_PRECURSOR: f64 = 1.007276466;

// Element / group monoisotopic masses.
pub const MASS_C: f64 = 12.00000;
pub const MASS_H: f64 = 1.007825035;
pub const MASS_O: f64 = 15.99491463;
pub const MASS_N: f64 = 14.003074;
pub const MASS_H2O: f64 = 18.01056;
pub const MASS_NH3: f64 = 17.026549;
pub const MASS_CO: f64 = 27.994915;

/// N-terminus group ('n' pseudo-residue) = a single H.
pub const GROUP_N_TERM: f64 = MASS_H; // 1.007825035
/// C-terminus group ('c' pseudo-residue) = OH = H + O.
pub const GROUP_C_TERM: f64 = MASS_H + MASS_O; // 17.00273963...

/// Monoisotopic *residue* masses (AA minus water).
/// NOTE: `X` (113.08406) is NOT bit-identical to `L`/`I` (113.084064) — keep distinct.
/// Returns `None` for an unrecognized character; validation happens at the boundary.
///
/// The residues that actually carry a modification in the supported configs (`S`/`T`/`Y` phospho,
/// `M` oxidation, `C` carbamidomethyl, `K` TMT) use PTMProphet's full precision so the emitted
/// `mass="…"` (residue + mod shift, formatted to 8 dp) is bit-identical — e.g.
/// `S(87.03202844) + 79.96633 = 166.99835844`. Truncating these to 6 dp shifted the 7th–8th decimal
/// of every `mod_aminoacid_mass`. Masses feed fragment m/z too, but there the 10 ppm tolerance makes
/// the extra digits irrelevant.
pub fn residue_mono_mass(aa: char) -> Option<f64> {
    let m = match aa {
        'A' => 71.037114,
        'B' => 114.53494, // avg Asp/Asn
        'C' => 103.00918450,
        'D' => 115.026943,
        'E' => 129.042593,
        'F' => 147.068414,
        'G' => 57.021464,
        'H' => 137.058912,
        'I' => 113.084064,
        'J' => 0.0, // fake AA
        'K' => 128.09496305,
        'L' => 113.084064,
        'M' => 131.04048464,
        'N' => 114.042927,
        'O' => 237.14773, // pyrrolysine
        'P' => 97.052764,
        'Q' => 128.058578,
        'R' => 156.101111,
        'S' => 87.03202844,
        'T' => 101.04767850,
        'U' => 150.95363, // selenocysteine
        'V' => 99.068414,
        'W' => 186.079313,
        'X' => 113.08406, // set to Leu/Ile but 5-dp value, intentionally distinct
        'Y' => 163.06332858,
        'Z' => 128.55059, // avg Glu/Gln
        _ => return None,
    };
    Some(m)
}

/// Backbone fragment-ion series. N-terminal (prefix) ions are `a`/`b`/`c`; C-terminal (suffix)
/// ions are `x`/`y`/`z`. Which series are actually produced is driven by `--n-ions`/`--c-ions`
/// (default `b` / `y`); the matcher iterates exactly the configured set (see
/// [`crate::mpx::Mpx::evaluate_placement`]).
///
/// Mass relationships (realized by the formulas in [`crate::model::peptide`]):
///   * `a = b − CO`,  `c = b + NH3`
///   * `z = y − (NH3 − proton)`  (z• radical, as seen in ETD)
///   * `x = y + CO − 2H`
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IonKind {
    A,
    B,
    C,
    X,
    Y,
    Z,
}

impl IonKind {
    /// Map an ion letter (case-insensitive) to its kind; `None` if unrecognized.
    pub fn from_char(c: char) -> Option<IonKind> {
        match c.to_ascii_lowercase() {
            'a' => Some(IonKind::A),
            'b' => Some(IonKind::B),
            'c' => Some(IonKind::C),
            'x' => Some(IonKind::X),
            'y' => Some(IonKind::Y),
            'z' => Some(IonKind::Z),
            _ => None,
        }
    }

    /// The canonical lowercase letter for this kind (inverse of [`IonKind::from_char`]).
    pub fn to_char(self) -> char {
        match self {
            IonKind::A => 'a',
            IonKind::B => 'b',
            IonKind::C => 'c',
            IonKind::X => 'x',
            IonKind::Y => 'y',
            IonKind::Z => 'z',
        }
    }

    /// True for N-terminal (prefix) series `a`/`b`/`c`; false for C-terminal `x`/`y`/`z`.
    pub fn is_nterm(self) -> bool {
        matches!(self, IonKind::A | IonKind::B | IonKind::C)
    }
}
