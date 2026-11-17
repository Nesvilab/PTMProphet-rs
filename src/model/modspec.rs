// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Modification-type definition — one per distinct shift string of the `<mods>` argument (tokens
//! with the same shift string are one type).
//!
//! Example: `STY:79.966331` → `aas=['S','T','Y']`, `shift=79.966331`, `ptm_label="STY:79.966331"`.

/// A single variable-modification *type* the localizer handles.
#[derive(Debug, Clone)]
pub struct ModType {
    /// Residues this type may sit on, including `'n'`/`'c'` for terminal mods. Order as in the input;
    /// a token that joins this type appends its residues.
    pub aas: Vec<char>,
    /// Mass shift (Δ), e.g. `79.966331`. Parsed from the shift substring with C-locale numeric parsing.
    pub shift: f64,
    /// `"AAs:shiftStr"` label → emitted as the `ptm="..."` attribute and ROC key. `shiftStr` is the
    /// *verbatim* shift substring from the input (string-equality merges types); a token that joins
    /// this type puts its residues at the front, as PTMProphet does (`nSTY:79.966331`).
    pub ptm_label: String,
}

impl ModType {
    /// True if `aa` is in this type's residue set.
    pub fn can_modify(&self, aa: char) -> bool {
        self.aas.contains(&aa)
    }

    /// True if this type targets the peptide **N-terminus** (`'n'` in its residue set), e.g. n-term
    /// acetylation `n:42.0106`.
    pub fn targets_nterm(&self) -> bool {
        self.aas.contains(&'n')
    }

    /// True if this type targets the peptide **C-terminus** (`'c'` in its residue set).
    pub fn targets_cterm(&self) -> bool {
        self.aas.contains(&'c')
    }
}
