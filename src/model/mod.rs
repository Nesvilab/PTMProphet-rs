// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Peptide / mass model: monoisotopic mass constants, modification-type definitions, and
//! the `Peptide` engine that parses modified-peptide strings and computes theoretical fragment m/z.

pub mod masses;
pub mod modspec;
pub mod peptide;
