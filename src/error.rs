// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Error type for PTMProphet-rs.
//!
//! PTMProphet-rs validates at boundaries and returns a `PtmProphetRsError` so the binary can print
//! a clear message and exit non-zero. Conditions that are fatal by design (e.g. static-mod mass
//! conflict, NaN probability) are mapped to explicit error variants.

use std::fmt;
use std::path::PathBuf;

pub type Result<T> = std::result::Result<T, PtmProphetRsError>;

#[derive(Debug)]
pub enum PtmProphetRsError {
    /// CLI usage error (e.g. no arguments supplied).
    Usage(String),
    /// A recognized option that is not implemented in this build.
    /// We fail loudly rather than silently mis-localize.
    UnsupportedOption(String),
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    /// Malformed pepXML (missing required attribute, etc.).
    PepXml(String),
    /// mzML parse/lookup failure (e.g. scan not found, bad <binary>).
    Mzml(String),
    /// Static-mod declared twice with conflicting mass.
    StaticModConflict { residue: char, m1: f64, m2: f64 },
    /// A NaN propagated into a probability computation.
    NaN(String),
}

impl fmt::Display for PtmProphetRsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PtmProphetRsError::Usage(m) => write!(f, "usage: {m}"),
            PtmProphetRsError::UnsupportedOption(m) => write!(f, "unsupported option in PTMProphet-rs build: {m}"),
            PtmProphetRsError::Io { path, source } => write!(f, "I/O error on {}: {source}", path.display()),
            PtmProphetRsError::PepXml(m) => write!(f, "pepXML error: {m}"),
            PtmProphetRsError::Mzml(m) => write!(f, "mzML error: {m}"),
            PtmProphetRsError::StaticModConflict { residue, m1, m2 } => {
                write!(f, "conflicting static modification on '{residue}': {m1} vs {m2}")
            }
            PtmProphetRsError::NaN(m) => write!(f, "NaN encountered: {m}"),
        }
    }
}

impl std::error::Error for PtmProphetRsError {}
