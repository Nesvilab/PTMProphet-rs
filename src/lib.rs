// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! PTMProphet-rs — a Rust tool for PTM site localization in shotgun proteomics.
//!
//! See `README.md` for the architecture and per-module responsibilities.
//!
//! A typical invocation:
//! `ptmprophet-rs --mods <mods> --frag-ppm-tol 10 --em 1 --n-ions b --min-prob 0.5 --no-stack
//! --static --filelist filelist_ptmloc.txt`.

// A few functions take many positional arguments (the per-placement matcher threads several
// accumulators); grouping them into structs would not make them clearer, so silence that one lint.
#![allow(clippy::too_many_arguments)]

pub mod config;
pub mod error;
pub mod fmt;
pub mod pipeline;
pub mod threads;
pub mod validate;

pub mod model;
pub mod pepxml;
pub mod spectrum;
pub mod mpx;
pub mod em;

pub use config::Config;
pub use error::{PtmProphetRsError, Result};
