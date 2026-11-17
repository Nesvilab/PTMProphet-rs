// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Spectrum I/O: an indexed mzML reader that returns a [`mzml::Spectrum`] by scan number.
//!
//! Only the spectrum fields needed for scoring are read: the MS2 peak list (m/z + intensity), MS
//! level, precursor m/z and charge, and the scan number. mzXML/mgf are out of scope.

pub mod mzml;
