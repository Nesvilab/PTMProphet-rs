// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! Expectation–Maximization recalibration. A 1-D non-parametric mixture model ([`kdmodel`])
//! over the per-site Oscores, fit by weighted Gaussian KDE ([`density`]) and iterated by
//! [`train::compute_ptm_model`]. When EM is enabled (EM=1) only the intensity Oscore model is
//! live; all other models are deleted.

pub mod kdmodel;
pub mod density;
pub mod train;
