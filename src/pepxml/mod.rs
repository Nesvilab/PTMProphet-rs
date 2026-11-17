// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! pepXML I/O: a `Tag` mini-parser/writer, the PASS-A index/scan
//! (`parseSpecProbs`), the PASS-C localizing writer (`parseWriteUpdateMpx` +
//! `writeUpdatedModTags`), and the ROC summary (`insertROC`).

pub mod tag;
pub mod index;
pub mod psm;
pub mod writer;
pub mod roc;
