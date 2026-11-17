// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetMain.cpp and
// src/Validation/PTMProphetParser/PTMProphetParser.cpp (LGPL-2.1-or-later); modified by the
// University of Michigan, 2026.

//! Command-line parsing (standard Unix-style flags via clap) and the resolved run configuration.
//!
//! [`Cli`] is the clap front end — `ptmprophet-rs --mods … [OPTIONS] --filelist <filelist>`.
//! [`Config::from_cli`] resolves it into the immutable [`Config`] the pipeline consumes and, in the
//! process, reconstructs the canonical `options="…"` string emitted verbatim into the output pepXML
//! header.

use std::path::{Path, PathBuf};

use clap::Parser;

use crate::error::{PtmProphetRsError, Result};
use crate::model::masses::IonKind;
use crate::model::modspec::ModType;

/// Standard Unix-style command line.
///
/// Input is always a file-list: `ptmprophet-rs --mods … --filelist filelist_ptmloc.txt`. Every
/// `interact-*.pep.xml` it names is localized in one PTMProphet-rs run, but each file trains its
/// own EM intensity model on its own evidence (matching PTMProphet's per-file invocation). Each is
/// written to a sibling `interact-*.mod.pep.xml`.
#[derive(Parser, Debug)]
#[command(
    name = "ptmprophet-rs",
    version,
    about = "PTM site localization: scores variable-modification placement on PSMs from MS/MS fragment evidence."
)]
pub struct Cli {
    /// Modification spec: comma-separated `AAs:shift` (e.g. `M:15.9949,STY:79.96633`). Neutral-loss
    /// suffixes (`AAs:shift:nl`) are not implemented and are rejected.
    #[arg(long, value_name = "MODS")]
    pub mods: String,

    /// Fragment-ion match tolerance, in ppm (default 15, PTMProphet's default under `STATIC`).
    #[arg(long, value_name = "PPM")]
    pub frag_ppm_tol: Option<f64>,

    /// EM mode (only `1` is implemented).
    #[arg(long, value_name = "N", default_value_t = 1)]
    pub em: u8,

    /// Minimum PSM probability to localize.
    #[arg(long, value_name = "P", default_value_t = 0.5)]
    pub min_prob: f64,

    /// N-terminal ion series: `a`, `b` or `c`, comma-separated (e.g. `b` or `a,b`).
    #[arg(long, value_name = "IONS", default_value = "b")]
    pub n_ions: String,

    /// C-terminal ion series: `x`, `y` or `z`, comma-separated (e.g. `y` or `y,z`).
    #[arg(long, value_name = "IONS", default_value = "y")]
    pub c_ions: String,

    /// Modification mass decimal precision (only `0` is implemented).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub mod_prec: u32,

    /// O-score pseudo-counts (only `0` is implemented: PTMProphet's MINO correction is not).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub min_o: i32,

    /// Maximum fragment charge: the precursor charge plus Z, at least 1. Only 0 (default) or below
    /// is implemented; give a negative value as `--max-frag-z=-1`.
    #[arg(long, value_name = "Z", default_value_t = 0)]
    pub max_frag_z: i32,

    /// Forbid two modifications on the same residue (PTMProphet `NOSTACK`).
    #[arg(long)]
    pub no_stack: bool,

    /// Match every PSM with the one `--frag-ppm-tol` (PTMProphet `STATIC`). Required: PTMProphet's
    /// per-PSM estimate of the fragment mass offset and tolerance is not implemented.
    #[arg(long = "static")]
    pub static_tol: bool,

    /// Worker-thread count (default: all logical CPUs). `PTMPROPHET_RS_THREADS` overrides this.
    #[arg(long, short = 'j', value_name = "N")]
    pub threads: Option<usize>,

    /// File-list: a text file naming the `interact-*.pep.xml` inputs and their mzML spectrum files to
    /// localize (each file trains its own model). One input per non-blank line (`#` comments); each
    /// line is the pepXML path followed by exactly one mzML path (tab- or whitespace-separated).
    #[arg(long, value_name = "FILE")]
    pub filelist: PathBuf,
}

/// One input pepXML to localize, resolved from a single file-list entry. All entries are processed
/// in a single PTMProphet-rs run, each with its own EM model trained on its own evidence.
#[derive(Debug, Clone)]
pub struct InputSpec {
    /// The `interact-*.pep.xml` to read.
    pub pep_xml: PathBuf,
    /// The mzML spectrum file for this input, taken verbatim from the file-list. Every PSM in
    /// `pep_xml` is localized against this file; the path is never read from inside the pepXML.
    pub mzml: PathBuf,
    /// Where the localized output is written.
    pub out: PathBuf,
}

/// Fully resolved configuration for one localization run.
#[derive(Debug, Clone)]
pub struct Config {
    // --- modifications ---
    /// One entry per distinct shift string of `--mods` (e.g. `STY:79.966331`, `M:15.9949`), in input
    /// order; tokens with the same shift string are one type. Built by [`parse_mod_string`].
    pub mods: Vec<ModType>,
    /// The `--mods` argument verbatim → header `mod_string="..."`.
    pub mod_string: String,
    /// The PTMProphet command-line tokens (`NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 ...`) that stand
    /// for the resolved options → header `options="..."`. Built by `build_opt_string`: `KEEPOLD`,
    /// `STATIC`, `FRAGPPMTOL`, `EM`, `NIONS`, the mods, `MINPROB` and `MAXTHREADS=1` are always
    /// written; `NOSTACK` only with `--no-stack`, `CIONS` only when it is not `y`, and `MAXFRAGZ`
    /// only when it is not 0. PTMProphet-rs's `--threads` and the file paths are not written.
    pub opt_string: String,

    // --- tolerances ---
    pub frag_tol: f64, // --frag-ppm-tol, in ppm

    // --- modes / knobs ---
    pub em: u8,          // --em (=1)
    pub min_prob: f64,   // --min-prob (=0.5)
    pub mod_prec: u32,   // --mod-prec (=0)
    pub min_o: i32,      // --min-o (=0)
    pub max_frag_z: i32, // --max-frag-z (=0)
    /// Explicit `--threads N` (None = all logical CPUs). Wired into the worker pool at startup.
    pub threads: Option<usize>,

    pub nostack: bool,    // --no-stack
    pub static_tol: bool, // --static

    // --- ion series ---
    pub nions: Vec<IonKind>, // --n-ions (default [B]); N-terminal a/b/c
    pub cions: Vec<IonKind>, // --c-ions (default [Y]); C-terminal x/y/z

    // --- files ---
    /// Every input pepXML named by the file-list. The pipeline iterates this list and localizes each
    /// input with its own EM model.
    pub inputs: Vec<InputSpec>,
    /// The file-list path. Retained for messages and the lower-level test/validation helpers (which
    /// set it directly to a pepXML); the pipeline drives off `inputs`.
    pub in_file: PathBuf,
    /// First input's output path. Retained for the test/validation helpers; the pipeline drives off
    /// `inputs`.
    pub out_file: PathBuf,
}

impl Config {
    /// Resolve a parsed [`Cli`] into the immutable run [`Config`].
    ///
    /// Steps: reject the options PTMProphet-rs does not implement (`--em` other than 1, a run
    /// without `--static`, `--min-o` or `--mod-prec` other than 0, a negative `--frag-ppm-tol`, a
    /// positive `--max-frag-z`); resolve the fragment tolerance (15 ppm unless given); split the
    /// ion-series strings; parse `--mods`; reconstruct the `options="…"` string; parse the
    /// `--filelist` into the input set.
    pub fn from_cli(cli: Cli) -> Result<Config> {
        if cli.em != 1 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--em {} (only --em 1 is implemented)",
                cli.em
            )));
        }
        // Without STATIC, PTMProphet estimates a fragment mass offset and tolerance for each PSM.
        if !cli.static_tol {
            return Err(PtmProphetRsError::UnsupportedOption(
                "a run without --static (PTMProphet's per-PSM estimate of the fragment mass offset \
                 and tolerance is not implemented; pass --static)"
                    .into(),
            ));
        }
        // A nonzero MINO also makes PTMProphet adjust mean_best_prob by its MINO factors (which it
        // computes only for MINO above 0); PTMProphet-rs has no MINO correction.
        if cli.min_o != 0 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--min-o {} (only --min-o 0 is implemented)",
                cli.min_o
            )));
        }
        // PTMProphet-rs matches and writes modification masses at PTMProphet's MODPREC=0 only.
        if cli.mod_prec != 0 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--mod-prec {} (only --mod-prec 0 is implemented)",
                cli.mod_prec
            )));
        }
        // A positive MAXFRAGZ makes PTMProphet cap the fragment charge at Z; PTMProphet-rs has no
        // such cap.
        if cli.max_frag_z > 0 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--max-frag-z {} (a positive --max-frag-z is not implemented in PTMProphet-rs)",
                cli.max_frag_z
            )));
        }
        // PTMProphet replaces a negative FRAGPPMTOL with 15 or 200 ppm, by where STATIC stands.
        if let Some(v) = cli.frag_ppm_tol {
            if v.is_nan() || v < 0.0 {
                return Err(PtmProphetRsError::UnsupportedOption(format!(
                    "--frag-ppm-tol {v} (the tolerance must be 0 or more)"
                )));
            }
        }

        // Fragment tolerance: explicit flag wins; else PTMProphet's STATIC default of 15 ppm.
        let frag_tol = cli.frag_ppm_tol.unwrap_or(15.0);

        // N-terminal series must be a/b/c, C-terminal x/y/z, each at least one letter; reject
        // anything else at the boundary.
        let nions = parse_ion_series(&cli.n_ions, true)?;
        let cions = parse_ion_series(&cli.c_ions, false)?;

        let mods = parse_mod_string(&cli.mods)?;
        let mod_string = cli.mods.clone();
        let opt_string = build_opt_string(&cli, frag_tol, &nions, &cions);

        // The input set is the file-list's entries (each localized with its own per-file model).
        let inputs = parse_filelist(&cli.filelist)?;
        // Back-compat scalars: the file-list path / first entry. The pipeline drives off `inputs`.
        let in_file = cli.filelist.clone();
        let out_file = inputs.first().map(|s| s.out.clone()).unwrap_or_default();

        Ok(Config {
            mods,
            mod_string,
            opt_string,
            frag_tol,
            em: cli.em,
            min_prob: cli.min_prob,
            mod_prec: cli.mod_prec,
            min_o: cli.min_o,
            max_frag_z: cli.max_frag_z,
            threads: cli.threads,
            nostack: cli.no_stack,
            static_tol: cli.static_tol,
            nions,
            cions,
            inputs,
            in_file,
            out_file,
        })
    }
}

/// Parse a file-list into one [`InputSpec`] per `interact-*.pep.xml`.
///
/// Format: one entry per non-blank line; `#` comment lines are skipped. A line is split into fields
/// on **tabs** (so paths may contain spaces) or, if it has no tab, on whitespace. Each field is
/// classified by extension: each line must name **exactly one** `*.pep.xml`/`*.pepXML` input **and
/// exactly one** `*.mzML` spectrum file — that mzML is where the input's spectra are read from (the
/// path is never taken from inside the pepXML). The output path defaults to the input with
/// `.pep.xml` → `.mod.pep.xml`. Example:
/// ```text
/// # one pep.xml + its mzML per line (tab- or space-separated)
/// /data/interact-f01.pep.xml    /raw/f01_calibrated.mzML
/// /data/interact-f02.pep.xml    /raw/f02_calibrated.mzML
/// ```
pub fn parse_filelist(path: &Path) -> Result<Vec<InputSpec>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| PtmProphetRsError::Io { path: path.to_path_buf(), source: e })?;

    let mut inputs = Vec::new();
    for (i, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = if line.contains('\t') {
            line.split('\t').map(str::trim).filter(|s| !s.is_empty()).collect()
        } else {
            line.split_whitespace().collect()
        };

        let mut pep_xml: Option<PathBuf> = None;
        let mut mzmls: Vec<PathBuf> = Vec::new();
        for tok in fields {
            let lower = tok.to_ascii_lowercase();
            if lower.ends_with(".pep.xml") || lower.ends_with(".pepxml") {
                if pep_xml.is_some() {
                    return Err(PtmProphetRsError::Usage(format!(
                        "file-list {} line {}: more than one pep.xml on a line",
                        path.display(),
                        i + 1
                    )));
                }
                pep_xml = Some(PathBuf::from(tok));
            } else if lower.ends_with(".mzml") {
                mzmls.push(PathBuf::from(tok));
            } else {
                return Err(PtmProphetRsError::Usage(format!(
                    "file-list {} line {}: unrecognized entry '{}' (expected *.pep.xml or *.mzML)",
                    path.display(),
                    i + 1,
                    tok
                )));
            }
        }
        let pep_xml = pep_xml.ok_or_else(|| {
            PtmProphetRsError::Usage(format!(
                "file-list {} line {}: no pep.xml file named",
                path.display(),
                i + 1
            ))
        })?;
        let mzml = match mzmls.len() {
            1 => mzmls.into_iter().next().unwrap(),
            0 => {
                return Err(PtmProphetRsError::Usage(format!(
                    "file-list {} line {}: no mzML file named (each pep.xml must be paired with exactly one *.mzML)",
                    path.display(),
                    i + 1
                )))
            }
            n => {
                return Err(PtmProphetRsError::Usage(format!(
                    "file-list {} line {}: {} mzML files named (each pep.xml takes exactly one *.mzML)",
                    path.display(),
                    i + 1,
                    n
                )))
            }
        };
        let out = default_out_path(&pep_xml);
        inputs.push(InputSpec { pep_xml, mzml, out });
    }

    if inputs.is_empty() {
        return Err(PtmProphetRsError::Usage(format!(
            "file-list {} names no pep.xml inputs",
            path.display()
        )));
    }
    Ok(inputs)
}

/// Default localized-output path for a batch input: `interact-X.pep.xml` → `interact-X.mod.pep.xml`
/// (the PTMProphet/oracle convention), preserving the original stem's case. Unknown suffixes get
/// `.mod.pep.xml` appended.
fn default_out_path(pep: &Path) -> PathBuf {
    let s = pep.to_string_lossy();
    let lower = s.to_ascii_lowercase();
    if lower.ends_with(".pep.xml") {
        PathBuf::from(format!("{}.mod.pep.xml", &s[..s.len() - ".pep.xml".len()]))
    } else if lower.ends_with(".pepxml") {
        PathBuf::from(format!("{}.mod.pepXML", &s[..s.len() - ".pepxml".len()]))
    } else {
        PathBuf::from(format!("{s}.mod.pep.xml"))
    }
}

/// Reconstruct PTMProphet's canonical `options="…"` header string from the resolved options.
///
/// This mirrors the string PTMProphet writes verbatim into its `<ptmprophet_summary>` header — NOT
/// PTMProphet-rs's own Unix CLI. PTMProphet echoes each of its arguments except the file names, in
/// the order given, each followed by one space: boolean flags upper-cased (`NOSTACK`), valued knobs
/// as `NAME=value`, and the mod string as given. PTMProphet-rs emits the tokens in the order of
/// FragPipe's PTMProphet command line
/// (`NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b <mods> MINPROB=0.5`), then a non-default
/// `MAXFRAGZ`, then the `MAXTHREADS=1 ` FragPipe appended to its per-file (single-thread)
/// PTMProphet invocation; PTMProphet-rs's own worker count does not affect results and is
/// deliberately not reflected here. Floats use Rust's shortest round-tripping form
/// (`10.0 → "10"`, `0.5 → "0.5"`).
///
/// `KEEPOLD` is always written: PTMProphet-rs never removes PTMProphet results already in the
/// input pepXML (the writer replaces only `<modification_info>`), which is what PTMProphet does
/// under `KEEPOLD`.
fn build_opt_string(cli: &Cli, frag_tol: f64, nions: &[IonKind], cions: &[IonKind]) -> String {
    let mut s = String::new();
    if cli.no_stack {
        s.push_str("NOSTACK ");
    }
    s.push_str("KEEPOLD ");
    if cli.static_tol {
        s.push_str("STATIC ");
    }
    s.push_str(&format!("FRAGPPMTOL={frag_tol} "));
    s.push_str(&format!("EM={} ", cli.em));
    // Ion series in PTMProphet's spelling: letters separated by commas (`NIONS=a,b`).
    let series = |kinds: &[IonKind]| -> String {
        kinds.iter().map(|k| k.to_char().to_string()).collect::<Vec<_>>().join(",")
    };
    s.push_str(&format!("NIONS={} ", series(nions)));
    // Non-default C-ion series (the default `y` is implicit and PTMProphet does not echo it).
    if cions != [IonKind::Y] {
        s.push_str(&format!("CIONS={} ", series(cions)));
    }
    // The mod string verbatim, as given to `--mods`.
    s.push_str(&format!("{} ", cli.mods));
    s.push_str(&format!("MINPROB={} ", cli.min_prob));
    // A non-default MAXFRAGZ, in PTMProphet's spelling. `--min-o` and `--mod-prec` accept only their
    // default 0 and are not written.
    if cli.max_frag_z != 0 {
        s.push_str(&format!("MAXFRAGZ={} ", cli.max_frag_z));
    }
    // FragPipe passed PTMProphet `MAXTHREADS=1` last; PTMProphet echoed it like every other token.
    s.push_str("MAXTHREADS=1 ");
    s
}

/// Parse the `--mods` string into [`ModType`]s.
///
/// Grammar: `AAs:shift[,AAs:shift]...`, read as PTMProphet's constructor reads it; `n`/`c` are
/// valid terminal "residues". Each type carries the `ptm_label` (`"AAs:shiftStr"`) used for the
/// `ptm=`/ROC attributes. A token whose shift *string* repeats an earlier token's joins that type,
/// as in PTMProphet: its residues are appended to the type's residues (duplicates kept) and put at
/// the front of its label, so `STY:79.966331,n:79.966331` gives the label `nSTY:79.966331`.
///
/// A repeat joins the earlier type in any token. PTMProphet merges it the same way but then, unless
/// the repeat is in the last token, reads on from that token's shift as one more type.
/// PTMProphet-rs does not follow it there: FragPipe builds `--mods` from MSFragger's variable-mod
/// table in table order, so one shift on two residue sets can sit in any token (e.g.
/// `M:15.9949,n:42.01056,K:42.01056,K:43.00581`).
///
/// Rejected rather than read differently from PTMProphet:
/// - a token that is not `AAs:shift`, such as PTMProphet's deprecated `STY,79.9663` form, a stray
///   `STY`, or an empty token;
/// - a shift that does not parse as a number;
/// - a neutral-loss suffix (`AAs:shift:nl1:...:nlN`), which PTMProphet-rs does not implement.
///
/// An empty (or all-blank) string means no modifications.
pub fn parse_mod_string(mod_string: &str) -> Result<Vec<ModType>> {
    let mut mods: Vec<ModType> = Vec::new();
    if mod_string.trim().is_empty() {
        return Ok(mods);
    }
    // Types that share the *string* form of the shift are one type.
    let mut seen: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
    for token in mod_string.split(',') {
        let parts: Vec<&str> = token.split(':').collect();
        if parts.len() < 2 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--mods token '{token}' is not of the form AAs:shift (e.g. STY:79.966331)"
            )));
        }
        if parts.len() > 2 {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "--mods token '{token}': neutral-loss suffixes in --mods are not implemented in PTMProphet-rs"
            )));
        }
        let (aa_str, shift_str) = (parts[0], parts[1]);
        let shift: f64 = shift_str.parse().map_err(|_| {
            PtmProphetRsError::UnsupportedOption(format!(
                "--mods token '{token}': '{shift_str}' is not a mass shift"
            ))
        })?;
        if let Some(&idx) = seen.get(shift_str) {
            mods[idx].aas.extend(aa_str.chars());
            mods[idx].ptm_label.insert_str(0, aa_str);
        } else {
            seen.insert(shift_str, mods.len());
            mods.push(ModType {
                aas: aa_str.chars().collect(),
                shift,
                ptm_label: format!("{aa_str}:{shift_str}"),
            });
        }
    }
    Ok(mods)
}

/// Parse an ion-series spec (`"b"`, `"a,b"`) into [`IonKind`]s, in input order.
///
/// PTMProphet splits `NIONS`/`CIONS` on commas and keeps the first character of each piece, so it
/// reads `ab` as `a` alone, and an empty value as its default series. To match it or fail, every
/// comma-separated piece must be exactly one lowercase ion letter of the requested terminus —
/// `a`/`b`/`c` when `nterm`, else `x`/`y`/`z` — and no letter may repeat (PTMProphet would match
/// a repeated series twice).
fn parse_ion_series(spec: &str, nterm: bool) -> Result<Vec<IonKind>> {
    let (flag, valid) = if nterm { ("--n-ions", "a/b/c") } else { ("--c-ions", "x/y/z") };
    let mut out: Vec<IonKind> = Vec::new();
    for piece in spec.split(',') {
        let mut chars = piece.chars();
        let kind = match (chars.next(), chars.next()) {
            (Some(c), None) if c.is_ascii_lowercase() => {
                IonKind::from_char(c).filter(|k| k.is_nterm() == nterm)
            }
            _ => None,
        };
        let kind = kind.ok_or_else(|| {
            PtmProphetRsError::UnsupportedOption(format!(
                "{flag} '{spec}': '{piece}' is not one {valid} ion letter (separate series with \
                 commas, e.g. {})",
                if nterm { "a,b" } else { "y,z" }
            ))
        })?;
        if out.contains(&kind) {
            return Err(PtmProphetRsError::UnsupportedOption(format!(
                "{flag} '{spec}': '{piece}' is given more than once"
            )));
        }
        out.push(kind);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ion_series_parses_and_validates() {
        // One letter per comma-separated piece, as PTMProphet reads NIONS/CIONS; order preserved.
        assert_eq!(parse_ion_series("b", true).unwrap(), vec![IonKind::B]);
        assert_eq!(parse_ion_series("a,b", true).unwrap(), vec![IonKind::A, IonKind::B]);
        assert_eq!(parse_ion_series("b,a", true).unwrap(), vec![IonKind::B, IonKind::A]);
        assert_eq!(parse_ion_series("y,z", false).unwrap(), vec![IonKind::Y, IonKind::Z]);
        assert_eq!(parse_ion_series("x", false).unwrap(), vec![IonKind::X]);
    }

    #[test]
    fn ion_series_rejects_bad_input() {
        assert!(parse_ion_series("q", true).is_err()); // unknown letter
        assert!(parse_ion_series("y", true).is_err()); // y is C-terminal, not an --n-ions ion
        assert!(parse_ion_series("b", false).is_err()); // b is N-terminal, not a --c-ions ion
        // Forms PTMProphet reads differently: `ab` is `a` alone, an empty value is its default.
        for spec in ["ab", "", "a,", "a, b", "B", "b,b"] {
            assert!(parse_ion_series(spec, true).is_err(), "{spec:?}");
        }
    }

    /// A repeated shift string joins the earlier type as in PTMProphet: residues appended (a
    /// duplicate kept), the new residues at the front of the label.
    #[test]
    fn mod_string_parses_and_merges_types() {
        let mods = parse_mod_string("STY:79.966331,M:15.9949,nS:79.966331").unwrap();
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].aas, vec!['S', 'T', 'Y', 'n', 'S']);
        assert_eq!(mods[0].shift, 79.966331);
        assert_eq!(mods[0].ptm_label, "nSSTY:79.966331");
        assert_eq!(mods[1].aas, vec!['M']);
        assert_eq!(mods[1].ptm_label, "M:15.9949");
        assert!(parse_mod_string("").unwrap().is_empty());
        assert!(parse_mod_string(" ").unwrap().is_empty());
    }

    /// FragPipe's `--mods` for an acetyl workflow with protein N-term acetyl enabled, in the order of
    /// MSFragger's variable-mod table: the repeated 42.01056 is not in the last token and still joins
    /// the n type, labeled as PTMProphet labels a merged type.
    #[test]
    fn mod_string_merges_a_repeat_in_any_token() {
        let mods =
            parse_mod_string("M:15.9949,n:42.01056,K:42.01056,K:43.00581,K:229.16293,N:0.984016")
                .unwrap();
        let labels: Vec<&str> = mods.iter().map(|m| m.ptm_label.as_str()).collect();
        assert_eq!(
            labels,
            ["M:15.9949", "Kn:42.01056", "K:43.00581", "K:229.16293", "N:0.984016"]
        );
        assert_eq!(mods[1].aas, vec!['n', 'K']);
        assert_eq!(mods[1].shift, 42.01056);
        assert!(mods[1].targets_nterm() && mods[1].can_modify('K'));
        assert_eq!(mods[2].aas, vec!['K']);
        assert_eq!(mods[2].shift, 43.00581);
        assert_eq!(mods[4].aas, vec!['N']);
        assert_eq!(mods[4].shift, 0.984016);
    }

    /// FragPipe's default `--mods` gives one type per token.
    #[test]
    fn mod_string_parses_fragpipe_default() {
        let mods = parse_mod_string("STY:79.966331,M:15.9949").unwrap();
        assert_eq!(mods.len(), 2);
        assert_eq!(mods[0].aas, vec!['S', 'T', 'Y']);
        assert_eq!(mods[0].shift, 79.966331);
        assert_eq!(mods[0].ptm_label, "STY:79.966331");
        assert_eq!(mods[1].aas, vec!['M']);
        assert_eq!(mods[1].shift, 15.9949);
        assert_eq!(mods[1].ptm_label, "M:15.9949");
    }

    #[test]
    fn mod_string_rejects_tokens_ptmprophet_reads_differently() {
        for (mods, expected) in [
            ("STY,79.9663", "'STY' is not of the form AAs:shift"),
            ("M:15.9949,STY", "'STY' is not of the form AAs:shift"),
            ("M:15.9949,,STY:79.966331", "'' is not of the form AAs:shift"),
            ("M:15.9949,", "'' is not of the form AAs:shift"),
            ("M:abc", "'abc' is not a mass shift"),
        ] {
            let err = parse_mod_string(mods).unwrap_err().to_string();
            assert!(err.contains(expected), "{mods}: {err}");
        }
    }

    #[test]
    fn mod_string_rejects_neutral_loss_suffix() {
        let err = parse_mod_string("M:15.9949,STY:79.966331:97.976896").unwrap_err().to_string();
        assert!(
            err.contains("neutral-loss suffixes in --mods are not implemented in PTMProphet-rs"),
            "unexpected error: {err}"
        );
        assert!(err.contains("STY:79.966331:97.976896"), "error does not name the token: {err}");
    }

    fn cli(args: &[&str]) -> std::result::Result<Cli, clap::Error> {
        let mut argv = vec!["ptmprophet-rs", "--filelist", "list.txt"];
        argv.extend_from_slice(args);
        Cli::try_parse_from(argv)
    }

    #[test]
    fn removed_flags_are_rejected() {
        for flag in ["--keep-old", "--no-update", "--no-min-o-factor"] {
            let err = cli(&["--mods", "M:15.9949", flag]).unwrap_err();
            assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument, "{flag}");
        }
        let err = cli(&["--mods", "M:15.9949", "--ppm-tol", "1"]).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument, "--ppm-tol");
    }

    fn from_cli_error(args: &[&str]) -> String {
        Config::from_cli(cli(args).unwrap()).unwrap_err().to_string()
    }

    /// Options PTMProphet-rs accepts but does not carry out as PTMProphet would are refused before
    /// any file is read.
    #[test]
    fn unimplemented_options_are_rejected() {
        let err = from_cli_error(&["--mods", "M:15.9949"]);
        assert!(err.contains("a run without --static"), "{err}");
        for (flag, value) in [("--min-o", "1"), ("--min-o=-1", ""), ("--mod-prec", "2")] {
            let mut args = vec!["--mods", "M:15.9949", "--static", flag];
            if !value.is_empty() {
                args.push(value);
            }
            let err = from_cli_error(&args);
            assert!(err.contains("is implemented"), "{flag} {value}: {err}");
        }
        let err = from_cli_error(&["--mods", "M:15.9949", "--static", "--frag-ppm-tol=-5"]);
        assert!(err.contains("--frag-ppm-tol -5"), "{err}");
        for z in ["1", "2"] {
            let err = from_cli_error(&["--mods", "M:15.9949", "--static", "--max-frag-z", z]);
            assert!(
                err.contains(&format!(
                    "--max-frag-z {z} (a positive --max-frag-z is not implemented in PTMProphet-rs)"
                )),
                "{err}"
            );
        }
    }

    /// `--max-frag-z` 0 or below passes the startup checks (the file-list is read next, so a
    /// missing one is the first error) and is echoed as PTMProphet's `MAXFRAGZ` when not 0.
    #[test]
    fn zero_or_negative_max_frag_z_is_accepted() {
        for z in ["--max-frag-z=0", "--max-frag-z=-1"] {
            let err = from_cli_error(&["--mods", "M:15.9949", "--static", z]);
            assert!(err.contains("list.txt"), "{z}: {err}");
        }
        let c = cli(&["--static", "--mods", "M:15.9949", "--max-frag-z=-1"]).unwrap();
        assert_eq!(
            build_opt_string(&c, 15.0, &[IonKind::B], &[IonKind::Y]),
            "KEEPOLD STATIC FRAGPPMTOL=15 EM=1 NIONS=b M:15.9949 MINPROB=0.5 MAXFRAGZ=-1 MAXTHREADS=1 "
        );
    }

    /// FragPipe's default PTMProphet-rs command line reproduces the `options` string PTMProphet
    /// printed for
    /// `NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b <mods> MINPROB=0.5 MAXTHREADS=1`.
    #[test]
    fn opt_string_matches_ptmprophet_for_fragpipe_default() {
        let c = cli(&[
            "--no-stack", "--static", "--frag-ppm-tol", "10", "--em", "1", "--n-ions", "b",
            "--mods", "STY:79.966331,M:15.9949", "--min-prob", "0.5",
        ])
        .unwrap();
        let nions = parse_ion_series(&c.n_ions, true).unwrap();
        let cions = parse_ion_series(&c.c_ions, false).unwrap();
        assert_eq!(
            build_opt_string(&c, 10.0, &nions, &cions),
            "NOSTACK KEEPOLD STATIC FRAGPPMTOL=10 EM=1 NIONS=b STY:79.966331,M:15.9949 MINPROB=0.5 MAXTHREADS=1 "
        );
    }

    /// Ion series are echoed comma-separated, as PTMProphet echoes NIONS/CIONS.
    #[test]
    fn opt_string_spells_ion_series_as_ptmprophet() {
        let c = cli(&[
            "--static", "--n-ions", "a,b", "--c-ions", "y,z", "--mods", "M:15.9949",
        ])
        .unwrap();
        let nions = parse_ion_series(&c.n_ions, true).unwrap();
        let cions = parse_ion_series(&c.c_ions, false).unwrap();
        assert_eq!(
            build_opt_string(&c, 15.0, &nions, &cions),
            "KEEPOLD STATIC FRAGPPMTOL=15 EM=1 NIONS=a,b CIONS=y,z M:15.9949 MINPROB=0.5 MAXTHREADS=1 "
        );
    }

    #[test]
    fn default_out_path_uses_mod_suffix() {
        assert_eq!(
            default_out_path(Path::new("/d/interact-f05.pep.xml")),
            PathBuf::from("/d/interact-f05.mod.pep.xml")
        );
        // Case of the stem is preserved; unknown suffixes get the full suffix appended.
        assert_eq!(
            default_out_path(Path::new("C:/d/Interact-F05.PEP.XML")),
            PathBuf::from("C:/d/Interact-F05.mod.pep.xml")
        );
        assert_eq!(default_out_path(Path::new("/d/run")), PathBuf::from("/d/run.mod.pep.xml"));
    }

    #[test]
    fn filelist_pairs_each_pepxml_with_one_mzml() {
        let dir = std::env::temp_dir();
        let path = dir.join("ptmprophet_rs_filelist_parse_test.txt");
        std::fs::write(
            &path,
            "# a batch of fractions\n\
             /data/interact-f01.pep.xml\t/raw/f01_calibrated.mzML\n\
             \n\
             /data/interact-f02.pep.xml\t/raw/f02_calibrated.mzML\n\
             /data/interact-f03.pep.xml /raw/f03.mzML\n",
        )
        .unwrap();

        let inputs = parse_filelist(&path).expect("parse");
        assert_eq!(inputs.len(), 3);
        assert_eq!(inputs[0].pep_xml, PathBuf::from("/data/interact-f01.pep.xml"));
        assert_eq!(inputs[0].mzml, PathBuf::from("/raw/f01_calibrated.mzML"));
        assert_eq!(inputs[0].out, PathBuf::from("/data/interact-f01.mod.pep.xml"));
        assert_eq!(inputs[1].mzml, PathBuf::from("/raw/f02_calibrated.mzML")); // tab-separated
        assert_eq!(inputs[2].mzml, PathBuf::from("/raw/f03.mzML")); // space-separated

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn filelist_rejects_bad_lines() {
        let dir = std::env::temp_dir();
        let bad = dir.join("ptmprophet_rs_filelist_bad_test.txt");
        std::fs::write(&bad, "/data/a.pep.xml /data/b.pep.xml\n").unwrap();
        assert!(parse_filelist(&bad).is_err()); // two pep.xml on one line
        std::fs::write(&bad, "/raw/only.mzML\n").unwrap();
        assert!(parse_filelist(&bad).is_err()); // no pep.xml
        std::fs::write(&bad, "/data/a.pep.xml\n").unwrap();
        assert!(parse_filelist(&bad).is_err()); // no mzML
        std::fs::write(&bad, "/data/a.pep.xml\t/x.mzML\t/y.mzML\n").unwrap();
        assert!(parse_filelist(&bad).is_err()); // more than one mzML
        std::fs::write(&bad, "# only comments\n\n").unwrap();
        assert!(parse_filelist(&bad).is_err()); // empty
        std::fs::remove_file(&bad).ok();
    }
}
