# PTMProphet-rs

PTMProphet-rs localizes variable post-translational modifications (PTMs) on peptide-spectrum
matches (PSMs). Given a `pep.xml` search result and its `mzML` run file, it scores, for each
ambiguous modification, which candidate residue carries it — using MS/MS fragment-ion evidence and a
1-D expectation–maximization (EM) mixture model — and writes a `mod.pep.xml` annotated with per-site
localization probabilities.

## Port of PTMProphet

PTMProphet-rs is a Rust port of PTMProphet, the PTM localization tool of the Trans-Proteomic
Pipeline (TPP), translated from the TPP 6.3.2 sources.
[Options and PTMProphet](#options-and-ptmprophet) maps each flag to its PTMProphet option, and
[Differences from PTMProphet](#differences-from-ptmprophet) lists where the two disagree.
PTMProphet is described in:

> Shteynberg DD, Deutsch EW, Campbell DS, Hoopmann MR, Kusebauch U, Lee D, Mendoza L, Midha MK,
> Sun Z, Whetton AD, Moritz RL. PTMProphet: Fast and Accurate Mass Modification Localization for the
> Trans-Proteomic Pipeline. *J Proteome Res* 2019;18(12):4262–4272.
> [doi:10.1021/acs.jproteome.9b00205](https://doi.org/10.1021/acs.jproteome.9b00205)

PTMProphet-rs is an independent Rust port of PTMProphet, part of the Trans-Proteomic Pipeline from
the Institute for Systems Biology; it is not affiliated with or endorsed by ISB or the TPP
developers.

## How it works

For every eligible PSM:

1. **Candidate enumeration** — find the residues each variable modification could occupy and
   enumerate the mass-consistent placements.
2. **Fragment matching** — generate the configured fragment-ion ladders for each placement (default
   `b`/`y`; selectable via `--n-ions`/`--c-ions`) and greedily match them to the observed peaks,
   claiming each peak at most once.
3. **Per-site scoring** — from the matched intensities, compute an O-score (intensity evidence) plus
   supporting m-/terminal scores for each candidate site.
4. **EM recalibration** — fit a weighted Gaussian kernel-density mixture model over all sites'
   O-scores and iterate to convergence, turning raw scores into calibrated probabilities.
5. **Output** — rewrite each PSM's `modification_info` to the localized sites and emit the per-site
   probability table, plus an ROC / error-rate summary in the file header.

## Usage

```
ptmprophet-rs --mods <MODS> --static [OPTIONS] --filelist <filelist_ptmloc.txt>
```

Key options: `--mods <spec>` (required, e.g. `M:15.9949,STY:79.96633`), `--static` (required),
`--filelist <file>` (required), `--frag-ppm-tol <ppm>`, `--em 1`, `--n-ions <ions>`,
`--c-ions <ions>`, `--min-prob <p>`, `--no-stack`, and `--threads`/`-j <n>` (default: all logical
CPUs). Run `ptmprophet-rs --help` for the full list.

### One run over many pep.xml files

`--filelist <file>` localizes every `interact-*.pep.xml` it names in a **single** PTMProphet-rs
run, but each file gets **its own** EM intensity model trained on **that file's** ambiguous-site
evidence — not one pooled model. PTMProphet is invoked once per pepXML, so its mixture model
(bandwidth + density curve) differs per file; training per file reproduces its per-site
probabilities. Each input is written to a sibling `interact-*.mod.pep.xml`.

The file-list has one input per non-blank line (`#` starts a comment line). Each line names one
pepXML and the one mzML its spectra are read from, separated by a **tab** (so paths may contain
spaces) or by whitespace. The mzML is required: PTMProphet-rs never takes the spectrum path from
inside the pepXML.

```text
# filelist_ptmloc.txt — five fractions, each pepXML with its mzML (tab- or space-separated)
/data/interact-f01.pep.xml    /raw/01CPTAC_..._f01_calibrated.mzML
/data/interact-f02.pep.xml    /raw/01CPTAC_..._f02_calibrated.mzML
/data/interact-f03.pep.xml    /raw/01CPTAC_..._f03_calibrated.mzML
/data/interact-f04.pep.xml    /raw/01CPTAC_..._f04_calibrated.mzML
/data/interact-f05.pep.xml    /raw/01CPTAC_..._f05_calibrated.mzML
```

```bash
./target/release/ptmprophet-rs --mods M:15.9949,STY:79.96633 --frag-ppm-tol 10 --em 1 --n-ions b \
    --min-prob 0.5 --no-stack --static --filelist filelist_ptmloc.txt
```

`--n-ions`/`--c-ions` select which backbone fragment series to match — N-terminal `a`/`b`/`c` and
C-terminal `x`/`y`/`z`, one letter per series, separated by commas as in PTMProphet (e.g.
`--n-ions a,b --c-ions y,z`). PTMProphet reads only the first letter of each comma-separated piece,
so PTMProphet-rs rejects `ab`, as it does an empty value or a repeated letter. Defaults are `b`/`y` (CID-style); for ETD spectra pass `--n-ions c --c-ions z`. (Series are not
auto-selected from the spectrum's activation method — unlike PTMProphet — so set them explicitly.)

```bash
cargo build --release
./target/release/ptmprophet-rs --mods M:15.9949,STY:79.96633 --frag-ppm-tol 10 --em 1 --n-ions b \
    --min-prob 0.5 --no-stack --static --filelist filelist_ptmloc.txt
```

## Options and PTMProphet

Every PTMProphet-rs flag and the `PTMProphetParser` option it stands for. The `options="…"`
attribute of the output's `<ptmprophet_summary>` echoes the PTMProphet spellings.

| PTMProphet-rs | PTMProphet | Notes |
|---|---|---|
| `--mods <MODS>` | `<AAs>:<shift>,…` | Required; PTMProphet defaults to `STY:79.9663,MW:15.9949`. A token may repeat an earlier shift string: it joins that type, and its residues go to the front of the type's label (`STY:79.966331,n:79.966331` → `nSTY:79.966331`), as in PTMProphet. |
| `--frag-ppm-tol <PPM>` | `FRAGPPMTOL=<PPM>` | Default 15, PTMProphet's default under `STATIC`. A negative value is rejected. |
| `--static` | `STATIC` | Required. Every PSM is matched with the one `--frag-ppm-tol`, as PTMProphet does under `STATIC`. |
| `--em 1` | `EM=1` | The only value accepted; others are rejected. PTMProphet defaults to `EM=2`. |
| `--min-prob <P>` | `MINPROB=<P>` | Default 0.5; PTMProphet defaults to 0.9. |
| `--n-ions <IONS>` | `NIONS=<IONS>` | `a`/`b`/`c`, comma-separated (`a,b`). Default `b`; PTMProphet defaults to `b,a`, or `c` for ETD spectra. |
| `--c-ions <IONS>` | `CIONS=<IONS>` | `x`/`y`/`z`, comma-separated (`y,z`). Default `y`; PTMProphet defaults to `y`, or `z` for ETD spectra. |
| `--no-stack` | `NOSTACK` | No two modifications on the same residue. |
| `--min-o 0` | `MINO=0` | The only value accepted; others are rejected. |
| `--max-frag-z <Z>` | `MAXFRAGZ=<Z>` | `0` (default) or below: the highest fragment charge is the precursor charge plus `Z`, at least 1. A positive value is rejected. Give a negative value as `--max-frag-z=-1`. |
| `--mod-prec 0` | `MODPREC=0` | The only value accepted; others are rejected. Modification masses are matched and written at integer precision. |
| `--threads`/`-j <N>` | `MAXTHREADS=<N>` | Default: all logical CPUs (PTMProphet: 1). Output does not depend on it, and `options` always shows FragPipe's `MAXTHREADS=1`. |
| `--filelist <FILE>` | `<input_file> [<output_file>]` | One pepXML and its mzML per line, each with its own EM model, as from one PTMProphet run per file. Output goes to `*.mod.pep.xml` beside the input. |
| *(always on)* | `KEEPOLD` | PTMProphet results already in the input pepXML are kept. |

PTMProphet options PTMProphet-rs does not implement. PTMProphet-rs rejects each of them:

| PTMProphet | PTMProphet-rs |
|---|---|
| `EM=0`, `EM=2`, `EM=3` | `--em` accepts only 1. |
| No `STATIC`: a fragment mass offset and tolerance estimated for each PSM | `--static` is required. |
| `MINO=` other than 0, with the MINO correction (`min_o_factors`, adjusted `mean_best_prob`) | `--min-o` accepts only 0. |
| `MODPREC=` other than 0 | `--mod-prec` accepts only 0. |
| `MAXFRAGZ=` above 0: the fragment charge capped at `Z` | `--max-frag-z` accepts only 0 or below. |
| `LABILITY`, `DIRECT`, `AUTODIRECT`, `IFRAGS`, `DENOISE`, `QUANTMODE`, `OSCOREMODE` | Not implemented. |
| `MASSDIFFMODE`, `MASSOFFSET=`, `EXCLUDEMASSDIFFMIN=`, `EXCLUDEMASSDIFFMAX=` | Not implemented. |
| `DALTONTOL=` | Not implemented; the fragment tolerance is in ppm only. |
| `VERBOSE`, `NORMAL`, `NOUPDATE`, `NOMINOFACTOR` | Not implemented. |
| `PPMTOL=` | Not implemented; the precursor check is fixed at 1 ppm. PTMProphet's is too: its check is a literal `ppm < 1`, and `PPMTOL` never reaches it. |
| Neutral-loss suffixes (`AAs:shift:nl1:…`) | `--mods` rejects them. |
| The deprecated `AAs,shift` form (`STY,79.9663`) | `--mods` rejects a token that is not `AAs:shift`, an empty token, and a shift that is not a number. |
| mzXML or mgf spectra | The file-list takes mzML only. |

## Differences from PTMProphet

- **x ions.** PTMProphet-rs computes x = y + CO − 2H. PTMProphet adds OH where the carbon belongs,
  so its x ions are about 5 Da heavy (`src/model/peptide.rs`).
- **Ion series.** PTMProphet switches to `c`/`z` ions for spectra whose activation method is ETD.
  PTMProphet-rs never reads the activation method; set `--n-ions`/`--c-ions` (e.g.
  `--n-ions c --c-ions z`).
- **Repeated masses in `--mods`.** PTMProphet misreads a repeated shift string in any `--mods`
  token but the last: it merges the token into the earlier type, then reads on from that token's
  shift as one more type. With `M:15.9949,n:42.01056,K:42.01056,K:43.00581,…`, which FragPipe passes
  when protein N-term acetyl is on in an acetyl workflow, PTMProphet adds a type with the residues
  `42.01056,K` and the shift 43.00581. PTMProphet-rs merges a repeat into the earlier type wherever
  it appears and labels the merged type as PTMProphet does: `M:15.9949`, `Kn:42.01056`,
  `K:43.00581`, ….
- **Fragment charge.** With `--max-frag-z` 0 or below, both match fragment charges 1 to the
  precursor charge plus `Z`, and charge 1 alone when that sum is 0. PTMProphet computes the sum in
  unsigned arithmetic (`initMaxCharge` in `PTMProphetMpx.h`), so the two differ in two cases:
  - `Z` below minus the precursor charge (e.g. `--max-frag-z=-3` for a 2+ precursor): the sum wraps
    around, and PTMProphet matches every charge up to the precursor charge. PTMProphet-rs matches
    charge 1.
  - `assumed_charge="0"`: PTMProphet matches no fragment ions. PTMProphet-rs matches charge 1.
    PTMProphet-rs also reads a missing `assumed_charge` as 0; PTMProphet does not handle a missing
    one.
- **Existing results.** PTMProphet-rs always keeps PTMProphet results already in the input pepXML,
  as PTMProphet does under `KEEPOLD`, so `options` always shows `KEEPOLD`.
- **Version attribute.** `<ptmprophet_summary version="…">` reads
  `PTMProphet-rs <version> (Rust port of TPP PTMProphet)`, not a TPP build. The element names and
  `analysis="ptmprophet"` stay, because Philosopher and other pepXML readers look for them.

## Performance

- **mzML** is read with [`mzdata`](https://crates.io/crates/mzdata) (configured for a pure-Rust
  build — no C toolchain needed).
- The pipeline is **multi-threaded**: spectra are decoded, PSMs are matched, the EM density fit is
  computed, and output blocks are built across all logical CPUs (`std::thread::scope`, no external
  runtime). Results are merged in input order, so the output is **identical regardless of thread
  count**. Set the worker count with `--threads`/`-j <n>` (or `PTMPROPHET_RS_THREADS=<n>`, which
  takes precedence); set `PTMPROPHET_RS_TIMING=1` for per-phase timings.
- [`mimalloc`](https://crates.io/crates/mimalloc) is the global allocator.

## Layout

```
src/
  config.rs  pipeline.rs  threads.rs   CLI parsing, orchestration, data-parallel helpers
  model/     masses, modspec, peptide  monoisotopic masses + theoretical fragment ions
  pepxml/    tag, index, psm, writer, roc   pep.xml parse, localized output, ROC summary
  spectrum/  mzml                      mzML reading
  mpx/       combinatorics             per-PSM localization engine + peak matching
  em/        kdmodel, density, train   KDE mixture model + EM training
```

## Build & test

```bash
cargo build --release
cargo test
```

## License

PTMProphet-rs is licensed under the GNU Lesser General Public License, version 3 or later
(`LGPL-3.0-or-later`). It is derived from TPP source files licensed under LGPL version 2.1 "or (at
your option) any later version", and that grant lets the derived work be distributed under version 3.
Version 3 is used because PTMProphet-rs links Apache-2.0 crates (`mzdata`, `mzpeaks`), which are
compatible with version 3.

- [`COPYING.LESSER`](COPYING.LESSER): the LGPL version 3.
- [`COPYING`](COPYING): the GPL version 3, which the LGPL supplements.
- [`NOTICE`](NOTICE): the TPP files PTMProphet-rs derives from, their copyright holders and the
  upstream license notice. Each source file's header names its own upstream files and holders.

Release zips also carry `THIRD-PARTY-LICENSES.txt`, the licenses of the Rust crates compiled into
the program, and `licenses/Rust-standard-library.html`, the licenses of the Rust standard library.
The build writes `THIRD-PARTY-LICENSES.txt` with cargo-about from [`about.toml`](about.toml) and
[`about.hbs`](about.hbs).
