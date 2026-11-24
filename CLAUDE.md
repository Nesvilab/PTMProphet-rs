# PTMProphet-rs

## Licensing (STRICT)

PTMProphet-rs is **LGPL-3.0-or-later**: it is translated from the Trans-Proteomic Pipeline's
PTMProphet sources, which are LGPL-2.1-or-later. `COPYING`, `COPYING.LESSER` and `NOTICE` at the
repo root carry the terms and the upstream attribution.

**Dependencies.** Every direct *and* transitive dependency must be compatible with LGPL-3.0/GPL-3.0:

- **Allowed:** `MIT`, `Apache-2.0`, `BSD-2-Clause`, `BSD-3-Clause`, `ISC`, `Zlib`, `MPL-2.0`,
  `Unicode-3.0`, `Unicode-DFS-2016`, `Unlicense`, `CC0-1.0`, `0BSD`, and LGPL/GPL version 3 or any
  "or later" version that reaches 3 (`LGPL-3.0-*`, `GPL-3.0-*`, `LGPL-2.1-or-later`,
  `GPL-2.0-or-later`). A dual license qualifies when one of its options is allowed.
- **Prohibited:** `GPL-2.0-only`, `LGPL-2.1-only`, and GPL-incompatible licenses (`EPL-1.0`, `CDDL`,
  `SSPL`, proprietary).
- **Before adding any crate**, check its license and its dependency tree, e.g. `cargo tree` plus the
  `license` field on crates.io, or `cargo deny check licenses`.
- **`about.toml` lists the licenses the release build accepts.** Its `accepted` holds only the licenses
  the tree uses now. A new crate under another allowed license needs that license added there, and a
  crate whose license files cargo-about cannot read needs a `clarify` entry, or the release build
  stops. Check with `cargo about generate --locked --fail about.hbs` (cargo-about 0.9.2): no
  `MISSING:` line and nothing on stderr.

**TPP-derived code.**

- Every `.rs` file starts with `// SPDX-License-Identifier: LGPL-3.0-or-later` and
  `// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan`; give a new file the
  same two lines.
- A file translated from TPP also carries one `SPDX-FileCopyrightText` line per upstream holder and a
  `// Translated into Rust from TPP 6.3.2 src/…` line, and `NOTICE` lists its upstream files. Keep
  those lines and entries. When TPP-derived code moves to another file, or more code is translated
  from TPP, carry the holder lines along and update `NOTICE`.
- Never copy PTMProphet-rs code into the lab's proprietary repositories (FragDIA, FragViz,
  FragSpecLib, FragCast, FragNovo, ...); they cannot take LGPL code.

## Release builds (FragPipe-tools)

PTMProphet-rs (GitHub: `Nesvilab/PTMProphet-rs`) has no CI of its own. Releases are built by
`../FragPipe-tools/.github/workflows/build-ptmprophet-rs.yml`, which starts in two ways:

- **Daily, on its own.** FragPipe-tools' `start-builds.yml` runs at 05:07 UTC and watches
  `ptmprophet-rs master`. When the `Cargo.toml` version at master's last commit has no release
  `ptmprophet-rs-<version>` yet, it starts `build-ptmprophet-rs.yml` at master's newest *bump
  commit* of that version: a commit whose subject starts with `Bump to `, `bump to `, `Release ` or
  `release ` and that changes `Cargo.toml` (e.g. `Bump to 0.2.0`). So a bump commit on master
  publishes that version at the next daily run (00:07 or 01:07 US Eastern), without the commits
  after it. The exact rules are under "Automatic builds" in `../FragPipe-tools/README.md`.
- **By hand**, at any tag, branch or commit:
  `gh workflow run build-ptmprophet-rs.yml -R Nesvilab/FragPipe-tools -f ref=<tag-or-commit>`.

Either way:

- It reads the version from `Cargo.toml` and fails when the release `ptmprophet-rs-<version>`
  already exists, so a new build needs a version bump.
- It builds with `--locked`, so `Cargo.lock` must be committed and current.
- It publishes `ptmprophet-rs-<version>-windows.zip` (`windows/ptmprophet-rs.exe`) and
  `ptmprophet-rs-<version>-linux.zip` (`linux/ptmprophet-rs`). **The binary keeps its plain
  `ptmprophet-rs[.exe]` name — never add the version to it.**
- It copies `COPYING`, `COPYING.LESSER` and `NOTICE` from the repo root to the root of each zip, and
  stops before building a commit that lacks one of them. Keep them at the root.
- It writes `THIRD-PARTY-LICENSES.txt` (cargo-about, from `about.toml` and `about.hbs`) and
  `licenses/Rust-standard-library.html` at the root of each zip, and stops at a commit without
  `about.toml` or `about.hbs` at the root, or whose `about.toml` `targets` lacks the target it builds
  (`x86_64-pc-windows-msvc`, `x86_64-unknown-linux-gnu`).

## Linux builds: glibc 2.28 floor

The Linux binary must not reference any glibc symbol newer than **2.28** (covers RHEL/Rocky 8 and
newer — many users are on old HPC clusters). Exceeding it means the binary won't start at all:
`version 'GLIBC_2.xx' not found`.

The floor is pinned by the version-suffixed target in the Linux matrix row of
`build-ptmprophet-rs.yml`, which is cargo-zigbuild's only interface for this:

```
cargo zigbuild --release --locked --target x86_64-unknown-linux-gnu.2.28 --bin ptmprophet-rs
```

- **Never switch to musl.**
- **Declare the version once** — the `target` field of the Linux matrix row. The `check-glibc`
  action (`../FragPipe-tools/.github/actions/check-glibc`) parses it back out; don't restate it.
- **`target` vs `dir` are intentionally separate.** cargo-zigbuild strips the suffix before calling
  cargo, so output lands in the unsuffixed `target/x86_64-unknown-linux-gnu/release/`. `dir` drives
  the rustup target and the binary path. Merging them breaks the build.
- `check-glibc` runs `objdump -T` on the binary and fails if any `GLIBC_*` symbol exceeds the floor
  or the binary needs `GLIBC_ABI_DT_RELR`. It runs before the upload, so a bad binary is never
  released.
- Plain `cargo build` on Linux loses the pin. Windows is unaffected (MSVC, plain `cargo build`).

To change the floor, edit that one suffix in `build-ptmprophet-rs.yml`.

## Windows builds: static C runtime

`.cargo/config.toml` links the C runtime statically for `x86_64-pc-windows-msvc` (`+crt-static`), so
the exe imports Windows system DLLs only and starts on a Windows without the Visual C++
redistributable. The release build checks the imports with FragPipe-tools' `check-windows-dlls`
action. A `RUSTFLAGS` environment variable replaces the config's flags, so a build that sets one
must add `-C target-feature=+crt-static` to it.

## Logging & message printing

All runtime messages go through the **`log` crate**, routed by the `env_logger` set up once in
`src/main.rs`. Follow these rules:

- **Use `log::{info,warn,error,debug,trace}!` for every runtime message.** Never use `println!`
  or `print!` (stdout is left clean), and do not call `eprintln!` outside the two places below.
- **`eprintln!` is reserved for `main.rs` only** — fatal startup/exit errors that occur outside
  the pipeline (CLI parse failure, thread-pool setup, final error exit). Print them as
  `eprintln!("error: {:#}", e)` so the full `anyhow` chain is shown. Tests may also use
  `eprintln!` for diagnostics.
- **Print full error chains with the `{:#}` alternate format** (`anyhow`), e.g.
  `log::error!("convert failed for {}: {:#}", path.display(), e)`.

### Level conventions

| Level   | Use for                                                                          |
|---------|----------------------------------------------------------------------------------|
| `info`  | Progress milestones and summary counts (parsing, annotating, "Wrote N …", timing). |
| `warn`  | Recoverable problems or skipped work (no PSMs matched a file, figure failed).     |
| `error` | Per-item failures inside a loop that do **not** abort the whole run.              |
| `debug` / `trace` | Deeper diagnostics, surfaced on demand via `RUST_LOG`.                  |

### Style

- Lead with a capitalized word: `"Parsing UniMod database"`, `"Wrote {n} transitions to {…}"`.
- Include counts **with units**: `"{} PSMs annotated"`, `"{} transition rows"`, `"{:.1}s"`.
- Render paths with `Path::display()`, never `{:?}` on the raw `Path`.
- Indent sub-items under a heading line with two leading spaces (e.g. `"  {}: {} PSMs used in fit"`).

### Format & verbosity

- Output goes to **stderr** as `[YYYY-MM-DD HH:MM:SS LEVEL] message`, with the level colorized.
  The format is defined once in `main.rs` — do not re-implement it elsewhere.
- Verbosity is controlled by the **`RUST_LOG`** environment variable and defaults to `info`
  (e.g. `RUST_LOG=debug` for more detail, `RUST_LOG=warn` for less).

## Release profile

`Cargo.toml` pins an aggressively optimized release profile. Keep it as-is unless there's a
measured reason to change it:

```toml
[profile.release]
opt-level = 3        # full optimization
lto = "fat"          # whole-program link-time optimization
codegen-units = 1    # one unit for best optimization (slower compile)
panic = "abort"      # no unwinding — smaller, faster; panics terminate the process
strip = true         # strip symbols from the binary
```

Note: `panic = "abort"` means panics cannot be caught/unwound — recoverable error paths must use
`Result`/`anyhow`, not `catch_unwind`.
