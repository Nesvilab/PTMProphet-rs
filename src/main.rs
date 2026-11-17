// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan

//! PTMProphet-rs binary entry point.
//!
//! Parses a standard Unix-style command line (see [`ptmprophet_rs::config::Cli`]) into a
//! [`ptmprophet_rs::Config`] and runs the localization pipeline. Input is a `--filelist` naming the
//! `interact-*.pep.xml` files to localize; each file trains its own EM model on its own evidence
//! (matching PTMProphet's per-file invocation) and is written to a sibling `interact-*.mod.pep.xml`:
//! ```text
//! ptmprophet-rs --mods M:15.9949,STY:79.96633 --frag-ppm-tol 10 --em 1 --n-ions b \
//!         --min-prob 0.5 --no-stack --static --filelist filelist_ptmloc.txt
//! ```

use std::io::Write;
use std::process::ExitCode;

use clap::Parser;

use ptmprophet_rs::config::{Cli, Config};
use ptmprophet_rs::pipeline;

// mimalloc as the global allocator — faster than the system allocator for PTMProphet-rs's
// allocation pattern (many short-lived per-PSM/per-spectrum buffers churned across worker threads).
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() -> ExitCode {
    // Wall-clock start, captured before any work so the closing message reflects the whole run.
    let start = std::time::Instant::now();

    // Logging mirrors FragSpecLib: the `log` facade over the `env_logger` backend, formatted with a
    // local timestamp + colored level and no module target, filtered at runtime via RUST_LOG (default
    // `info`). All records go to stderr.
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format(|buf, record| {
            let ts = chrono::Local::now().format("%Y-%m-%d %H:%M:%S");
            let level_style = buf.default_level_style(record.level());
            writeln!(
                buf,
                "[{ts} {level_style}{}{level_style:#}] {}",
                record.level(),
                record.args(),
            )
        })
        .init();

    // clap handles --help / --version and prints usage errors itself (exiting non-zero).
    let cli = Cli::parse();

    let config = match Config::from_cli(cli) {
        Ok(c) => c,
        Err(e) => {
            // Startup/config failure, before the pipeline: print the full error chain straight to
            // stderr (eprintln, per the logging policy) and exit non-zero.
            eprintln!("error: {:#}", e);
            return ExitCode::FAILURE;
        }
    };

    // Wire an explicit --threads into the worker pool (PTMPROPHET_RS_THREADS still overrides).
    if let Some(n) = config.threads {
        ptmprophet_rs::threads::set_thread_override(n);
    }

    match pipeline::run(&config) {
        Ok(()) => {
            // Run time in seconds as the final message on a successful run.
            log::info!("Run time: {:.2} s", start.elapsed().as_secs_f64());
            ExitCode::SUCCESS
        }
        Err(e) => {
            // Final error exit: surface the whole error chain on stderr (eprintln, per the logging
            // policy — the pipeline has already returned) and exit non-zero.
            eprintln!("error: {:#}", e);
            ExitCode::FAILURE
        }
    }
}
