// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// Ported from batmass-io batmass-io-java/batmass-io/src/main/java/umich/ms/util/ScanIdNum.java.

//! Scan numbers from mzML native IDs, by the rules MSFragger uses for pepXML `start_scan`.
//!
//! MSFragger reads spectra through batmass-io, whose `ScanIdNum.getNumFromId` turns a spectrum's
//! native ID into its scan number. [`scan_number`] ports that function and [`is_positional`] ports
//! its `isPattern1`; their tests copy batmass-io's `ScanIdNumTest`.

use std::sync::LazyLock;

use regex::Regex;

/// Native-ID formats that hold more than one number. MSFragger numbers these spectra by their
/// position in the file, not by a number in the ID (batmass-io `ScanIdNum.pattern1`).
static POSITIONAL_IDS: LazyLock<Vec<Regex>> = LazyLock::new(|| {
    [
        r"controllerType=[1-9] controllerNumber=1 scan=[0-9]+", // Thermo, non-MS controllers
        r"controllerType=[0-9] controllerNumber=1 scan=[0-9]+ demux=[0-9]+", // Thermo, demultiplexed
        r"function=[0-9]+ process=[0-9]+ scan=[0-9]+", // Waters
        r"sample=[0-9]+ period=[0-9]+ cycle=[0-9]+ experiment=[0-9]+", // SCIEX WIFF
        r"declaration=[0-9]+ collection=[0-9]+ scan=[0-9]+", // Bruker U2
        r"source=[^ ]+ start=[0-9]+ end=[0-9]+", // Shimadzu Biotech
        r"jobRun=[0-9]+ spotLabel=[^ ]+ spectrum=[0-9]+", // SCIEX TOF/TOF
        r"frame=[0-9]+ scan=[0-9]+ frameType=[0-9]+", // UIMF
        r"frame=[0-9]+ scan=[0-9]+", // Bruker TDF
    ]
    .iter()
    .map(|p| whole_id(p))
    .collect()
});

/// Native-ID formats that hold the scan number, highest priority first; the number is capture
/// group 1 (batmass-io `ScanIdNum.pattern2`).
static NUMBERED_IDS: LazyLock<Vec<Regex>> = LazyLock::new(|| {
    [
        r"controllerType=0 controllerNumber=1 scan=([0-9]+)", // Thermo
        r"scan=([0-9]+)", // scan number only: Bruker/Agilent YEP, Bruker BAF, Shimadzu Biotech QTOF
        r"frame=([0-9]+)", // Bruker TSF; diaTracer pseudo-MS2 spectra
        r"([0-9]+)",
        r"spectrum=([0-9]+)",
        r"scanId=([0-9]+)", // Agilent MassHunter
        r"index=([0-9]+)",
    ]
    .iter()
    .map(|p| whole_id(p))
    .collect()
});

/// Compile `pattern` to match a whole native ID, as Java's `Matcher.matches()` does.
fn whole_id(pattern: &str) -> Regex {
    Regex::new(&format!("^(?:{pattern})$")).expect("native-ID pattern is valid")
}

/// The scan number MSFragger gives the spectrum with native ID `native_id`, which is the number it
/// writes to pepXML `start_scan`. `fallback` is the number for an ID that carries none MSFragger
/// reads; MSFragger passes the spectrum's 0-based mzML `index` plus 1.
///
/// After trimming whitespace, an ID that [`is_positional`] gets `fallback`, and any other ID gets
/// the number from the first format it matches whole, or `fallback` if it matches none. An empty
/// ID, on which batmass-io throws, gets `fallback`. Numbers are read as `i64`, so a number too large
/// for batmass-io's `int`, on which it throws, is returned as is; only one too large for `i64` gets
/// `fallback`.
pub(crate) fn scan_number(native_id: &str, fallback: i64) -> i64 {
    let id = native_id.trim();
    if is_positional(id) {
        return fallback;
    }
    NUMBERED_IDS
        .iter()
        .find_map(|re| re.captures(id))
        .map_or(fallback, |caps| caps[1].parse().unwrap_or(fallback))
}

/// Whether `native_id`, after trimming whitespace, is in a format that holds more than one number,
/// so MSFragger numbers its spectrum by position (batmass-io `ScanIdNum.isPattern1`). An empty ID
/// is not; batmass-io throws on it.
pub(crate) fn is_positional(native_id: &str) -> bool {
    let id = native_id.trim();
    POSITIONAL_IDS.iter().any(|re| re.is_match(id))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// batmass-io `ScanIdNumTest.getNumFromId`, case for case.
    #[test]
    fn get_num_from_id() {
        assert_eq!(23, scan_number("controllerType=0 controllerNumber=1 scan=2 demux=0", 23));
        assert_eq!(23, scan_number("sample=2 period=3 cycle=4 experiment=5", 23));
        assert_eq!(23, scan_number("frame=2 scan=3", 23));
        assert_eq!(3, scan_number("controllerType=1 controllerNumber=1 scan=23", 3));
        assert_eq!(23, scan_number("controllerType=0 controllerNumber=1 scan=23", 3));
        assert_eq!(2, scan_number("function=0 process=1 scan=23", 2));
        assert_eq!(3, scan_number("jobRun=0 spotLabel=asw spectrum=23", 3));
        assert_eq!(23, scan_number("23", 2));
        assert_eq!(23, scan_number("scan=23", 2));
        assert_eq!(23, scan_number("spectrum=23", 2));
        assert_eq!(23, scan_number("scanId=23", 3));
        assert_eq!(23, scan_number("index=23", 2));
        assert_eq!(2, scan_number("frame=2", 23));
    }

    /// batmass-io `ScanIdNumTest.isPattern1`, case for case.
    #[test]
    fn is_pattern1() {
        assert!(is_positional("controllerType=0 controllerNumber=1 scan=2 demux=0"));
        assert!(is_positional("sample=2 period=3 cycle=4 experiment=5"));
        assert!(is_positional("frame=2 scan=3"));
        assert!(is_positional("controllerType=1 controllerNumber=1 scan=23"));
        assert!(!is_positional("controllerType=0 controllerNumber=1 scan=23"));
        assert!(is_positional("function=0 process=1 scan=23"));
        assert!(is_positional("jobRun=0 spotLabel=asw spectrum=23"));
        assert!(!is_positional("23"));
        assert!(!is_positional("scan=23"));
        assert!(!is_positional("spectrum=23"));
        assert!(!is_positional("scanId=23"));
        assert!(!is_positional("index=23"));
        assert!(!is_positional("frame=2"));
    }

    /// Cases batmass-io's test does not cover.
    #[test]
    fn scan_number_more_formats() {
        // IDs that hold the scan number.
        assert_eq!(scan_number("frame=10090", 10090), 10090);
        assert_eq!(scan_number("frame=5", 1), 5, "frame=N gives N, not the fallback");
        assert_eq!(scan_number("index=3", 1), 3, "index=N gives N, without adding 1");
        assert_eq!(scan_number("  scan=7 ", 1), 7, "surrounding whitespace is ignored");

        // IDs with several numbers, or none MSFragger reads, give the fallback.
        assert_eq!(scan_number("controllerType=4 controllerNumber=1 scan=5", 100), 100);
        assert_eq!(scan_number("frame=12 scan=345 frameType=1", 7), 7);
        assert_eq!(scan_number("declaration=1 collection=2 scan=3", 4), 4);
        assert_eq!(scan_number("source=a.raw start=1 end=2", 5), 5);
        assert_eq!(scan_number("merged=0 frame=3 scanStart=1 scanEnd=900", 3), 3);
        assert_eq!(scan_number("nonsense", 5), 5);

        // Where batmass-io throws.
        assert_eq!(scan_number("", 1), 1, "an empty ID gives the fallback");
        assert_eq!(scan_number("scan=3000000000", 6), 3_000_000_000, "a number beyond Java's int is read as i64");
        assert_eq!(scan_number("scan=99999999999999999999", 6), 6, "a number too large for i64 gives the fallback");
        assert!(!is_positional(""));
        assert!(is_positional(" frame=2 scan=3 "), "is_positional trims, as isPattern1 does");
    }
}
