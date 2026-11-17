// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// SPDX-FileCopyrightText: 2006 Henry Lam
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetMpx.cpp,
// src/Validation/PTMProphetParser/PTMProphetMpx.h and src/Search/SpectraST/SpectraSTPeakList.hpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! Per-PSM worker. Identifies candidate sites, enumerates mod placements, matches b/y fragments to
//! the spectrum, computes per-site O/M/direct/terminal scores, and — once an EM model is available —
//! the final per-site probabilities.
//!
//! Candidate sites are indexed by **0-based stripped residue index** and output `position = idx + 1`
//! (e.g. TMT-n-term peptide `ALRASESGI` → positions 5,7). A terminus never becomes a candidate site:
//! a mod type on `n`/`c` (e.g. `n:42.0106`) has one possible site, so when the matching terminus
//! carries its shift the type is flagged [`TypeResult::terminus`] and reported at probability 1,
//! unscored and without EM evidence.
//!
//! Masses use the **exact** static-mod masses from the search-summary tables plus the command-line
//! variable shifts — never the integer-rounded brackets, which are used only to *identify* which
//! input positions carry which mod.

pub mod combinatorics;

use std::collections::BTreeMap;

use crate::config::Config;
use crate::em::kdmodel::KdModel;
use crate::em::train::EmEvidence;
use crate::model::masses::{self, IonKind};
use crate::model::modspec::ModType;
use crate::model::peptide::Peptide;
use crate::pepxml::index::ModTables;
use crate::spectrum::mzml::Spectrum;

/// All computed scores for one candidate site. Keyed in [`TypeResult::scores`] by 0-based residue
/// index; the output `position` is that key + 1.
#[derive(Debug, Clone, Default)]
pub struct SiteScore {
    pub oscore: f64,
    pub mscore: f64,
    pub direct_oscore: f64,
    pub direct_mscore: f64,
    pub cterm_score: f64,
    pub nterm_score: f64,
    pub shift: char, // 'n' | 'c'
    pub prob: f64,   // final localization probability (set by the EM)
}

/// Per-mod-type localization result.
#[derive(Debug, Clone, Default)]
pub struct TypeResult {
    /// Candidate residue indices (0-based stripped), ascending.
    pub sites: Vec<usize>,
    /// Number of mods of this type present on the PSM.
    pub nmods: usize,
    /// Per candidate site (key = 0-based residue index) → scores.
    pub scores: BTreeMap<usize, SiteScore>,
    /// Set to `Some('n')`/`Some('c')` when this is a **terminal** mod type (e.g. n-term acetyl) and
    /// the matching terminus carries it on this PSM. A terminus has exactly one candidate site, so the
    /// mod is always fully localized at probability 1 — it is reported directly (no residue scoring,
    /// no EM evidence), via `<mod_terminal_probability>` rather than `<mod_aminoacid_probability>`.
    pub terminus: Option<char>,
}

/// One enumerated placement and the fragment evidence it produced.
struct Placement {
    /// Modified residue indices per type.
    mods_by_type: Vec<Vec<usize>>,
    /// Claimed observed peaks → intensity (each peak claimed by at most one fragment).
    matched: BTreeMap<u64, f64>,
    total: f64,
    count: i32,
    direct_total: f64,
    direct_count: i32,
    nterm_run: i32,
    cterm_run: i32,
}

/// Mutable accumulators threaded through one placement's ion-ladder matching. Shared across every
/// configured ion series so a peak is claimed at most once per placement.
#[derive(Default)]
struct MatchAcc {
    /// Peaks already claimed, keyed by m/z bit pattern.
    claimed: std::collections::HashSet<u64>,
    /// Claimed observed peaks → intensity.
    matched: BTreeMap<u64, f64>,
    for_evidence: f64,
    for_count: i32,
    direct_evidence: f64,
    direct_count: i32,
}

/// Per-PSM localization engine.
pub struct Mpx<'a> {
    residues: Vec<char>,
    /// Base peptide: stripped sequence + all *static* mods (residue + termini) applied; variable
    /// mods are added per placement.
    base: Peptide,
    charge: i32,
    calc_neut_mass: f64,
    mod_types: &'a [ModType],
    /// Per type: candidate sites + nmods + (filled by `process`) scores.
    pub types: Vec<TypeResult>,
    frag_tol: f64,
    nostack: bool,
    min_o: i32,
    max_frag_z: i32,
    /// N-terminal series to generate (`--n-ions`, a/b/c), in user-specified order.
    nions: Vec<IonKind>,
    /// C-terminal series to generate (`--c-ions`, x/y/z), in user-specified order.
    cions: Vec<IonKind>,
}

impl<'a> Mpx<'a> {
    /// Build an Mpx from the parsed `modification_info`.
    ///
    /// * `pep_seq` — stripped peptide.
    /// * `pep_mods` — 1-based position → absolute residue+mod mass; key `0` = N-term, `len+1` = C-term.
    ///   Full-precision masses from the input.
    /// * `mods`/`charge`/`calc_neut_mass` — config + PSM precursor info.
    pub fn build(
        config: &Config,
        mod_types: &'a [ModType],
        _tables: &ModTables,
        pep_seq: &str,
        pep_mods: &BTreeMap<usize, f64>,
        charge: i32,
        calc_neut_mass: f64,
    ) -> Mpx<'a> {
        let residues: Vec<char> = pep_seq.chars().map(|c| c.to_ascii_uppercase()).collect();
        let len = residues.len();

        // Base peptide from the bare sequence (residue cache built by parse); static mods added below.
        let mut base = Peptide::parse(&pep_seq.to_uppercase(), charge, config.mod_prec)
            .unwrap_or_default();

        // Which (type, residue-index) positions are variable mods present in the input.
        let mut input_modded: Vec<Vec<bool>> = mod_types.iter().map(|_| vec![false; len]).collect();

        for (&pos, &abs_mass) in pep_mods.iter() {
            if pos == 0 {
                base.nterm_mod = Some(abs_mass - masses::GROUP_N_TERM);
                continue;
            }
            if pos == len + 1 {
                base.cterm_mod = Some(abs_mass - masses::GROUP_C_TERM);
                continue;
            }
            let idx = pos - 1;
            if idx >= len {
                continue;
            }
            let aa = residues[idx];
            let res_mono = masses::residue_mono_mass(aa).unwrap_or(0.0);
            // Variable mod of some type? (integer-rounded match, mod_prec=0.)
            let mut classified = false;
            for (t, mt) in mod_types.iter().enumerate() {
                if mt.can_modify(aa)
                    && (abs_mass).round() == (res_mono + mt.shift).round()
                {
                    input_modded[t][idx] = true;
                    classified = true;
                    break;
                }
            }
            if !classified {
                // Static residue mod (e.g. C+57, K+TMT): apply to base, never moved.
                base.set_mod(idx, abs_mass - res_mono);
            }
        }

        // Candidate sites + nmods per type. Terminal mod types (`n`/`c`) carry no *residue* sites; if
        // the matching terminus holds this type's shift (kept on `base.{n,c}term_mod` for fragment
        // masses), the type is flagged `terminus` — a single, unambiguous, fully-localized site.
        let mut types = Vec::with_capacity(mod_types.len());
        for (t, mt) in mod_types.iter().enumerate() {
            let mut sites = Vec::new();
            for (i, &aa) in residues.iter().enumerate() {
                if mt.can_modify(aa) {
                    sites.push(i);
                }
            }
            let nmods = (0..len).filter(|&i| input_modded[t][i]).count();
            let mut terminus = None;
            if mt.targets_nterm() {
                if let Some(nm) = base.nterm_mod {
                    if (nm - mt.shift).abs() < 0.05 {
                        terminus = Some('n');
                    }
                }
            }
            if terminus.is_none() && mt.targets_cterm() {
                if let Some(cm) = base.cterm_mod {
                    if (cm - mt.shift).abs() < 0.05 {
                        terminus = Some('c');
                    }
                }
            }
            types.push(TypeResult { sites, nmods, scores: BTreeMap::new(), terminus });
        }

        Mpx {
            residues,
            base,
            charge,
            calc_neut_mass,
            mod_types,
            types,
            frag_tol: config.frag_tol,
            nostack: config.nostack,
            min_o: config.min_o,
            max_frag_z: config.max_frag_z,
            nions: config.nions.clone(),
            cions: config.cions.clone(),
        }
    }

    /// True if any type has a localization to report: a residue mod present (`nmods > 0`) or a
    /// terminal mod on this PSM.
    pub fn has_work(&self) -> bool {
        self.types.iter().any(|t| t.nmods > 0 || t.terminus.is_some())
    }

    /// Harvest the per-site O-scores of every **ambiguous** type (nmods>0 and sites>nmods) into the
    /// EM evidence. Fully-localized and absent types contribute nothing (their output prob string is
    /// empty, so they are skipped). Call after [`Mpx::process`].
    pub fn harvest_evidence(&self, ev: &mut EmEvidence) {
        let len = self.residues.len();
        for t in &self.types {
            if t.nmods == 0 || t.sites.len() <= t.nmods {
                continue;
            }
            // The harvest loops over EVERY residue index (0..len) for an ambiguous type, not just
            // candidate sites — non-candidate positions contribute O-score 0 (their reset value).
            // Those zeros anchor the negative class of the EM density at O-score 0.
            for idx in 0..len {
                let oscore = t.scores.get(&idx).map(|s| s.oscore).unwrap_or(0.0);
                ev.push(oscore);
            }
            ev.total_mods += t.nmods;
            ev.total_sites += t.sites.len();
        }
    }

    /// Information gain: returns `(norm_info_gain, localized_mods)` for a type from its per-site
    /// probabilities. `S` = #candidate sites, `N` = #mods. raw entropy `H = -Σ p·ln p`; normalized
    /// `H' = H/ln(S/N)` (0 if S==N); `localized_mods = N − H'`; `norm_info_gain = max(0, 1 − H'/N)`.
    pub fn info_gain(&self, t: usize) -> (f64, f64) {
        let tr = &self.types[t];
        let s = tr.sites.len() as f64;
        let n = tr.nmods as f64;
        let mut h = 0.0;
        for sc in tr.scores.values() {
            let p = sc.prob;
            if p > 0.0 {
                h += (0.0 - p) * p.ln();
            }
        }
        let norm_entropy = if (s - n).abs() < 1e-9 || n == 0.0 { 0.0 } else { h / (s / n).ln() };
        let localized = n - norm_entropy;
        let nig = (1.0 - norm_entropy / n).max(0.0);
        (nig, localized)
    }

    /// Compute the final per-site localization probabilities from the trained intensity model
    /// (EM=1 → only the O model contributes). Fills `SiteScore::prob`.
    ///
    /// Per type: skip if nmods==0; if nmods==sites → all prob 1; else per-site Bayesian posterior
    /// with the flat prior `Pv = nmods/|sites|`, then normalize so the sites sum to `nmods`.
    pub fn compute_probs(&mut self, model: &KdModel) {
        for t in 0..self.types.len() {
            let sites = self.types[t].sites.clone();
            let nmods = self.types[t].nmods;
            if sites.is_empty() || nmods == 0 {
                continue;
            }
            if nmods == sites.len() {
                for &pos in &sites {
                    if let Some(s) = self.types[t].scores.get_mut(&pos) {
                        s.prob = 1.0;
                    }
                }
                continue;
            }
            let pv = nmods as f64 / sites.len() as f64;
            let mut prs: Vec<(usize, f64)> = Vec::with_capacity(sites.len());
            let mut sum = 0.0;
            for &pos in &sites {
                let o = self.types[t].scores.get(&pos).map(|s| s.oscore).unwrap_or(0.0);
                let pos_d = model.get_pos_prob(o);
                let neg_d = model.get_neg_prob(o);
                let denom = pv * pos_d + (1.0 - pv) * neg_d;
                let pr = if denom > 0.0 { pv * pos_d / denom } else { pv };
                prs.push((pos, pr));
                sum += pr;
            }
            for (pos, pr) in prs {
                let mut prob = if sum > 0.0 { nmods as f64 * pr / sum } else { 0.0 };
                if prob > 1.0 {
                    prob = 1.0;
                }
                if let Some(s) = self.types[t].scores.get_mut(&pos) {
                    s.prob = prob;
                }
            }
        }
    }

    fn len(&self) -> usize {
        self.residues.len()
    }

    /// The stripped residue characters (read-only) — used by the output writer.
    pub fn residues_view(&self) -> &[char] {
        &self.residues
    }

    /// The localized residue positions per type after [`Mpx::compute_probs`]: the top-`nmods`
    /// candidate sites by probability, `--no-stack` ⇒ at most one mod per residue across types.
    /// Positions are 0-based residue indices (sorted ascending per type).
    pub fn localized_positions(&self) -> Vec<Vec<usize>> {
        let mut taken: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut out = Vec::with_capacity(self.types.len());
        for t in &self.types {
            let mut ranked: Vec<(usize, f64)> =
                t.scores.iter().map(|(&p, s)| (p, s.prob)).collect();
            // descending probability; stable tie-break by ascending position
            ranked.sort_by(|a, b| {
                b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal).then(a.0.cmp(&b.0))
            });
            let mut chosen = Vec::new();
            for (pos, _) in ranked {
                if chosen.len() >= t.nmods {
                    break;
                }
                if self.nostack && taken.contains(&pos) {
                    continue;
                }
                chosen.push(pos);
                taken.insert(pos);
            }
            chosen.sort_unstable();
            out.push(chosen);
        }
        out
    }

    fn init_max_charge(&self) -> u32 {
        max_fragment_charge(self.charge, self.max_frag_z)
    }

    /// Build the modified peptide for a placement (base static mods + variable shifts).
    fn placement_peptide(&self, mods_by_type: &[Vec<usize>]) -> Peptide {
        let mut p = self.base.clone();
        for (t, positions) in mods_by_type.iter().enumerate() {
            let shift = self.mod_types[t].shift;
            for &pos in positions {
                // add to any existing (candidate sites have no static mod, so this sets it)
                let cur = p.mods.get(&pos).copied().unwrap_or(0.0);
                p.set_mod(pos, cur + shift);
            }
        }
        p
    }

    /// Match a placement's fragment ladders to the spectrum with peak **claiming**.
    ///
    /// Generation order (load-bearing for the greedy claim): every configured N-terminal series
    /// (`--n-ions`, in order) first, then every C-terminal series (`--c-ions`); within a series, for
    /// each charge `1..=max_ch`, each cleavage site `i = 1..len`, each mod type. Each fragment claims
    /// the most-intense **unclaimed** peak in its ppm window; claimed peaks cannot be reused by later
    /// fragments — and the claim set is shared across all series and termini. The per-type inner loop
    /// lets a second peak in the same window be claimed for a second mod type (rare).
    fn evaluate_placement(
        &self,
        spectrum: &Spectrum,
        pep: &Peptide,
        mods_by_type: &[Vec<usize>],
        min_int: &mut f64,
        pref: &mut Vec<f64>,
        suf: &mut Vec<f64>,
    ) -> Placement {
        // Precompute this placement's cumulative residue+mod masses once, so every fragment m/z in the
        // ladders below is an O(1) lookup instead of an O(num_aa) re-sum from scratch.
        fill_frag_sums(pep, pref, suf);
        let max_ch = self.init_max_charge();

        let mut acc = MatchAcc::default();

        // N-terminal series first, then C-terminal — series order matters for the greedy claim.
        let mut nterm_run = 0i32;
        for &kind in &self.nions {
            let run = self.match_ion_ladder(spectrum, kind, max_ch, mods_by_type, &mut acc, min_int, pref, suf);
            nterm_run = nterm_run.max(run);
        }
        let mut cterm_run = 0i32;
        for &kind in &self.cions {
            let run = self.match_ion_ladder(spectrum, kind, max_ch, mods_by_type, &mut acc, min_int, pref, suf);
            cterm_run = cterm_run.max(run);
        }

        Placement {
            mods_by_type: mods_by_type.to_vec(),
            matched: acc.matched,
            total: acc.for_evidence,
            count: acc.for_count,
            direct_total: acc.direct_evidence,
            direct_count: acc.direct_count,
            nterm_run,
            cterm_run,
        }
    }

    /// Match one fragment-ion series (`kind`) across all charges and cleavage sites, claiming peaks
    /// into `acc`. Returns the longest run of consecutive singly-charged fragments matched (the
    /// terminal run-length feeding the N-/C-term scores). N-terminal series take the first `i`
    /// residues with `has_mod` = "a mod lies in the prefix `[0, i)`"; C-terminal series take the last
    /// `len - i` residues with `has_mod` = "a mod lies in the suffix `[i, len)`".
    #[allow(clippy::too_many_arguments)]
    fn match_ion_ladder(
        &self,
        spectrum: &Spectrum,
        kind: IonKind,
        max_ch: u32,
        mods_by_type: &[Vec<usize>],
        acc: &mut MatchAcc,
        min_int: &mut f64,
        pref: &[f64],
        suf: &[f64],
    ) -> i32 {
        let len = self.len();
        let nterm = kind.is_nterm();
        let mut best_run = 0i32;

        for ch in 1..=max_ch {
            let zf = ch as f64;
            let mut last_i = 0usize;
            let mut curr = 0i32;
            for i in 1..len {
                // a1 ions (the a-ion at the first cleavage) are uninformative; PTMProphet skips them.
                if kind == IonKind::A && i == 1 {
                    continue;
                }
                let num_aa = if nterm { i } else { len - i };
                // Fragment m/z in O(1) from the placement's cumulative sums. The add sequence matches
                // `Peptide::monoisotopic_mz_fragment` exactly for N-terminal ions (bit-identical) and
                // to within ~1 ULP for C-terminal ions (`suf` accumulates from the C-terminus inward).
                let mut sum = if nterm { pref[num_aa] } else { suf[len - num_aa] };
                if !nterm {
                    sum += masses::MASS_H2O;
                }
                sum += zf * masses::PROTON_FRAG;
                match kind {
                    IonKind::A => sum -= masses::MASS_CO,
                    IonKind::C => sum += masses::MASS_NH3,
                    IonKind::Z => sum -= masses::MASS_NH3 - masses::PROTON_FRAG,
                    IonKind::X => sum += masses::MASS_CO - 2.0 * masses::MASS_H,
                    _ => {}
                }
                let mz = sum / zf;
                if mz <= 0.0 {
                    continue;
                }
                for mods_t in mods_by_type {
                    // has_mod: a mod of this type lies on the side this fragment covers.
                    let has_mod = if nterm {
                        mods_t.iter().any(|&m| m < i)
                    } else {
                        mods_t.iter().any(|&m| m >= i)
                    };
                    if let Some((peak_mz, inten)) =
                        find_unclaimed_peak(spectrum, mz, self.frag_tol, &acc.claimed)
                    {
                        acc.claimed.insert(peak_mz.to_bits());
                        acc.matched.insert(peak_mz.to_bits(), inten);
                        acc.for_evidence += inten;
                        acc.for_count += 1;
                        // PTMProphet's minInt_ update; never fires here (see `process`).
                        if inten > 0.0 && inten < *min_int {
                            *min_int = inten;
                        }
                        if has_mod {
                            acc.direct_evidence += inten;
                            acc.direct_count += 1;
                        }
                        if ch == 1 {
                            // consecutive-fragment run at charge 1, updated per claim.
                            if i as i64 - last_i as i64 == 1 {
                                last_i = i;
                                curr += 1;
                            } else {
                                last_i = i;
                                curr = 1;
                            }
                            if curr > best_run {
                                best_run = curr;
                            }
                        }
                    }
                }
            }
        }
        best_run
    }

    /// Run the full per-PSM localization: enumerate placements, evaluate, score.
    /// Fills `self.types[*].scores` with O/M/term scores (probabilities come later via EM).
    pub fn process(&mut self, spectrum: &Spectrum) {
        let n_types = self.mod_types.len();

        // Per-type combinations, then the Cartesian product.
        let per_type: Vec<Vec<Vec<usize>>> = (0..n_types)
            .map(|t| combinatorics::combinations(self.types[t].sites.len(), self.types[t].nmods))
            .collect();
        let global = combinatorics::cartesian(&per_type);

        // min_int is the smallest POSITIVE peak intensity in the whole spectrum — NOT the min of
        // matched peaks. `match_ion_ladder` keeps PTMProphet's update that lowers it to a weaker
        // matched peak, but every matched peak is a peak of this spectrum and so is never weaker:
        // min_int stays fixed for the PSM. It is the O-score denominator and the `<1` guard threshold.
        let mut min_int = spectrum.min_intensity();
        let mut placements: Vec<Placement> = Vec::new();
        // Scratch cumulative-mass buffers, reused across this PSM's placements (refilled per placement
        // by `fill_frag_sums`) so the O(1)-per-fragment m/z path allocates nothing per placement.
        let mut pref_buf: Vec<f64> = Vec::new();
        let mut suf_buf: Vec<f64> = Vec::new();

        for combo in &global {
            // Map per-type chosen indices → actual residue positions.
            let mut mods_by_type: Vec<Vec<usize>> = vec![Vec::new(); n_types];
            for t in 0..n_types {
                for &idx in &combo[t] {
                    mods_by_type[t].push(self.types[t].sites[idx]);
                }
            }
            // --no-stack: reject if two types land on the same residue.
            if self.nostack {
                let mut all: Vec<usize> = mods_by_type.iter().flatten().copied().collect();
                let n = all.len();
                all.sort_unstable();
                all.dedup();
                if all.len() != n {
                    continue;
                }
            }
            // Build the concrete placement peptide once: it feeds both the mass gate and the fragment
            // cumulative sums below (previously it was rebuilt separately for each).
            let pep = self.placement_peptide(&mods_by_type);
            // Mass gate: rebuilt neutral mass within <1 ppm of calc_neut_mass. The 1 ppm is fixed, as
            // in PTMProphet, whose check is a literal `ppm < 1` (its PPMTOL option never reaches it).
            if self.calc_neut_mass > 0.0 {
                let ppm = (self.calc_neut_mass - pep.monoisotopic_neutral_mass()).abs() * 1e6
                    / self.calc_neut_mass;
                if ppm >= 1.0 {
                    continue;
                }
            }
            placements.push(self.evaluate_placement(
                spectrum,
                &pep,
                &mods_by_type,
                &mut min_int,
                &mut pref_buf,
                &mut suf_buf,
            ));
        }

        if !min_int.is_finite() || min_int <= 0.0 {
            min_int = spectrum.min_intensity();
        }

        // Accumulate per-site maxima + score, per type.
        for t in 0..n_types {
            self.score_type(t, &placements, min_int);
        }
    }

    /// Per-site evidence accumulation + O/M/term scores for one type.
    fn score_type(&mut self, t: usize, placements: &[Placement], min_int: f64) {
        let sites = self.types[t].sites.clone();
        let nmods = self.types[t].nmods;
        if sites.is_empty() {
            return;
        }

        // argmax placement index for mod-here / no-mod-here, per candidate residue.
        let mut mod_arg: BTreeMap<usize, usize> = BTreeMap::new();
        let mut nomod_arg: BTreeMap<usize, usize> = BTreeMap::new();
        let mut max_ev: BTreeMap<usize, f64> = BTreeMap::new();
        let mut nomod_max_ev: BTreeMap<usize, f64> = BTreeMap::new();
        let mut max_cnt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut nomod_max_cnt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut dmax_ev: BTreeMap<usize, f64> = BTreeMap::new();
        let mut nomod_dmax_ev: BTreeMap<usize, f64> = BTreeMap::new();
        let mut dmax_cnt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut nomod_dmax_cnt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut max_nt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut nomod_max_nt: BTreeMap<usize, i32> = BTreeMap::new();
        let mut max_ct: BTreeMap<usize, i32> = BTreeMap::new();
        let mut nomod_max_ct: BTreeMap<usize, i32> = BTreeMap::new();

        for (pi, p) in placements.iter().enumerate() {
            let modset = &p.mods_by_type[t];
            // mod-here updates (every placement; argmax recorded for the evidence max).
            for &pos in modset {
                if p.total > *max_ev.get(&pos).unwrap_or(&-1.0) {
                    max_ev.insert(pos, p.total);
                    mod_arg.insert(pos, pi);
                }
                if p.count > *max_cnt.get(&pos).unwrap_or(&0) {
                    max_cnt.insert(pos, p.count);
                }
                if p.direct_total > *dmax_ev.get(&pos).unwrap_or(&-1.0) {
                    dmax_ev.insert(pos, p.direct_total);
                }
                if p.direct_count > *dmax_cnt.get(&pos).unwrap_or(&0) {
                    dmax_cnt.insert(pos, p.direct_count);
                }
                if p.nterm_run > *max_nt.get(&pos).unwrap_or(&0) {
                    max_nt.insert(pos, p.nterm_run);
                }
                if p.cterm_run > *max_ct.get(&pos).unwrap_or(&0) {
                    max_ct.insert(pos, p.cterm_run);
                }
            }
            // no-mod-here updates (fully-modded placements only; every combo places exactly nmods).
            if modset.len() == nmods {
                for &pos in &sites {
                    if modset.contains(&pos) {
                        continue;
                    }
                    if p.total > *nomod_max_ev.get(&pos).unwrap_or(&-1.0) {
                        nomod_max_ev.insert(pos, p.total);
                        nomod_arg.insert(pos, pi);
                    }
                    if p.count > *nomod_max_cnt.get(&pos).unwrap_or(&0) {
                        nomod_max_cnt.insert(pos, p.count);
                    }
                    if p.direct_total > *nomod_dmax_ev.get(&pos).unwrap_or(&-1.0) {
                        nomod_dmax_ev.insert(pos, p.direct_total);
                    }
                    if p.direct_count > *nomod_dmax_cnt.get(&pos).unwrap_or(&0) {
                        nomod_dmax_cnt.insert(pos, p.direct_count);
                    }
                    if p.nterm_run > *nomod_max_nt.get(&pos).unwrap_or(&0) {
                        nomod_max_nt.insert(pos, p.nterm_run);
                    }
                    if p.cterm_run > *nomod_max_ct.get(&pos).unwrap_or(&0) {
                        nomod_max_ct.insert(pos, p.cterm_run);
                    }
                }
            }
        }

        let fully_localized = nmods == sites.len();
        let mut scores = BTreeMap::new();
        for &pos in &sites {
            let mut s = SiteScore { shift: 'c', ..Default::default() };
            if fully_localized {
                s.oscore = 1.0;
                s.mscore = 1.0;
                s.direct_oscore = 1.0;
                s.direct_mscore = 1.0;
                s.cterm_score = 1.0;
                s.nterm_score = 1.0;
                scores.insert(pos, s);
                continue;
            }

            // Overlap (common matched intensity/count) between the mod-here and no-mod-here argmax.
            let (ov, cc) = match (mod_arg.get(&pos), nomod_arg.get(&pos)) {
                (Some(&a), Some(&b)) => common_intensity(&placements[a], &placements[b]),
                _ => (0.0, 0),
            };

            let me = *max_ev.get(&pos).unwrap_or(&0.0);
            let ne = *nomod_max_ev.get(&pos).unwrap_or(&0.0);
            s.oscore = oscore(me - ov, ne - ov, min_int, self.min_o);

            let mc = *max_cnt.get(&pos).unwrap_or(&0) - cc;
            let nc = *nomod_max_cnt.get(&pos).unwrap_or(&0) - cc;
            s.mscore = ratio_with_guards(mc as f64, nc as f64);

            // Direct (degenerate ~0.5 in our config; computed for completeness).
            let dme = *dmax_ev.get(&pos).unwrap_or(&0.0);
            let dne = *nomod_dmax_ev.get(&pos).unwrap_or(&0.0);
            s.direct_oscore = oscore(dme, dne, min_int, 0);
            let dmc = *dmax_cnt.get(&pos).unwrap_or(&0);
            let dnc = *nomod_dmax_cnt.get(&pos).unwrap_or(&0);
            s.direct_mscore = ratio_with_guards(dmc as f64, dnc as f64);

            // Terminal run-length scores (no pseudo-count; 0 if both runs 0).
            let (cm, cu) = (*max_ct.get(&pos).unwrap_or(&0), *nomod_max_ct.get(&pos).unwrap_or(&0));
            s.cterm_score = if cm + cu > 0 { cm as f64 / (cm + cu) as f64 } else { 0.0 };
            let (nm, nu) = (*max_nt.get(&pos).unwrap_or(&0), *nomod_max_nt.get(&pos).unwrap_or(&0));
            s.nterm_score = if nm + nu > 0 { nm as f64 / (nm + nu) as f64 } else { 0.0 };

            s.shift = shift_direction(pos, mod_arg.get(&pos), nomod_arg.get(&pos), placements, t);
            scores.insert(pos, s);
        }
        self.types[t].scores = scores;
    }
}

/// Highest fragment charge matched: the precursor charge plus `max_frag_z` (0 or below; a positive
/// `--max-frag-z` is rejected at startup), at least 1 and at most the precursor charge (1 for a
/// precursor charge of 0). PTMProphetMpx.h `initMaxCharge` computes it in unsigned arithmetic and
/// differs in two cases: a `max_frag_z` below −charge wraps around and leaves the precursor charge
/// (PTMProphet-rs: 1), and a precursor charge of 0 gives 0, so no fragment ion is matched
/// (PTMProphet-rs: 1).
fn max_fragment_charge(charge: i32, max_frag_z: i32) -> u32 {
    let mut m = charge.saturating_add(max_frag_z).max(1);
    if m > charge {
        m = charge.max(1);
    }
    m as u32
}

/// O-score `Om/(Om+Ou)` with overlap-subtracted, min_int-normalized evidence and the
/// clamp → `+= min_o` → cumulative `<1` add-one guard sequence.
fn oscore(mod_ev: f64, nomod_ev: f64, min_int: f64, min_o: i32) -> f64 {
    let mut om = (mod_ev / min_int).max(0.0) + min_o as f64;
    let mut ou = (nomod_ev / min_int).max(0.0) + min_o as f64;
    if om < 1.0 {
        om += 1.0;
        ou += 1.0;
    }
    if ou < 1.0 {
        ou += 1.0;
        om += 1.0;
    }
    if om + ou <= 0.0 {
        return 0.0;
    }
    om / (om + ou)
}

/// Count-style ratio with the same cumulative `<1` add-one guards.
fn ratio_with_guards(mut cm: f64, mut cu: f64) -> f64 {
    if cm < 0.0 {
        cm = 0.0;
    }
    if cu < 0.0 {
        cu = 0.0;
    }
    if cm < 1.0 {
        cm += 1.0;
        cu += 1.0;
    }
    if cu < 1.0 {
        cu += 1.0;
        cm += 1.0;
    }
    if cm + cu <= 0.0 {
        return 0.0;
    }
    cm / (cm + cu)
}

/// Common matched intensity (and count) of peaks present in both placements.
/// Ascending-m/z iteration via the BTreeMap.
fn common_intensity(a: &Placement, b: &Placement) -> (f64, i32) {
    let mut sum = 0.0;
    let mut cnt = 0;
    for (k, &ai) in &a.matched {
        if let Some(&bi) = b.matched.get(k) {
            if ai > 0.0 && bi > 0.0 {
                sum += bi; // sum the second peptide's intensity
                cnt += 1;
            }
        }
    }
    (sum, cnt)
}

/// Mod shift direction: 'n' if more no-mod-here mods are N-terminal to `pos` than mod-here mods,
/// else 'c'. (Only the N-term comparison is observable.)
fn shift_direction(
    pos: usize,
    mod_arg: Option<&usize>,
    nomod_arg: Option<&usize>,
    placements: &[Placement],
    t: usize,
) -> char {
    let here_n = mod_arg
        .map(|&a| placements[a].mods_by_type[t].iter().filter(|&&m| m < pos).count())
        .unwrap_or(0);
    let not_n = nomod_arg
        .map(|&b| placements[b].mods_by_type[t].iter().filter(|&&m| m < pos).count())
        .unwrap_or(0);
    if not_n > here_n {
        'n'
    } else {
        'c'
    }
}

/// Fill `pref`/`suf` with one placement's cumulative residue+mod masses, enabling O(1) fragment m/z.
///
/// `pref[k]` = `nterm_mod + Σ_{j<k}(residue_j + mod_j)`, accumulated in ascending residue order — the
/// exact same f64 add sequence as the old N-terminal per-fragment sum, so b/a/c m/z stay
/// **bit-identical**. `suf[i]` = `cterm_mod + Σ_{j≥i}(residue_j + mod_j)`, accumulated from the
/// C-terminus inward so all suffix sums build in one pass; this reorders the adds versus the old
/// ascending C-terminal sum, so y/z/x m/z may differ by ~1 ULP — orders of magnitude below the ppm
/// match window (the full localized output is byte-identical on the test set). Both buffers are reused
/// across a PSM's placements; they are cleared and resized to `len + 1` here.
fn fill_frag_sums(pep: &Peptide, pref: &mut Vec<f64>, suf: &mut Vec<f64>) {
    let res = pep.residues();
    let len = res.len();
    pref.clear();
    pref.resize(len + 1, 0.0);
    suf.clear();
    suf.resize(len + 1, 0.0);

    // N-terminal prefix: start at nterm_mod, add residue then mod ascending (same order as before).
    pref[0] = pep.nterm_mod.unwrap_or(0.0);
    for k in 1..=len {
        let mut v = pref[k - 1] + masses::residue_mono_mass(res[k - 1]).unwrap_or(0.0);
        if let Some(&d) = pep.mods.get(&(k - 1)) {
            v += d;
        }
        pref[k] = v;
    }

    // C-terminal suffix: start at cterm_mod, accumulate residues from the C-terminus inward.
    suf[len] = pep.cterm_mod.unwrap_or(0.0);
    for i in (0..len).rev() {
        let mut v = suf[i + 1] + masses::residue_mono_mass(res[i]).unwrap_or(0.0);
        if let Some(&d) = pep.mods.get(&i) {
            v += d;
        }
        suf[i] = v;
    }
}

/// Among peaks in the OPEN window `(theor_mz-tol, theor_mz+tol)`, return the most-intense
/// **unclaimed** peak. `max` is initialized to the lowest-m/z peak in the window (claimed or not)
/// and only advanced to strictly-more-intense *unclaimed* peaks; if the resulting `max` is itself
/// claimed (or non-positive), return `None` (so a claimed lowest-peak that dominates all unclaimed
/// peaks blocks the match). Returns `(peak_mz, intensity)`.
fn find_unclaimed_peak(
    spectrum: &Spectrum,
    theor_mz: f64,
    ppm: f64,
    claimed: &std::collections::HashSet<u64>,
) -> Option<(f64, f64)> {
    let tol = theor_mz * ppm / 1e6;
    let lo = theor_mz - tol;
    let hi = theor_mz + tol;
    let mz = &spectrum.mz;
    // low = first peak strictly > lo (upper_bound); window is [low, first peak >= hi).
    let start = mz.partition_point(|&m| m <= lo);
    if start >= mz.len() || mz[start] >= hi {
        return None;
    }
    let mut max_i = start; // initialized to lowest-m/z peak (claimed or not)
    let mut i = start;
    while i < mz.len() && mz[i] < hi {
        if !claimed.contains(&mz[i].to_bits()) && spectrum.intensity[i] > spectrum.intensity[max_i] {
            max_i = i;
        }
        i += 1;
    }
    if claimed.contains(&mz[max_i].to_bits()) || spectrum.intensity[max_i] <= 0.0 {
        return None;
    }
    Some((mz[max_i], spectrum.intensity[max_i]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectrum::mzml::MzmlReader;
    use std::path::{Path, PathBuf};

    fn cfg() -> Config {
        Config {
            mods: Vec::new(),
            mod_string: String::new(),
            opt_string: String::new(),
            frag_tol: 10.0,
            em: 1,
            min_prob: 0.5,
            mod_prec: 0,
            min_o: 0,
            max_frag_z: 0,
            threads: None,
            nostack: true,
            static_tol: true,
            nions: vec![IonKind::B],
            cions: vec![IonKind::Y],
            inputs: Vec::new(),
            in_file: PathBuf::new(),
            out_file: PathBuf::new(),
        }
    }

    fn fragpipe_mods() -> Vec<ModType> {
        vec![
            ModType { aas: vec!['M'], shift: 15.9949, ptm_label: "M:15.9949".into() },
            ModType { aas: vec!['S', 'T', 'Y'], shift: 79.96633, ptm_label: "STY:79.96633".into() },
        ]
    }

    fn tmt_tables() -> ModTables {
        let mut t = ModTables::default();
        t.stat_mods.insert('C', 57.02146);
        t.stat_mods.insert('K', 229.16293);
        t.stat_mods.insert('n', 229.16293);
        t.var_mods.insert('M', vec![15.9949]);
        t.var_mods.insert('S', vec![79.96633]);
        t.var_mods.insert('T', vec![79.96633]);
        t.var_mods.insert('Y', vec![79.96633]);
        t
    }

    /// `--max-frag-z` 0 or below: the precursor charge plus Z, from 1 up to the precursor charge.
    /// Where PTMProphet's unsigned `initMaxCharge` differs (Z below −charge; precursor charge 0),
    /// PTMProphet-rs gives 1.
    #[test]
    fn max_fragment_charge_adds_z_to_the_precursor_charge() {
        assert_eq!(max_fragment_charge(3, 0), 3);
        assert_eq!(max_fragment_charge(3, -1), 2);
        assert_eq!(max_fragment_charge(3, -2), 1);
        assert_eq!(max_fragment_charge(3, -3), 1);
        assert_eq!(max_fragment_charge(3, -4), 1); // PTMProphet: 3
        assert_eq!(max_fragment_charge(3, i32::MIN), 1); // PTMProphet: 3
        assert_eq!(max_fragment_charge(1, 0), 1);
        assert_eq!(max_fragment_charge(0, 0), 1); // PTMProphet: 0
        assert_eq!(max_fragment_charge(0, -1), 1); // PTMProphet: 0
    }

    #[test]
    fn sites_and_nmods_scan3185() {
        // No spectrum needed: just the site identification.
        let mut pep_mods = BTreeMap::new();
        pep_mods.insert(0usize, 230.17076); // n-term TMT
        pep_mods.insert(5usize, 166.99835); // S5 phospho (1-based)
        let mt = fragpipe_mods();
        let tb = tmt_tables();
        let m = Mpx::build(&cfg(), &mt, &tb, "ALRASESGI", &pep_mods, 2, 0.0);
        // STY is type index 1: candidate S at residue indices 4 and 6; one mod present.
        assert_eq!(m.types[1].sites, vec![4, 6]);
        assert_eq!(m.types[1].nmods, 1);
        // M (type 0): no candidate residues.
        assert!(m.types[0].sites.is_empty());
        assert_eq!(m.types[0].nmods, 0);
        // n-term TMT classified static.
        assert!((m.base.nterm_mod.unwrap() - 229.16293).abs() < 1e-4);
    }

    #[test]
    fn oscore_scan3185_vs_oracle() {
        let mzml = Path::new("F:/data/ccrcc_dda_phospho/mzml/01CPTAC_CCRCC_Phosphoproteome_JHU_20171106/01CPTAC_CCRCC_P_JHU_20171106_LUMOS_f05_calibrated.mzML");
        if !mzml.exists() {
            eprintln!("[skip] mzML not present");
            return;
        }
        let reader = MzmlReader::open(mzml).unwrap();
        let spec = reader.spectrum_for_scan(3185).unwrap();

        let mut pep_mods = BTreeMap::new();
        pep_mods.insert(0usize, 230.17076);
        pep_mods.insert(5usize, 166.99835);
        let mt = fragpipe_mods();
        let tb = tmt_tables();
        // calc_neut_mass=0 disables the mass gate; both STY placements are mass-equal & valid.
        let mut m = Mpx::build(&cfg(), &mt, &tb, "ALRASESGI", &pep_mods, 2, 0.0);
        m.process(&spec);

        let sty = &m.types[1];
        let o5 = sty.scores[&4].oscore; // residue index 4 → output position 5
        let o7 = sty.scores[&6].oscore; // residue index 6 → output position 7
        eprintln!("[scan3185] oscore pos5={o5:.4} pos7={o7:.4} (oracle 0.2724 / 0.7276)");
        // First-pass tolerance; tighten once the peak-matching/overlap nuances are pinned.
        assert!((o5 - 0.2724).abs() < 0.05, "pos5 oscore {o5}");
        assert!((o7 - 0.7276).abs() < 0.05, "pos7 oscore {o7}");
        assert!((o5 + o7 - 1.0).abs() < 1e-6, "should sum to 1");
    }
}
