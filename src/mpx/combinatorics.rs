// SPDX-License-Identifier: LGPL-3.0-or-later
// SPDX-FileCopyrightText: 2026 The Regents of the University of Michigan
// SPDX-FileCopyrightText: 2011 David Shteynberg
// Translated into Rust from TPP 6.3.2 src/Validation/PTMProphetParser/PTMProphetMpx.cpp
// (LGPL-2.1-or-later); modified by the University of Michigan, 2026.

//! Mod-placement combination enumeration.
//!
//! `combinations(n, k)` returns every k-subset of `{0..n}` in **lexicographic order**,
//! e.g. C(4,2) = `[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]`. The full
//! per-PSM placement set is the Cartesian product of each type's `combinations(|sites|, nmods)`
//! (an odometer with the last type varying fastest) — built in `mpx::Mpx`.

/// All `k`-subsets of `{0,1,…,n-1}` in lexicographic order. `k==0` → one empty subset; `k>n` → none.
pub fn combinations(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![Vec::new()];
    }
    if k > n {
        return Vec::new();
    }
    let mut result = Vec::new();
    let mut combo: Vec<usize> = (0..k).collect();
    loop {
        result.push(combo.clone());
        // advance to the lexicographic successor
        let mut t = k - 1;
        while combo[t] == n - k + t {
            if t == 0 {
                return result;
            }
            t -= 1;
        }
        combo[t] += 1;
        for j in (t + 1)..k {
            combo[j] = combo[j - 1] + 1;
        }
    }
}

/// Cartesian product of per-type index lists (odometer; last list varies fastest).
/// `lists[t]` are the candidate combination-indices for type `t`; returns each global tuple.
pub fn cartesian(lists: &[Vec<Vec<usize>>]) -> Vec<Vec<Vec<usize>>> {
    let mut out: Vec<Vec<Vec<usize>>> = vec![Vec::new()];
    for list in lists {
        let mut next = Vec::with_capacity(out.len() * list.len().max(1));
        for prefix in &out {
            if list.is_empty() {
                // a type with no combinations contributes an empty placement
                let mut p = prefix.clone();
                p.push(Vec::new());
                next.push(p);
            } else {
                for combo in list {
                    let mut p = prefix.clone();
                    p.push(combo.clone());
                    next.push(p);
                }
            }
        }
        out = next;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c_4_2_lexicographic() {
        assert_eq!(
            combinations(4, 2),
            vec![vec![0, 1], vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3], vec![2, 3]]
        );
    }

    #[test]
    fn c_2_1() {
        assert_eq!(combinations(2, 1), vec![vec![0], vec![1]]);
    }

    #[test]
    fn edge_cases() {
        assert_eq!(combinations(3, 0), vec![Vec::<usize>::new()]);
        assert!(combinations(2, 3).is_empty());
        assert_eq!(combinations(3, 3), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn cartesian_two_types() {
        // type0: C(2,1)=[[0],[1]]; type1: C(1,1)=[[0]]  → 2 placements, last type fastest.
        let lists = vec![combinations(2, 1), combinations(1, 1)];
        let prod = cartesian(&lists);
        let expected: Vec<Vec<Vec<usize>>> = vec![vec![vec![0], vec![0]], vec![vec![1], vec![0]]];
        assert_eq!(prod, expected);
    }
}
