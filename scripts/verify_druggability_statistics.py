#!/usr/bin/env python3
"""
Verify the experimental-druggability statistics reported in the manuscript.

Recomputes, from data/druggability_by_fold_592.csv alone:

  - hit recovery for GPCRs and for each SLC structural fold class
  - Fisher exact test, GPCRs vs accessible-cavity SLC transporters
  - chi-square test across the four SLC fold classes

and prints each result next to the value stated in the manuscript, so that any
discrepancy is visible without reading the code.

Expected input schema (one row per target):

    uniprot        UniProt accession                        e.g. Q9NQQ7
    gene_symbol    HGNC-approved gene symbol                e.g. SLC35H1
    family         "GPCR" or "SLC"
    fold_class     for SLC: one of
                     accessible_cavity
                     compact_gated
                     elevator_like
                     ion_channel_like
                   for GPCR: leave empty
    hit            1 if the screen returned at least one binder, else 0

Usage:
    python scripts/verify_druggability_statistics.py
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

DATA = Path(__file__).resolve().parent.parent / "data" / "druggability_by_fold_592.csv"

FOLD_ORDER = [
    ("accessible_cavity", "Accessible cavity"),
    ("compact_gated", "Compact gated"),
    ("elevator_like", "Elevator-like"),
    ("ion_channel_like", "Ion channel-like"),
]

# values as stated in the manuscript
REPORTED = {
    "gpcr": (213, 296, 72.0),
    "accessible_cavity": (129, 181, 71.3),
    "compact_gated": (30, 66, 45.5),
    "elevator_like": (15, 36, 41.7),
    "ion_channel_like": (2, 13, 15.4),
    "fisher_or": 1.03,
    "fisher_p": 0.92,
    "chi2": 31.0,
    "chi2_dof": 3,
    "chi2_p": 8.3e-7,
}

failures = []


def check(label, computed, reported, fmt="{:.1f}", tol=0.05):
    ok = abs(computed - reported) <= tol
    mark = "OK " if ok else "!! "
    print(
        f"  {mark}{label:<44} computed {fmt.format(computed):>10}"
        f"   manuscript {fmt.format(reported):>10}"
    )
    if not ok:
        failures.append(label)


def main():
    if not DATA.exists():
        sys.exit(
            f"Input file not found: {DATA}\n"
            "See the docstring of this script for the expected column schema."
        )

    df = pd.read_csv(DATA)

    required = {"family", "fold_class", "hit"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"Missing required column(s): {', '.join(sorted(missing))}")

    df["fold_class"] = df["fold_class"].fillna("").astype(str).str.strip()
    df["hit"] = df["hit"].astype(int)

    gpcr = df[df["family"] == "GPCR"]
    slc = df[df["family"] == "SLC"]

    print("\nTarget counts")
    print(f"  GPCR targets                  {len(gpcr)}")
    print(f"  SLC targets                   {len(slc)}")
    print(f"  total                         {len(df)}")

    classified = slc[slc["fold_class"] != ""]
    if len(classified) != len(slc):
        print(
            f"\n  note: {len(slc) - len(classified)} SLC target(s) have no fold class "
            "and are excluded from the fold analysis"
        )

    # ---- hit recovery -------------------------------------------------
    print("\nHit recovery (%)")
    gp_h, gp_n = int(gpcr["hit"].sum()), len(gpcr)
    check(f"GPCR  ({gp_h}/{gp_n})", 100 * gp_h / gp_n, REPORTED["gpcr"][2])

    counts = {}
    for key, label in FOLD_ORDER:
        sub = classified[classified["fold_class"] == key]
        h, n = int(sub["hit"].sum()), len(sub)
        counts[key] = (h, n)
        if n == 0:
            print(f"  !! {label:<44} no rows found")
            failures.append(label)
            continue
        check(f"{label}  ({h}/{n})", 100 * h / n, REPORTED[key][2])

    # ---- Fisher exact: GPCR vs accessible-cavity SLC ------------------
    ac_h, ac_n = counts.get("accessible_cavity", (0, 0))
    print("\nFisher exact test: GPCR vs accessible-cavity SLC transporters")
    if ac_n:
        table = [[gp_h, gp_n - gp_h], [ac_h, ac_n - ac_h]]
        odds, p = fisher_exact(table)
        check("odds ratio", odds, REPORTED["fisher_or"], "{:.2f}", tol=0.005)
        check("P", p, REPORTED["fisher_p"], "{:.2f}", tol=0.005)

    # ---- chi-square across fold classes -------------------------------
    print("\nChi-square test across SLC fold classes")
    table = [[h, n - h] for h, n in (counts[k] for k, _ in FOLD_ORDER) if n]
    if len(table) == len(FOLD_ORDER):
        chi2, p, dof, _ = chi2_contingency(table, correction=False)
        check("chi-square", chi2, REPORTED["chi2"], "{:.1f}", tol=0.05)
        check("degrees of freedom", dof, REPORTED["chi2_dof"], "{:.0f}", tol=0)
        check("P", p, REPORTED["chi2_p"], "{:.1e}", tol=1e-7)

    # ---- summary ------------------------------------------------------
    print()
    if failures:
        print(f"{len(failures)} value(s) did not match the manuscript:")
        for f in failures:
            print(f"  - {f}")
        sys.exit(1)
    print("All druggability statistics match the values reported in the manuscript.")


if __name__ == "__main__":
    main()
