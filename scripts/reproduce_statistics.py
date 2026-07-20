#!/usr/bin/env python3
"""
Recompute every statistic reported in the manuscript from the aggregate
per-target metrics in ../data/, and print it next to the value stated in
the text so that any discrepancy is immediately visible.

Usage:
    python scripts/reproduce_statistics.py

Requires: pandas, numpy, scipy
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path(__file__).resolve().parent.parent / "data"

G_VS = pd.read_csv(DATA / "gpcr_VS_FINAL_20.csv")
G_TAIL = pd.read_csv(DATA / "gpcr_tail_FINAL_20.csv")
S_VS = pd.read_csv(DATA / "slc_VS_FINAL_23.csv")
CLS = pd.read_csv(DATA / "clsprob_FINAL.csv")
G_SEP = pd.read_csv(DATA / "sep_decomp_gpcr.csv")
S_SEP = pd.read_csv(DATA / "sep_decomp_slc.csv")
MECH = pd.read_csv(DATA / "mechanism_slc.csv")

_checks = []


def report(label, computed, stated, fmt="{:.3f}", tol=0.01):
    """Print a computed value beside the value stated in the manuscript."""
    ok = stated is None or abs(computed - stated) <= tol
    _checks.append(ok)
    c = fmt.format(computed)
    s = "—" if stated is None else fmt.format(stated)
    flag = "" if ok else "   <-- DIFFERS"
    print(f"  {label:<52} computed {c:>10}   manuscript {s:>10}{flag}")


def cohens_d(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    s = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                / (len(a) + len(b) - 2))
    return (a.mean() - b.mean()) / s


print(__doc__.split("Usage")[0].strip())
print()
print("Datasets loaded:")
print(f"  GPCR virtual screening      n = {len(G_VS)}")
print(f"  SLC  virtual screening      n = {len(S_VS)}")
print(f"  Two-output comparison       n = {len(CLS)}")
print()

# ---------------------------------------------------------------- Table 1 / Fig 3
print("Virtual-screening outcome (Fig. 3, Table 1)")
report("GPCR successes (AUROC >= 0.7)", int((G_VS.AUROC >= 0.7).sum()), 7, "{:.0f}")
report("SLC successes (AUROC >= 0.7)", int((S_VS.AUROC >= 0.7).sum()), 2, "{:.0f}")
report("GPCR mean AUROC", G_VS.AUROC.mean(), 0.64, "{:.2f}")
report("SLC mean AUROC", S_VS.AUROC.mean(), 0.56, "{:.2f}")
report("GPCR targets with EF@1% = 0", int((G_VS["EF@1%"] == 0).sum()), 8, "{:.0f}")
report("SLC targets with EF@1% = 0", int((S_VS["EF@1%"] == 0).sum()), 15, "{:.0f}")
print()

# ---------------------------------------------------------------- Fig 4
print("Target-side descriptors versus screening performance (Fig. 4)")
print("  (the 'r' column is per-target ranking accuracy from the binder-only")
print("   prediction run; see Methods and slc_ranking_accuracy_82.csv)")
r, p = stats.pearsonr(G_VS.r.dropna(), G_VS.loc[G_VS.r.notna(), "AUROC"])
report("ranking accuracy -> AUROC (GPCR), r", r, 0.60, "{:.2f}")
report("ranking accuracy -> AUROC (GPCR), P", p, 0.006, "{:.3f}")

sub = G_VS[["r", "EF@1%"]].dropna()
r, p = stats.pearsonr(sub.r, sub["EF@1%"])
report("ranking accuracy -> EF@1% (GPCR), P", p, 0.055, "{:.3f}")

sub = S_VS[["r", "AUROC"]].dropna()
r, p = stats.pearsonr(sub.r, sub.AUROC)
report("ranking accuracy -> AUROC (SLC), n", len(sub), 23, "{:.0f}")
report("ranking accuracy -> AUROC (SLC), r", r, 0.17, "{:.2f}")
report("ranking accuracy -> AUROC (SLC), P", p, 0.44, "{:.2f}")

sub = S_VS[["r", "EF@1%"]].dropna()
r, p = stats.pearsonr(sub.r, sub["EF@1%"])
report("ranking accuracy -> EF@1% (SLC), r", r, 0.27, "{:.2f}")
report("ranking accuracy -> EF@1% (SLC), P", p, 0.21, "{:.2f}")

for col, label, stated_p in [
    ("fpocket_drug_score", "fpocket score -> AUROC (SLC), P", None),
    ("p2rank_top_score", "P2Rank score -> AUROC (SLC), P", None),
    ("plddt_mean", "pLDDT -> AUROC (SLC), P", None),
]:
    if col in S_VS.columns:
        sub = S_VS[[col, "AUROC"]].dropna()
        if len(sub) >= 5:
            r, p = stats.pearsonr(sub[col], sub.AUROC)
            report(label, p, stated_p, "{:.2f}")
print("  (structural descriptors are distributed with Supplementary Table S1)")
print()

# ---------------------------------------------------------------- Fig 5
print("Score separation (Fig. 5)")
pool = pd.concat([
    G_TAIL[["Target", "Sep"]].merge(G_VS[["Target", "AUROC"]], on="Target"),
    S_VS[["Target", "Sep", "AUROC"]],
], ignore_index=True)
r, p = stats.pearsonr(pool.Sep, pool.AUROC)
report("Sep -> AUROC, pooled (n = 43), r", r, 0.97, "{:.2f}")
report("Sep -> AUROC, pooled, R^2", r ** 2, 0.94, "{:.2f}")
report("GPCR Sep maximum", G_TAIL.Sep.max(), 1.47, "{:.2f}")
report("SLC Sep maximum", S_VS.Sep.max(), 0.52, "{:.2f}")
report("GPCR targets with Sep > 0.55", int((G_TAIL.Sep > 0.55).sum()), 5, "{:.0f}")
report("SLC targets with Sep > 0.55", int((S_VS.Sep > 0.55).sum()), 0, "{:.0f}")
report("Sep variance, Levene P", stats.levene(G_TAIL.Sep, S_VS.Sep).pvalue, 0.018, "{:.3f}")
print()

# ---------------------------------------------------------------- Fig 6
print("Predicted scores of experimental binders (Fig. 6)")
report("GPCR targets with binder_mean > 0", int((G_SEP.binder_mean > 0).sum()), 13, "{:.0f}")
report("SLC targets with binder_mean > 0", int((S_SEP.binder_mean > 0).sum()), 3, "{:.0f}")
report("GPCR mean binder score", G_SEP.binder_mean.mean(), 0.30, "{:.2f}")
report("SLC mean binder score", S_SEP.binder_mean.mean(), -0.57, "{:.2f}")
report("binder score, Cohen's d", cohens_d(G_SEP.binder_mean, S_SEP.binder_mean), 1.48, "{:.2f}")

gap_g = G_SEP.binder_mean - G_SEP.nb_mean
gap_s = S_SEP.binder_mean - S_SEP.nb_mean
report("GPCR mean binder-non-binder gap", gap_g.mean(), 0.55, "{:.2f}")
report("SLC mean binder-non-binder gap", gap_s.mean(), 0.14, "{:.2f}")
report("gap difference, Mann-Whitney P",
       stats.mannwhitneyu(gap_g, gap_s).pvalue, 0.07, "{:.2f}")

if "iptm" in MECH.columns:
    sub = MECH[["iptm", "Sep"]].dropna()
    r, p = stats.pearsonr(sub.iptm, sub.Sep)
    report("ipTM -> Sep (SLC), r", r, 0.26, "{:.2f}")
    report("ipTM -> Sep (SLC), P", p, 0.24, "{:.2f}")
    report("SLC mean ipTM", sub.iptm.mean(), 0.888, "{:.3f}")
print()

# ---------------------------------------------------------------- Fig 7
print("Robustness to the choice of scoring output (Fig. 7, Table 1)")
g = CLS[CLS.Class == "GPCR"]
s = CLS[CLS.Class == "SLC"]
report("GPCR successes, affinity value", int((g.AUROC_val >= 0.7).sum()), 7, "{:.0f}")
report("GPCR successes, classification prob", int((g.AUROC_prob >= 0.7).sum()), 10, "{:.0f}")
report("SLC successes, affinity value", int((s.AUROC_val >= 0.7).sum()), 2, "{:.0f}")
report("SLC successes, classification prob", int((s.AUROC_prob >= 0.7).sum()), 2, "{:.0f}")
report("GPCR two outputs, Wilcoxon P",
       stats.wilcoxon(g.AUROC_val, g.AUROC_prob).pvalue, 0.002, "{:.3f}")
report("SLC two outputs, Wilcoxon P",
       stats.wilcoxon(s.AUROC_val, s.AUROC_prob).pvalue, 0.45, "{:.2f}")
report("GPCR vs SLC (affinity value), Mann-Whitney P",
       stats.mannwhitneyu(g.AUROC_val, s.AUROC_val).pvalue, 0.17, "{:.2f}")
report("GPCR vs SLC (classification prob), Mann-Whitney P",
       stats.mannwhitneyu(g.AUROC_prob, s.AUROC_prob).pvalue, 0.03, "{:.2f}")
print()

n_ok, n_all = sum(_checks), len(_checks)
print(f"{n_ok} of {n_all} values match the manuscript within tolerance.")
if n_ok < n_all:
    print("Entries marked DIFFERS should be checked against the manuscript text.")
