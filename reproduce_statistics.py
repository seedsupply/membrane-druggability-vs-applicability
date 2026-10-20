#!/usr/bin/env python3
"""
Recompute the statistics reported for the virtual-screening evaluation from
the aggregate per-target metrics in ../data/, and print it next to the value
stated in the manuscript so that any discrepancy is immediately visible.

Usage:
    python scripts/reproduce_statistics.py

Requires: pandas, numpy, scipy
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.isotonic import IsotonicRegression

DATA = Path(__file__).resolve().parent.parent / "data"

MET = pd.read_csv(DATA / "per_target_metrics.csv")
SD = pd.read_csv(DATA / "score_distributions.csv")
CLS = pd.read_csv(DATA / "clsprob_FINAL.csv")
DESC = pd.read_csv(DATA / "target_descriptors.csv")
TRAIN = pd.read_csv(DATA / "training_representation.csv")
CC = pd.read_csv(DATA / "contact_consistency_24.csv")
HELD_T = pd.read_csv(DATA / "sep_heldout_target.csv")

G = MET[MET.family == "GPCR"]
S = MET[MET.family == "SLC"]

_checks = []


def report(label, computed, stated, fmt="{:.3f}", tol=0.01):
    ok = stated is None or abs(computed - stated) <= tol
    _checks.append(ok)
    c = fmt.format(computed)
    s = "—" if stated is None else fmt.format(stated)
    flag = "" if ok else "   <-- DIFFERS"
    print(f"  {label:<54} computed {c:>10}   manuscript {s:>10}{flag}")


def bh(pvals):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for k in range(n - 1, -1, -1):
        j = order[k]
        prev = min(prev, p[j] * n / (k + 1))
        q[j] = prev
    return q


print(__doc__.split("Usage")[0].strip())
print()
print("Datasets loaded:")
print(f"  per-target metrics          n = {len(MET)}  "
      f"(GPCR {len(G)}, SLC {len(S)})")
print(f"  score distributions         n = {len(SD)}")
print(f"  two-output comparison       n = {len(CLS)}")
print(f"  target-side descriptors     n = {len(DESC)}")
print(f"  training-domain proxies     n = {len(TRAIN)}")
print(f"  interface consistency       n = {len(CC)}")
print()

# ================================================================ Results 3
print("Early enrichment under the evaluated predictor (Results 3, Fig. 3, Table 1)")
report("GPCR mean EF@0.5%", G.EF05.mean(), 20.2, "{:.1f}", tol=0.05)
report("SLC  mean EF@0.5%", S.EF05.mean(), 4.3, "{:.1f}", tol=0.05)
report("GPCR median EF@0.5%", G.EF05.median(), 13.1, "{:.1f}", tol=0.05)
report("SLC  median EF@0.5%", S.EF05.median(), 0.0, "{:.1f}", tol=0.05)
report("EF@0.5%, Mann-Whitney P",
       stats.mannwhitneyu(G.EF05, S.EF05, alternative="two-sided").pvalue,
       0.011, "{:.3f}")
report("GPCR mean EF@1%", G.EF1.mean(), 13.9, "{:.1f}", tol=0.05)
report("SLC  mean EF@1%", S.EF1.mean(), 3.6, "{:.1f}", tol=0.05)
report("EF@1%, Mann-Whitney P",
       stats.mannwhitneyu(G.EF1, S.EF1, alternative="two-sided").pvalue,
       0.035, "{:.3f}")
report("GPCR mean BEDROC (a = 80.5)", G["BEDROC80.5"].mean(), 0.204, "{:.3f}")
report("SLC  mean BEDROC (a = 80.5)", S["BEDROC80.5"].mean(), 0.056, "{:.3f}")
report("BEDROC (a = 80.5), Mann-Whitney P",
       stats.mannwhitneyu(G["BEDROC80.5"], S["BEDROC80.5"],
                          alternative="two-sided").pvalue, 0.082, "{:.3f}")
report("GPCR mean AUROC", G.AUROC.mean(), 0.64, "{:.2f}")
report("SLC  mean AUROC", S.AUROC.mean(), 0.56, "{:.2f}")
report("AUROC, Mann-Whitney P",
       stats.mannwhitneyu(G.AUROC, S.AUROC, alternative="two-sided").pvalue,
       0.17, "{:.2f}")
print()

print("  BEDROC sensitivity (Suppl. Fig. S1)")
report("GPCR mean BEDROC (a = 20)", G.BEDROC20.mean(), 0.244, "{:.3f}")
report("SLC  mean BEDROC (a = 20)", S.BEDROC20.mean(), 0.097, "{:.3f}")
report("BEDROC (a = 20), Mann-Whitney P",
       stats.mannwhitneyu(G.BEDROC20, S.BEDROC20,
                          alternative="two-sided").pvalue, 0.111, "{:.3f}")
report("GPCR mean BEDROC (a = 321.9)", G["BEDROC321.9"].mean(), 0.224, "{:.3f}")
report("SLC  mean BEDROC (a = 321.9)", S["BEDROC321.9"].mean(), 0.043, "{:.3f}")
report("BEDROC (a = 321.9), Mann-Whitney P",
       stats.mannwhitneyu(G["BEDROC321.9"], S["BEDROC321.9"],
                          alternative="two-sided").pvalue, 0.056, "{:.3f}")
print()

print("  Per-target bootstrap intervals (Fig. 3b)")
report("GPCR targets, EF@0.5% lower bound > 1",
       int((G.EF05_lo > 1).sum()), 4, "{:.0f}")
report("SLC  targets, EF@0.5% lower bound > 1",
       int((S.EF05_lo > 1).sum()), 0, "{:.0f}")
report("GPCR targets, EF@1% lower bound > 1",
       int((G.EF1_lo > 1).sum()), 4, "{:.0f}")
report("SLC  targets, EF@1% lower bound > 1",
       int((S.EF1_lo > 1).sum()), 0, "{:.0f}")
report("GPCR targets, AUROC lower bound > 0.5",
       int((G.AUROC_lo > 0.5).sum()), 9, "{:.0f}")
report("SLC  targets, AUROC lower bound > 0.5",
       int((S.AUROC_lo > 0.5).sum()), 6, "{:.0f}")
print()

print("  Binder prevalence and attainable enrichment")
report("minimum prevalence (%)", MET.prevalence.min() * 100, 0.7, "{:.1f}", tol=0.05)
report("maximum prevalence (%)", MET.prevalence.max() * 100, 20.6, "{:.1f}", tol=0.05)
print()

print("  Threshold sensitivity (Fig. 3d, Suppl. Table S5)")
print("  The AUROC >= 0.7 criterion is retained for comparability only.")
tbl = []
for cut in np.round(np.arange(0.55, 0.8501, 0.01), 2):
    a = int((G.AUROC >= cut).sum())
    b = int((S.AUROC >= cut).sum())
    p = stats.fisher_exact([[a, len(G) - a], [b, len(S) - b]])[1]
    tbl.append((cut, a, b, p))
t = {c: (a, b, p) for c, a, b, p in tbl}
report("targets above 0.70 (GPCR)", t[0.70][0], 7, "{:.0f}")
report("targets above 0.70 (SLC)", t[0.70][1], 2, "{:.0f}")
report("Fisher P at cutoff 0.70", t[0.70][2], 0.059, "{:.3f}")
report("Fisher P at cutoff 0.68", t[0.68][2], 0.039, "{:.3f}")
report("Fisher P at cutoff 0.69", t[0.69][2], 0.148, "{:.3f}")
report("Fisher P at cutoff 0.78", t[0.78][2], 0.016, "{:.3f}")
cross = sum(1 for i in range(1, len(tbl))
            if (tbl[i - 1][3] < 0.05) != (tbl[i][3] < 0.05))
report("times the P = 0.05 level is crossed", cross, 3, "{:.0f}")
print()

print("  Targets below the cutoff that still enrich (Results 3)")
low = S[(S.AUROC < 0.7) & (S.EF1 > 0)]
report("SLC targets, AUROC < 0.7 but EF@1% > 0", len(low), 6, "{:.0f}")
for name, val in [("SLC6A12", 9.8), ("SLC7A4", 8.9), ("SLC7A14", 8.2)]:
    row = S[S.target == name]
    if len(row):
        report(f"  {name} EF@1%", float(row.EF1.iloc[0]), val, "{:.1f}", tol=0.05)
print()

# ================================================================ Results 4
print("Target-side descriptors (Results 4, Fig. 4, Suppl. Table S6)")
desc_cols = [("fpocket_druggability", "fpocket druggability score"),
             ("fpocket_pockets_n", "fpocket pocket count"),
             ("fpocket_top_volume", "top pocket volume"),
             ("p2rank_top_score", "P2Rank top score"),
             ("plddt_mean", "mean pLDDT"),
             ("n_binders", "number of binders")]
d = MET.merge(
    DESC.rename(columns={"gene": "target"}).drop(columns=["family"],
                                                 errors="ignore"),
    on="target", how="left")
rows = []
for col, label in desc_cols:
    if col not in d.columns:
        print(f"  (missing column: {col})")
        continue
    for fam in ("GPCR", "SLC"):
        sub = d[d.family == fam].dropna(subset=[col])
        for out in ("AUROC", "EF1"):
            r, p = stats.pearsonr(sub[col], sub[out])
            rows.append((label, fam, out, r, p))
for fam in ("GPCR", "SLC"):
    sub = MET[MET.family == fam]
    for out in ("AUROC", "EF1"):
        r, p = stats.pearsonr(sub.ranking_accuracy, sub[out])
        rows.append(("ranking accuracy", fam, out, r, p))

if rows:
    qs = bh([r[4] for r in rows])
    report("number of exploratory correlations", len(rows), 28, "{:.0f}")
    ts = [(r, q) for r, q in zip(rows, qs) if r[0] != "ranking accuracy"]
    report("minimum uncorrected P, target-side",
           min(r[4] for r, _ in ts), 0.038, "{:.3f}")
    report("minimum q, target-side", min(q for _, q in ts), 0.512, "{:.3f}")
    report("target-side associations surviving BH",
           sum(1 for _, q in ts if q < 0.05), 0, "{:.0f}")
    for r, q in zip(rows, qs):
        if r[0] == "ranking accuracy" and r[1] == "GPCR" and r[2] == "AUROC":
            report("ranking accuracy -> AUROC (GPCR), r", r[3], 0.60, "{:.2f}")
            report("ranking accuracy -> AUROC (GPCR), P", r[4], 0.006, "{:.3f}")
            report("ranking accuracy -> AUROC (GPCR), q", q, 0.159, "{:.3f}")
        if r[0] == "ranking accuracy" and r[1] == "GPCR" and r[2] == "EF1":
            report("ranking accuracy -> EF@1% (GPCR), r", r[3], 0.44, "{:.2f}")
            report("ranking accuracy -> EF@1% (GPCR), P", r[4], 0.055, "{:.3f}")
        if r[0] == "ranking accuracy" and r[1] == "SLC" and r[2] == "EF1":
            report("ranking accuracy -> EF@1% (SLC), r", r[3], 0.27, "{:.2f}")
            report("ranking accuracy -> EF@1% (SLC), P", r[4], 0.21, "{:.2f}")
        if r[0] == "P2Rank top score" and r[2] == "EF1":
            report(f"P2Rank top score -> EF@1% ({r[1]}), r", r[3],
                   0.47 if r[1] == "GPCR" else -0.36, "{:.2f}", tol=0.02)
print()

# ================================================================ Results 5
print("Training-domain representation (Results 5, Fig. 5, Suppl. Table S8)")
tr = TRAIN.rename(columns={"Target": "target"})
gg = tr[tr.family == "GPCR"]
ss = tr[tr.family == "SLC"]
LB, CH = "pdb_ligand_bound_pre_cutoff", "chembl_compounds"
report("median ligand-bound PDB entries, GPCR", gg[LB].median(), 4, "{:.0f}")
report("median ligand-bound PDB entries, SLC", ss[LB].median(), 0, "{:.0f}")
report("ligand-bound PDB, Mann-Whitney P",
       stats.mannwhitneyu(gg[LB], ss[LB], alternative="two-sided").pvalue,
       9.5e-6, "{:.2e}", tol=5e-6)
report("median ChEMBL compounds, GPCR", gg[CH].median(), 1720, "{:.0f}", tol=1)
report("median ChEMBL compounds, SLC", ss[CH].median(), 2, "{:.0f}", tol=1)
report("ChEMBL compounds, Mann-Whitney P",
       stats.mannwhitneyu(gg[CH], ss[CH], alternative="two-sided").pvalue,
       1.9e-5, "{:.2e}", tol=1e-5)
report("SLC targets with no ligand-bound structure",
       int((ss[LB] == 0).sum()), 21, "{:.0f}")
report("GPCR targets with no ligand-bound structure",
       int((gg[LB] == 0).sum()), 4, "{:.0f}")
rho, p = stats.spearmanr(gg[LB], gg["EF@1%"])
report("ligand-bound PDB -> EF@1% (GPCR), rho", rho, 0.477, "{:.3f}")
report("ligand-bound PDB -> EF@1% (GPCR), P", p, 0.033, "{:.3f}")
rho, p = stats.spearmanr(gg["chembl_activities"], gg["EF@1%"])
report("ChEMBL activities -> EF@1% (GPCR), rho", rho, 0.535, "{:.3f}")
report("ChEMBL activities -> EF@1% (GPCR), P", p, 0.015, "{:.3f}")
print()

# ================================================================ Results 6
print("Score separation (Results 6, Fig. 6)")
r, p = stats.pearsonr(MET.Sep, MET.AUROC)
report("Sep -> AUROC, pooled (n = 43), r", r, 0.97, "{:.2f}")
report("Sep -> AUROC, pooled, R^2", r ** 2, 0.94, "{:.2f}")
report("GPCR mean Sep", G.Sep.mean(), 0.355, "{:.3f}")
report("SLC  mean Sep", S.Sep.mean(), 0.115, "{:.3f}")
report("GPCR Sep minimum", G.Sep.min(), -0.38, "{:.2f}")
report("GPCR Sep maximum", G.Sep.max(), 1.47, "{:.2f}")
report("SLC  Sep minimum", S.Sep.min(), -0.33, "{:.2f}")
report("SLC  Sep maximum", S.Sep.max(), 0.52, "{:.2f}")
report("Sep, Mann-Whitney P",
       stats.mannwhitneyu(G.Sep, S.Sep, alternative="two-sided").pvalue,
       0.13, "{:.2f}")
if "iptm" in SD.columns:
    gi = SD[(SD.family == "GPCR") & SD.iptm.notna()]
    si = SD[(SD.family == "SLC") & SD.iptm.notna()]
    report("GPCR mean ipTM", gi.iptm.mean(), 0.93, "{:.2f}")
    report("SLC  mean ipTM", si.iptm.mean(), 0.89, "{:.2f}")
    m = si.merge(MET[["target", "Sep"]], on="target", suffixes=("_sd", ""))
    r, p = stats.pearsonr(m.iptm, m.Sep)
    report("ipTM -> Sep (SLC), r", r, 0.26, "{:.2f}")
    report("ipTM -> Sep (SLC), P", p, 0.24, "{:.2f}")
print()

print("  Held-out evaluation (Results 6)")
# sep_heldout_target.csv stores the unconstrained linear fit, which returns
# negative predicted enrichment for eight targets. The manuscript uses an
# isotonic fit, recomputed here by the same leave-one-target-out procedure.
ht = HELD_T
pred = []
for i in range(len(ht)):
    train = ht.drop(ht.index[i])
    iso = IsotonicRegression(out_of_bounds="clip").fit(train.Sep, train.EF1_true)
    pred.append(float(iso.predict([ht.Sep.iloc[i]])[0]))
pred = np.asarray(pred)
r, p = stats.pearsonr(ht.EF1_true, pred)
report("leave-one-target-out (isotonic), r", r, 0.889, "{:.3f}")
report("leave-one-target-out (isotonic), R^2", r ** 2, 0.789, "{:.3f}")
report("leave-one-target-out (isotonic), MAE",
       np.mean(np.abs(ht.EF1_true - pred)), 4.66, "{:.2f}")
report("isotonic fit, negative predictions", int((pred < 0).sum()), 0, "{:.0f}")
report("linear fit, negative predictions",
       int((ht.EF1_pred < 0).sum()), 8, "{:.0f}")
print()

# ================================================================ Results 7
print("Model-generated interface consistency (Results 7, Fig. 7, Suppl. Table S9)")
cg = CC[CC.family == "GPCR"]
cs = CC[CC.family == "SLC"]
report("n targets with predicted complexes", len(CC), 24, "{:.0f}")
report("GPCR mean binder contact overlap (8 A)", cg.jb8.mean(), 0.596, "{:.3f}")
report("SLC  mean binder contact overlap (8 A)", cs.jb8.mean(), 0.629, "{:.3f}")
report("binder contact overlap, Mann-Whitney P",
       stats.mannwhitneyu(cg.jb8, cs.jb8, alternative="two-sided").pvalue,
       0.68, "{:.2f}")
report("GPCR mean centroid spread (A)", cg.centroid_spread.mean(), 4.51, "{:.2f}")
report("SLC  mean centroid spread (A)", cs.centroid_spread.mean(), 3.38, "{:.2f}")
report("centroid spread, Mann-Whitney P",
       stats.mannwhitneyu(cg.centroid_spread, cs.centroid_spread,
                          alternative="two-sided").pvalue, 0.48, "{:.2f}")

cc = CC.merge(MET[["target", "Sep", "EF1", "AUROC", "BEDROC80.5"]],
              on="target", how="left", suffixes=("", "_m"))
pairs, ps = [], []
for meas in ("jb8", "jr8", "centroid_spread"):
    for out in ("Sep", "EF1", "AUROC", "BEDROC80.5"):
        col = out + "_m" if out + "_m" in cc.columns else out
        rho, p = stats.spearmanr(cc[meas], cc[col])
        pairs.append((meas, out, rho, p))
        ps.append(p)
qs = bh(ps)
report("number of correlations examined", len(pairs), 12, "{:.0f}")
jb = [(x, q) for x, q in zip(pairs, qs) if x[0] == "jb8"]
report("binder overlap, minimum rho", min(x[2] for x, _ in jb), 0.50, "{:.2f}", tol=0.02)
report("binder overlap, maximum rho", max(x[2] for x, _ in jb), 0.59, "{:.2f}", tol=0.02)
report("binder overlap, maximum q", max(q for _, q in jb), 0.031, "{:.3f}", tol=0.005)
cs_ = [(x, q) for x, q in zip(pairs, qs) if x[0] == "centroid_spread"]
report("centroid spread, minimum q", min(q for _, q in cs_), 0.74, "{:.2f}", tol=0.05)
print()

# ================================================================ Results 8
print("Comparison of the two model outputs (Results 8, Fig. 8, Suppl. Table S7)")
g = CLS[CLS.Class == "GPCR"] if "Class" in CLS.columns else CLS[CLS.family == "GPCR"]
s = CLS[CLS.Class == "SLC"] if "Class" in CLS.columns else CLS[CLS.family == "SLC"]
report("GPCR mean AUROC, affinity value", g.AUROC_val.mean(), 0.64, "{:.2f}")
report("GPCR mean AUROC, classification prob", g.AUROC_prob.mean(), 0.70, "{:.2f}")
report("SLC  mean AUROC, affinity value", s.AUROC_val.mean(), 0.56, "{:.2f}")
report("SLC  mean AUROC, classification prob", s.AUROC_prob.mean(), 0.58, "{:.2f}")
report("GPCR two outputs, Wilcoxon P (exact)",
       stats.wilcoxon(g.AUROC_val, g.AUROC_prob, method="exact").pvalue,
       0.002, "{:.3f}")
report("SLC  two outputs, Wilcoxon P (exact)",
       stats.wilcoxon(s.AUROC_val, s.AUROC_prob, method="exact").pvalue,
       0.45, "{:.2f}")
report("GPCR vs SLC (affinity value), Mann-Whitney P",
       stats.mannwhitneyu(g.AUROC_val, s.AUROC_val,
                          alternative="two-sided").pvalue, 0.17, "{:.2f}")
report("GPCR vs SLC (classification prob), Mann-Whitney P",
       stats.mannwhitneyu(g.AUROC_prob, s.AUROC_prob,
                          alternative="two-sided").pvalue, 0.03, "{:.2f}")
for name, a, b in [("SLC18A2", 0.741, 0.626), ("SLC18B1", 0.648, 0.777)]:
    row = CLS[CLS.Target == name] if "Target" in CLS.columns else CLS[CLS.target == name]
    if len(row):
        report(f"  {name}, affinity value", float(row.AUROC_val.iloc[0]), a, "{:.3f}")
        report(f"  {name}, classification prob", float(row.AUROC_prob.iloc[0]), b, "{:.3f}")
print()

# ================================================================ summary
n_ok, n_all = sum(_checks), len(_checks)
print("=" * 78)
print(f"{n_ok} of {n_all} values match the manuscript within tolerance.")
if n_ok < n_all:
    print("Values marked DIFFERS above should be checked.")
print("=" * 78)
