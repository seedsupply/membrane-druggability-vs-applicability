#!/usr/bin/env python3
"""
Recompute the experimental-ligandability statistics reported in the manuscript
from the per-target screening outcomes in ../data/, and print each value next to
the one stated in the text.

These analyses are independent of any model prediction: they describe which
targets yielded at least one validated binder in the affinity-selection screen.

Usage:
    python scripts/verify_ligandability_statistics.py

Requires: pandas, numpy, scipy
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path(__file__).resolve().parent.parent / "data"

SLC = pd.read_csv(DATA / "slc_targets_296.csv")
GPCR = pd.read_csv(DATA / "gpcr_targets_295.csv")

_checks = []


def report(label, computed, stated, fmt="{:.3f}", tol=0.01):
    ok = stated is None or abs(computed - stated) <= tol
    _checks.append(ok)
    c = fmt.format(computed)
    s = "—" if stated is None else fmt.format(stated)
    flag = "" if ok else "   <-- DIFFERS"
    print(f"  {label:<54} computed {c:>10}   manuscript {s:>10}{flag}")


def wilson(k, n, z=1.959963984540054):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def chi_square(table):
    """Chi-square test of independence, with the minimum expected count."""
    chi2, p, dof, exp = stats.chi2_contingency(table)
    return chi2, p, dof, exp.min()


hit_col_s = "hit" if "hit" in SLC.columns else SLC.columns[-1]
hit_col_g = "hit" if "hit" in GPCR.columns else GPCR.columns[-1]


def as_bool(series):
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(("yes", "1", "true"))


SLC["_hit"] = as_bool(SLC[hit_col_s])
GPCR["_hit"] = as_bool(GPCR[hit_col_g])

print(__doc__.split("Usage")[0].strip())
print()
print("Datasets loaded:")
print(f"  GPCR targets                n = {len(GPCR)}")
print(f"  SLC transporters            n = {len(SLC)}")
print()

# ---------------------------------------------------------------- GPCRs
print("GPCR hit recovery (Results 2, Fig. 2)")
k, n = int(GPCR._hit.sum()), len(GPCR)
report("GPCR targets screened", n, 295, "{:.0f}")
report("GPCR targets with at least one binder", k, 212, "{:.0f}")
report("GPCR hit recovery (%)", 100 * k / n, 71.9, "{:.1f}", tol=0.05)
lo, hi = wilson(k, n)
report("  Wilson 95% CI lower (%)", 100 * lo, 66.5, "{:.1f}", tol=0.1)
report("  Wilson 95% CI upper (%)", 100 * hi, 76.7, "{:.1f}", tol=0.1)

cw = next((c for c in GPCR.columns if "cwxp" in c.lower()), None)
if cw:
    pos = GPCR[GPCR[cw].astype(str).str.contains(r"\(\+\)|present|yes", case=False,
                                                 regex=True, na=False)]
    neg = GPCR[GPCR[cw].astype(str).str.contains(r"\(-\)|absent|no", case=False,
                                                 regex=True, na=False)]
    if len(pos) and len(neg):
        kp, np_ = int(pos._hit.sum()), len(pos)
        kn, nn = int(neg._hit.sum()), len(neg)
        report("  CWxP present, hit recovery (%)", 100 * kp / np_, 73.0, "{:.1f}", tol=0.1)
        report("  CWxP absent, hit recovery (%)", 100 * kn / nn, 71.2, "{:.1f}", tol=0.1)
        p = stats.fisher_exact([[kp, np_ - kp], [kn, nn - kn]])[1]
        report("  CWxP comparison, Fisher P", p, 0.79, "{:.2f}")
print()

# ---------------------------------------------------------------- SLCs
print("SLC hit recovery by mechanistic class (Results 2, Fig. 2, Suppl. Table S3)")
mc = next((c for c in SLC.columns if "mechanistic" in c.lower()), None)
if mc is None:
    print("  (no mechanistic-class column found)")
else:
    order = ["Accessible-cavity", "Compact gated", "Elevator", "Ion-channel-like"]
    rows, table = [], []
    for cls, grp in SLC.groupby(mc):
        kk, nn = int(grp._hit.sum()), len(grp)
        rows.append((cls, nn, kk, 100 * kk / nn))
        table.append([kk, nn - kk])
    rows.sort(key=lambda r: -r[3])
    for cls, nn, kk, pct in rows:
        print(f"    {str(cls):<34} {kk:3d}/{nn:3d}   {pct:5.1f}%")
    print()
    exp = {"accessible": (181, 129, 71.3), "compact": (66, 30, 45.5),
           "elevator": (36, 15, 41.7), "channel": (13, 2, 15.4)}
    for cls, nn, kk, pct in rows:
        key = next((k for k in exp if k in str(cls).lower()), None)
        if key:
            en, ek, ep = exp[key]
            report(f"  {str(cls)[:28]:<28} n", nn, en, "{:.0f}")
            report(f"  {str(cls)[:28]:<28} hits", kk, ek, "{:.0f}")
            report(f"  {str(cls)[:28]:<28} recovery (%)", pct, ep, "{:.1f}", tol=0.1)
    chi2, p, dof, min_exp = chi_square(table)
    report("  across mechanistic classes, chi-square", chi2, 31.05, "{:.2f}", tol=0.05)
    report("  across mechanistic classes, df", dof, 3, "{:.0f}")
    report("  across mechanistic classes, P", p, 8.3e-7, "{:.2e}", tol=5e-8)
    report("  minimum expected count", min_exp, 5.3, "{:.1f}", tol=0.05)
print()

# ------------------------------------------------- GPCR vs accessible-cavity
print("GPCRs versus accessible-cavity SLC transporters (Results 2)")
if mc is not None:
    acc = SLC[SLC[mc].astype(str).str.contains("accessible", case=False, na=False)]
    ka, na = int(acc._hit.sum()), len(acc)
    report("accessible-cavity targets", na, 181, "{:.0f}")
    report("accessible-cavity with at least one binder", ka, 129, "{:.0f}")
    report("accessible-cavity hit recovery (%)", 100 * ka / na, 71.3, "{:.1f}", tol=0.05)
    lo, hi = wilson(ka, na)
    report("  Wilson 95% CI lower (%)", 100 * lo, 64.3, "{:.1f}", tol=0.1)
    report("  Wilson 95% CI upper (%)", 100 * hi, 77.4, "{:.1f}", tol=0.1)
    odds, p = stats.fisher_exact([[k, n - k], [ka, na - ka]])
    report("GPCR vs accessible-cavity, Fisher P", p, 0.92, "{:.2f}")
    report("GPCR vs accessible-cavity, odds ratio", odds, 1.03, "{:.2f}", tol=0.02)
    print()
    print("  A non-significant test of this kind does not establish that the two")
    print("  values are the same, and the manuscript does not treat it as doing so.")
print()

n_ok, n_all = sum(_checks), len(_checks)
print("=" * 78)
print(f"{n_ok} of {n_all} values match the manuscript within tolerance.")
if n_ok < n_all:
    print("Values marked DIFFERS above should be checked.")
print("=" * 78)
