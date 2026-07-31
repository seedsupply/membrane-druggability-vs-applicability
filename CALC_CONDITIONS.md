# Calculation conditions

This file records the conventions under which the per-target metrics in `data/`
were generated, the provenance of each file, and analyses that were performed
and subsequently found not to be supported. It is included so that the reported
values can be interpreted without access to the per-compound data.

Values verified 2026-07-13; unchanged since.

---

## 1. Sign convention and definitions

**Predicted score.** Boltz-2 returns `affinity_pred_value` on a scale where
lower values indicate stronger predicted binding. All files in `data/` use

```
predicted score = -affinity_pred_value
```

so that higher values indicate stronger predicted binding. Where a column is
named `pred_pKd` or `predicted_pKd`, it holds this sign-flipped quantity, not
an experimentally measured pKd.

**Separation index.**

```
Sep = (binder_mean - nonbinder_mean) / (binder_SD + nonbinder_SD)
```

Standard deviations are sample standard deviations (`ddof = 1`).

**Binder top margin.**

```
binder_top_margin = max(binder scores) - 95th percentile(non-binder scores)
```

**Virtual-screening success.** A target is counted as successful when
AUROC >= 0.7.

**Multi-binder targets.** Targets with more than one experimental binder are
not excluded from any analysis.

---

## 2. Provenance of the per-target metrics

**GPCR virtual screening (20 targets).** Metrics were computed from the
full-library prediction runs in which each target's binder set was scored
together with the 1,200-compound GPCR library in a single batch. For one
target (CB1R) the corresponding run was stored separately and was used in
place of the batch output; the protocol was otherwise identical.

Prediction runs in which binders were scored on their own, without the
non-binder library, were **not** used for any virtual-screening metric. Those
runs give different values because the non-binder distribution against which
binders are ranked is absent. Files derived from them do not appear in `data/`.

Outputs: `gpcr_tail_FINAL_20.csv`, `gpcr_VS_FINAL_20.csv`

**SLC virtual screening (23 targets).** Metrics were computed from the
equivalent full-library runs against the SLC compound library. Column naming
in the source runs was not uniform (`pred_pKd` and `predicted_pKd` both
occur); the published files use a single consistent naming.

Output: `slc_VS_FINAL_23.csv`

**SLC ranking accuracy (82 targets).** Per-target ranking accuracy is a
different quantity with a different source. It is the Pearson correlation
between predicted score and experimental pKd computed over the experimental
binder set only, for every SLC transporter with at least 10 binder pairs
(1,391 pairs across 82 targets). It therefore derives from the binder-only
prediction runs, not from the full-library runs used for virtual screening.
The 23 transporters evaluated by virtual screening are a subset of these 82.

Output: `slc_ranking_accuracy_82.csv`

Because these two families of metrics come from different runs, per-target
ranking accuracy and per-target virtual-screening metrics should not be
compared value-for-value across files.

**Experimental druggability (592 targets).** Screening outcomes and fold
classifications are independent of any model prediction and derive directly
from the affinity-selection screen.

The SLC set comprises the 296 transporters that could be assigned to one of the
four mechanistic classes. One additional transporter (SLC53A1/XPR1) was screened
but could not be assigned to a fold, and is therefore not included in the
classification or in `druggability_by_fold_592.csv`.

One entry required correction before the set was finalised. UniProt accession
Q9NQQ7 initially appeared under two gene symbols, SLC35C2 and SLC35H1. These are
the same gene; the duplicate row was removed and the target is listed under its
HGNC-approved symbol, SLC35H1.

Output: `druggability_by_fold_592.csv`

---

## 3. Confirmed results

Virtual-screening success (AUROC >= 0.7):

| Family | Successful | Evaluated |
|---|---|---|
| GPCR | 7 | 20 |
| SLC (accessible cavity) | 2 | 23 |

Correlations between the separation index and screening performance:

| Relationship | GPCR | SLC |
|---|---|---|
| Sep to AUROC | 0.972 | 0.991 |
| Sep to EF@1% | 0.922 | 0.573 |
| binder_top_margin to EF@1% | 0.798 | 0.625 |

The weaker Sep-to-EF@1% correlation among SLC transporters reflects the
compressed range of Sep in that family rather than a different underlying
relationship: restricting GPCRs to the same range (Sep < 0.55) reduces their
correlation to 0.655.

Representative per-target ranking accuracies (binder-only runs): SLC2A7
+0.729, SLC18A2 +0.709, SLC6A4 -0.216.

---

## 4. Analyses performed and not supported

The following analyses were carried out and are reported here because they
bound what the data support. None of them is used in the manuscript.

**Lower-tail hypothesis.** The 5th percentile of the non-binder score
distribution (`nb_5pct`) was tested as a predictor of EF@1%. It was not
significant in either family. The hypothesis that screening failure is driven
by the lower tail of the non-binder distribution is therefore not supported.

**Per-target ranking accuracy as a predictor of screening success.** Ranking
accuracy over the binder set does not predict virtual-screening performance.
The separation between binder and non-binder score distributions does. These
are distinct quantities and only the latter is used in the manuscript.

---

## 5. Corrections made during analysis

Recorded for completeness, since earlier internal versions of these numbers
circulated before the files in `data/` were finalised.

- An early note recorded OX1R as a virtual-screening failure with AUROC 0.088.
  This was an error: the correct value is AUROC 0.912 with EF 41.7, a success.
  The figure 0.088 belonged to a different target and a different metric.
- An early set of Sep values was computed from binder-only prediction runs and
  did not match the full-library runs (for example ADRB1 1.501 rather than the
  correct 1.466). The files in `data/` are computed exclusively from the
  full-library runs. Target rankings and all reported conclusions are
  unaffected by this correction.
