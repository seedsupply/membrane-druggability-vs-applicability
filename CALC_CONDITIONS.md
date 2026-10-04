# Calculation conditions

This file records the conventions under which the per-target metrics in `data/`
were generated, the provenance of each file, and analyses that were performed
and subsequently found not to be supported. It is included so that the reported
values can be interpreted without access to the per-compound data.

Values verified 2026-10-02.

---

## 1. Experimental conventions

**Binder criterion.** Compounds retained by target-expressing membranes but not
by mock-transfected controls are carried into a concentration-dependent binding
analysis using a seven-point concentration series spanning 7–30,000 nM, with Kd
determined by nonlinear regression. A compound is designated a binder when a Kd
of 10 μM or lower is determined.

**Non-binder criterion.** Compounds that do not yield concentration-dependent
binding under this analysis, including those retained in the initial screen
without concentration dependence, are designated non-binders for that target.
This records a measured outcome rather than an absence of data, but it is
bounded by the concentration range and conditions tested and does not exclude
binding weaker than the threshold applied.

**Screening libraries.** Both family libraries were assembled by the same
procedure, from compounds identified as binders of targets within the same
family. The GPCR library comprises 1,200 compounds drawn from the binders of
100 GPCRs; the SLC library 1,371 compounds drawn from the binders of 160 SLC
transporters. Each compound was measured against the targets screened during
library construction and did not show concentration-dependent binding to targets
other than those for which it was identified as a binder. Selectivity is
therefore defined with respect to that target set rather than to the family as a
whole. The targets evaluated by virtual screening are contained within this set.

---

## 2. Sign convention and definitions

**Predicted score.** Boltz-2 returns `affinity_pred_value` on a scale where
lower values indicate stronger predicted binding. All files in `data/` use

```
predicted score = -affinity_pred_value
```

so that higher values indicate stronger predicted binding. **Zero on this scale
is a convention of the model's output and not a calibrated boundary between
binding and non-binding.** Where a column is named `pred_pKd` or `predicted_pKd`
in a source run, it holds this sign-flipped quantity, not an experimentally
measured pKd.

**Separation index.**

```
Sep = (binder_mean - nonbinder_mean) / (binder_SD + nonbinder_SD)
```

Standard deviations are sample standard deviations (`ddof = 1`).

**Binder top margin.**

```
binder_top_margin = max(binder scores) - 95th percentile(non-binder scores)
```

**Enrichment factor.**

```
EF@X% = (h / n) / (k / N),   k = ceil(N x X%)
```

The attainable maximum is bounded by binder prevalence, which ranges from 0.7%
to 20.6% across the 43 evaluated targets. Per-target prevalence and attainable
maxima are given in `per_target_metrics.csv`.

**BEDROC.** Truchon and Bayly (2007), Eq. 36. Primary setting α = 80.5;
sensitivity analyses at α = 20 and α = 321.9. The primary α was chosen because
the principal metrics of this study concern early enrichment, so a weighting
that emphasises early recognition is the appropriate companion to EF@0.5% and
EF@1%. The attainable range and random expectation depend on α, library size
and prevalence; the per-target random expectation at the primary α is given in
`per_target_metrics.csv`.

**Confidence intervals.** Two-sided 95% percentile intervals from 2,000
bootstrap resamples of the screened compounds **within each target**, without
stratification by label. Resamples in which one class was absent were discarded.
These intervals describe the precision of each target-level estimate and were
not used to classify targets.

**Interface confidence.** The mean ipTM across each target's predicted
complexes.

**Ties.** Ties in predicted score are assigned average ranks throughout.

---

## 3. The AUROC = 0.7 threshold

The original analysis classified a target as a virtual-screening success when
AUROC reached 0.7. This threshold was used without external justification and is
**not used as the primary result** in the published analysis. It is retained in
`per_target_metrics.csv` and in Supplementary Tables S4 and S7 for comparability
with the original submission.

Scanning cutoffs from 0.55 to 0.85 in steps of 0.01 shows that the Fisher *P*
value for the family comparison is not a monotonic function of the threshold: it
crosses the 0.05 level three times, reaching 0.039 at a cutoff of 0.68,
returning to 0.148 at 0.69, and falling to 0.016 from 0.78 upward. Adjacent
cutoffs differing by 0.01 therefore yield opposite conclusions about statistical
significance.

---

## 4. Provenance of the per-target metrics

**Virtual screening (43 targets).** Metrics were computed from the full-library
prediction runs in which each target's binder set was scored together with the
family library in a single batch. For one target (CB1R) the corresponding run
was stored separately and was used in place of the batch output; the protocol
was otherwise identical.

Prediction runs in which binders were scored on their own, without the
non-binder library, were **not** used for any virtual-screening metric. Those
runs give different values because the non-binder distribution against which
binders are ranked is absent.

Outputs: `per_target_metrics.csv`, `score_distributions.csv`

**Ranking accuracy (82 SLC targets).** Per-target ranking accuracy is a
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

**Experimental ligandability (591 targets).** Screening outcomes and structural
classifications are independent of any model prediction and derive directly from
the affinity-selection screen.

The SLC set comprises the 296 transporters that could be assigned to one of the
four mechanistic classes. One additional transporter (SLC53A1/XPR1) was screened
but could not be assigned to a fold, and is therefore not included.

One SLC entry required correction before the set was finalised. UniProt
accession Q9NQQ7 initially appeared under two gene symbols, SLC35C2 and
SLC35H1. These are the same gene; the duplicate row was removed and the target
is listed under its HGNC-approved symbol, SLC35H1.

One GPCR entry (GAST) was removed on re-checking, giving 295 rather than 296
GPCRs. All affected values were updated.

Outputs: `slc_targets_296.csv`, `gpcr_targets_295.csv`

**Training-domain proxies (43 targets).** The composition of the processed
Boltz-2 training set is not public. Representation was quantified using proxies
computable from the same public sources: PDB entries released before the
structure cutoff of 1 June 2023, retrieved through the RCSB Search API by
UniProt accession, with those containing a non-polymer entity of at least 200 Da
counted separately as ligand-bound entries; and ChEMBL compound and
activity-record counts restricted to records with Ki, Kd or IC50. **ChEMBL
counts reflect the database at the time of analysis rather than at the training
cutoff.**

Output: `training_representation.csv`

**Interface consistency (24 targets).** Boltz-2 does not output an explicit
pocket, so the predicted ligand placements were used to characterise the
model-generated interface. Contact residues are Cα atoms within the stated
cutoff of any ligand heavy atom; consistency is the mean pairwise Jaccard
overlap of contact-residue sets across binders, with the same quantity computed
across a random sample of 300 non-binders and their ratio taken. Cutoffs of 6, 8
and 10 Å were evaluated; 8 Å is used in the main text. Predicted complexes were
available for all 20 GPCRs and for four SLC transporters of the LeuT fold.

**This analysis quantifies the consistency of model-generated placements, not
their correctness.** Recovery of experimentally supported contact residues could
not be assessed because 21 of the 23 SLC transporters have no ligand-bound
structure.

Output: `contact_consistency_24.csv`

---

## 5. Multiple-comparison correction

Benjamini–Hochberg correction is applied across the exploratory set specified in
each analysis.

**Categorical comparisons.** Fisher's exact test for 2 x 2 tables; the chi-square
test of independence for the four-class comparison, where the minimum expected
count is 5.3 and the usual requirement is satisfied.

**Target-side descriptors: 28 correlations.** Six descriptors (fpocket
druggability score, fpocket pocket count, top pocket volume, P2Rank top score,
mean pLDDT, number of experimental binders) tested against AUROC and EF@1% in
each family, together with per-target ranking accuracy tested against the same
two metrics in each family. No target-side association survives correction
(minimum uncorrected *P* = 0.038, minimum *q* = 0.512).

**Interface consistency: 12 correlations.** Three measures (binder contact
overlap, binder-to-non-binder ratio, centroid spread) against Sep, EF@1%, AUROC
and BEDROC.

---

## 6. Analyses performed and not supported

The following analyses were carried out and are reported here because they bound
what the data support.

**The separation index as a prospective indicator.** Sep tracks screening
performance closely, but it is computed from the same predicted scores and the
same experimental labels used to assess that performance. Held-out evaluation in
two designs (compound splitting within targets; leave-one-target-out) showed
that the relationship between separation and enrichment is stable across
targets, but **both designs require experimentally defined binders for the
target being evaluated**. Sep is therefore a retrospective diagnostic of where
the predictor performed, not a prospective indicator of where it will perform.

**Target-side descriptors as predictors of screening performance.** None of the
six descriptors examined showed a consistent association with screening
performance in either family. The one descriptor reaching nominal significance
before correction, the P2Rank top score against EF@1% in GPCRs (*r* = +0.47),
was inversely correlated with the same metric in SLC transporters (*r* = −0.36),
so its sign is not consistent across families.

**Training-domain representation as a within-family indicator.** The two
families differ markedly on every proxy, but within families the proxies did not
consistently track performance. Among GPCRs the nominal associations were
sensitive to two targets with essentially no representation (GPR85, GPR37L1);
after excluding them, none of the nine associations examined reached
significance. Among SLC transporters, ChEMBL representation spans three orders
of magnitude with no association to any performance metric. Family and
representation are themselves strongly correlated (*r* = −0.64), so a joint
regression cannot resolve their contributions. **Training-domain representation
remains a plausible contributor to the family-level difference and is not
excluded.**

**Lower-tail hypothesis.** The 5th percentile of the non-binder score
distribution (`nb_5pct`) was tested as a predictor of EF@1%. It was not
significant in either family.

**Per-target ranking accuracy as a predictor of screening performance.** Ranking
accuracy over the binder set was not associated with enrichment in either family
(GPCRs *r* = 0.44, *P* = 0.055; SLC transporters *r* = 0.27, *P* = 0.21). Its
association with AUROC among GPCRs (*r* = 0.60, *P* = 0.006) did not survive
correction across the 28 exploratory correlations (*q* = 0.159).

---

## 7. Corrections made during analysis

Recorded for completeness, since earlier internal versions of these numbers
circulated before the files in `data/` were finalised.

**Transcription errors found on re-verification.** Three errors were found when
every supplementary table was re-checked against the underlying data: the binder
count for SLC18B1 (20 → 279), the binder count for CB1R (13 → 12), and the
non-binder count for SLC6A12 (1,343 → 1,361) with its corresponding enrichment
values. None affects any reported conclusion; the virtual-screening analyses
were run on the correct underlying data in every case.

**Bootstrap method.** An intermediate version of the confidence intervals used
stratified resampling of binders and non-binders separately. The published
intervals use simple resampling of the screened compounds within each target,
without stratification. The counts of targets whose lower bound exceeds the
random expectation differ between the two methods; the published values are from
the simple resampling.

**BEDROC weighting.** An intermediate version used α = 160.9 as the primary
setting. This value was not among the three specified for the analysis and gave
the smallest *P* value of those tried; it was withdrawn. The published primary
setting is α = 80.5.

**OX1R.** An early note recorded OX1R as a virtual-screening failure with AUROC
0.088. This was an error: the correct value is AUROC 0.912 with EF@1% 41.7. The
figure 0.088 belonged to a different target and a different metric.

**Early Sep values.** An early set of Sep values was computed from binder-only
prediction runs and did not match the full-library runs (for example ADRB1 1.501
rather than the correct 1.466). The files in `data/` are computed exclusively
from the full-library runs.

---

## 8. Known limitations of the deposited data

**Predictions were run without a fixed random seed.** Target-level estimates
therefore carry run-to-run variation that was not characterised systematically.
This is one reason the published analysis does not assign target-level outcomes
on the basis of a fixed threshold.

**The two libraries derive from different source target sets** (100 GPCRs and
160 SLC transporters), so their compositions are not matched between families.
Absolute enrichment values would shift with library size and composition.

**Interface consistency covers 24 of the 43 targets**, limited by the
availability of predicted complexes.
