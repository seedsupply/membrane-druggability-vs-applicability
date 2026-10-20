# membrane-druggability-vs-applicability

Aggregate data and analysis code for:

**Evaluating the applicability domain of an AI affinity predictor against standardized membrane-protein binding measurements**
Naoki Tarui, Thuy Duong Nguyen, Masaharu Nakayama
SEEDSUPPLY INC., Japan

---

## What this repository contains

Target-level aggregate metrics supporting every figure and table in the manuscript, together with the scripts that reproduce the reported statistics from those metrics.

**This repository does not contain per-compound data.** Individual compound structures and per-compound binding measurements from the Binder2030 platform are proprietary and available under a Data Use Agreement from the corresponding author for academic research; per-compound model predictions are available under the same agreement. A PubChem-overlap subset of the experimental data has been released separately (Tarui, N., Nakayama, M. & Nguyen, T. D. *SLAS Discov.* **39**, 100299, 2026).

The deposited files are sufficient to reproduce the reported downstream statistical
analyses and the associated quantitative figures and tables without access to the
proprietary compound-level measurements. The compound-level evaluation itself cannot
be independently repeated, because the underlying measurements and the per-compound
model predictions are proprietary.

---

## Data files

All files are one row per target.

### Experimental ligandability

| File | n | Contents | Used in |
|---|---|---|---|
| `data/slc_targets_296.csv` | 296 | Per-target screening outcome, structural fold, mechanistic class, subcellular location and structure-derived descriptors for all SLC transporters | Fig. 2; Results §2; Suppl. Table S1 |
| `data/gpcr_targets_295.csv` | 295 | Per-target screening outcome and CWxP toggle-switch status for all GPCRs | Fig. 2; Results §2; Suppl. Table S2 |

### Virtual-screening evaluation

| File | n | Contents | Used in |
|---|---|---|---|
| `data/per_target_metrics.csv` | 43 | Consolidated per-target metrics for the 20 GPCRs and 23 accessible-cavity SLC transporters evaluated: AUROC, EF at 0.5%, 1% and 5%, BEDROC at three α values, the separation index, per-target ranking accuracy, binder prevalence, the attainable EF maximum, the random expectation of BEDROC, and bootstrap 95% confidence intervals for every metric | Fig. 3, 6; Table 1; Suppl. Tables S4, S5 |
| `data/score_distributions.csv` | 43 | Predicted-score distribution statistics per target: binder and non-binder means, standard deviations, percentiles, the separation index, the binder top margin, and interface confidence (ipTM) | Results §6; Suppl. Table S4 |
| `data/clsprob_FINAL.csv` | 43 | AUROC, EF@1% and Sep computed from each of the two Boltz-2 scalar outputs (predicted affinity value; binary classification probability) | Fig. 8; Suppl. Table S7 |
| `data/slc_ranking_accuracy_82.csv` | 82 | Per-target ranking accuracy over the binder set only, for all SLC transporters with at least 10 binder pairs | Results §4; Methods |

### Analyses of candidate explanations

| File | n | Contents | Used in |
|---|---|---|---|
| `data/target_descriptors.csv` | 592 | Structure-derived descriptors for all GPCRs and SLC transporters: fpocket druggability score, pocket count and top pocket volume, P2Rank top score, mean pLDDT | Fig. 4; Suppl. Table S6 |
| `data/training_representation.csv` | 43 | Public-data proxies for representation in the model's training domain: PDB entries before the structure cutoff, the ligand-bound subset, ChEMBL compound and activity-record counts | Fig. 5; Suppl. Table S8 |
| `data/contact_consistency_24.csv` | 24 | Consistency of model-generated ligand placements at three distance cutoffs, with the spatial spread of ligand centroids | Fig. 7; Suppl. Table S9 |
| `data/sep_heldout_compound.csv` | 8,538 | Compound-split held-out evaluation of the separation index (200 repetitions per target) | Results §6 |
| `data/sep_heldout_target.csv` | 43 | Leave-one-target-out held-out evaluation of the separation index | Results §6 |

---

## Definitions

**Predicted score.** Boltz-2 returns `affinity_pred_value` on a scale where lower values indicate stronger predicted binding. Throughout, we use

```
predicted score = -affinity_pred_value
```

so that higher values indicate stronger predicted binding. Zero on this scale is a convention of the model's output and not a calibrated boundary between binding and non-binding.

**Separation index.**

```
Sep = (mean_binder - mean_nonbinder) / (SD_binder + SD_nonbinder)
```

with sample standard deviations (`ddof = 1`). Sep and AUROC are computed from the same predicted scores and the same experimental labels, so their close association is not independent evidence of predictive performance. Sep is a retrospective diagnostic, not a prospective indicator.

**Binder top margin.**

```
binder_top_margin = max(binder scores) - 95th percentile(non-binder scores)
```

**Enrichment factor.**

```
EF@X% = (h / n) / (k / N)
```

where *N* is the number of compounds screened for that target, *n* the number of binders, *k* = ceil(*N* × X%) the size of the top-ranked subset, and *h* the number of binders within it. Because *k* is finite and *n* varies across targets, the attainable maximum is bounded by binder prevalence; the per-target maximum is given in `per_target_metrics.csv`.

**BEDROC.** Computed following Truchon and Bayly (2007), Eq. 36. α = 80.5 is the primary setting; α = 20 and α = 321.9 are reported as sensitivity analyses. The attainable range and random expectation depend on α, library size and binder prevalence; the per-target random expectation at the primary α is given in `per_target_metrics.csv`.

**Confidence intervals.** Two-sided 95% percentile intervals from 2,000 bootstrap resamples of the screened compounds within each target. These describe the precision of each target-level estimate and are not used to classify targets.

**Hit recovery.** The fraction of targets for which the affinity-selection screen returned at least one compound meeting the binder criterion. This is a per-target measure of experimental ligandability and is independent of any model prediction.

**Ranking accuracy.** Pearson's *r* between predicted score and experimental pKd over the experimental binder set only (minimum 10 binder pairs), where pKd = −log₁₀(Kd) with Kd in molar units. It derives from a separate prediction run to the virtual-screening metrics, which score binders against the full non-binder library.

---

## Reproducing the reported statistics

```bash
pip install -r requirements.txt
python scripts/reproduce_statistics.py
python scripts/verify_ligandability_statistics.py
```

`reproduce_statistics.py` recomputes the correlations, group comparisons and threshold analyses reported for the virtual-screening evaluation from the CSV files in `data/`, and prints each one alongside the value stated in the manuscript so that any discrepancy is immediately visible.

`verify_ligandability_statistics.py` does the same for the experimental ligandability analysis: hit recovery by protein family and by mechanistic class, the Fisher exact test comparing GPCRs with accessible-cavity transporters, and the chi-square test across mechanistic classes.

Note that these scripts reproduce the *statistics*, not the underlying predictions. Regenerating the per-target metrics themselves requires the per-compound experimental data and model predictions, which are available under the Data Use Agreement described above.

---

## Calculation conditions

`CALC_CONDITIONS.md` records the conditions under which the per-target metrics were generated, including the prediction runs used, the sign convention, the assay concentration range, and analyses that were performed and subsequently found not to hold. It is included for transparency about how the reported values were arrived at.

---

## License

- **Code** (`scripts/`): MIT License — see `LICENSE`
- **Data** (`data/`): Creative Commons Attribution 4.0 International (CC BY 4.0) — see `data/LICENSE`

No registration, login, or data-use agreement is required to access any material in this repository.

---

## Related resources

- Boltz-2: https://github.com/jwohlwend/boltz
- AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk
- fpocket: https://github.com/Discngine/fpocket
- P2Rank: https://github.com/rdk/p2rank
- Binder2030 (PubChem-overlap subset): *SLAS Discov.* **39**, 100299 (2026)

---

## Contact

Correspondence: naoki.tarui@seedsupply.co.jp

Requests for access to per-compound data under a Data Use Agreement should be directed to the corresponding author.
