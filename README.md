# membrane-druggability-vs-applicability

Aggregate data and analysis code for:

**Defining the applicability domain of AI-based affinity prediction using experimentally standardized membrane-protein benchmarks**
Naoki Tarui, Thuy Duong Nguyen, Masaharu Nakayama
SEEDSUPPLY INC., Japan

---

## What this repository contains

Target-level aggregate metrics supporting all figures and tables in the manuscript, together with the scripts that reproduce the reported statistics from those metrics.

**This repository does not contain per-compound data.** Individual compound structures and per-compound binding measurements from the Binder2030 platform are proprietary and available under a Data Use Agreement from the corresponding author for academic research; per-compound model predictions are available under the same agreement. A PubChem-overlap subset of the experimental data has been released separately (Tarui, N., Nakayama, M. & Nguyen, T. D. *SLAS Discov.* **39**, 100299, 2026).

---

## Data files

All files are one row per target.

### Experimental druggability

| File | n | Contents | Used in |
|---|---|---|---|
| `data/druggability_by_fold_592.csv` | 592 | Per-target screening outcome (hit / no hit) for 296 GPCRs and 296 SLC transporters, with protein family and, for SLC transporters, structural fold class | Fig. 2; Results (hit recovery, Fisher, chi-square) |

### AI virtual-screening evaluation

| File | n | Contents | Used in |
|---|---|---|---|
| `data/gpcr_VS_FINAL_20.csv` | 20 | GPCR virtual-screening metrics: AUROC, EF at 0.5% and 1%, Sep, ranking accuracy, outcome | Fig. 3, 4, 5; Table 1; Suppl. Table S3 |
| `data/gpcr_tail_FINAL_20.csv` | 20 | GPCR predicted-score distribution statistics: binder/non-binder means, SDs, percentiles, Sep, binder top margin | Fig. 5, 6; Suppl. Table S3 |
| `data/slc_VS_FINAL_23.csv` | 23 | As above for accessible-cavity SLC transporters, with distribution statistics included | Fig. 3, 4, 5, 6; Table 1; Suppl. Table S3 |
| `data/clsprob_FINAL.csv` | 43 | AUROC, EF@1% and Sep computed from each of the two Boltz-2 scalar outputs (predicted affinity value; binary classification probability) | Fig. 7; Table 1; Suppl. Table S5 |
| `data/sep_decomp_gpcr.csv` | 20 | Decomposition of Sep into its components for GPCRs | Fig. 6a |
| `data/sep_decomp_slc.csv` | 23 | Decomposition of Sep into its components for accessible-cavity transporters | Fig. 6a |
| `data/slc_ranking_accuracy_82.csv` | 82 | Per-target ranking accuracy (Pearson's *r* between predicted score and experimental pKd) over the binder set only, for all SLC transporters with at least 10 binder pairs | Results (individual cases); Methods |
| `data/mechanism_slc.csv` | 23 | Per-target interface predicted TM-score (ipTM) and confidence alongside Sep and binder scores | Fig. 6b |

Structural descriptors used in the target-side analysis (fpocket druggability score and pocket count, P2Rank top score, pLDDT) are provided in the virtual-screening files above and, for the full target set, as Supplementary Tables S2 and S7 with the manuscript.

---

## Definitions

**Predicted score.** Boltz-2 returns `affinity_pred_value` on a scale where lower values indicate stronger predicted binding. Throughout, we use

```
predicted score = -affinity_pred_value
```

so that higher values indicate stronger predicted binding. A target with `binder_mean > 0` is one where experimental binders are placed on the binding side of the model's own scale.

**Separation index.**

```
Sep = (mean_binder - mean_nonbinder) / (SD_binder + SD_nonbinder)
```

with sample standard deviations (`ddof = 1`).

**Binder top margin.**

```
binder_top_margin = max(binder scores) - 95th percentile(non-binder scores)
```

**Virtual-screening success.** AUROC >= 0.7.

**Hit recovery.** The fraction of targets for which the affinity-selection screen returned at least one compound meeting the binder criterion. This is a per-target measure of experimental druggability and is independent of any model prediction.

**Ranking accuracy.** Per-target ranking accuracy is computed over the experimental binder set only (minimum 10 binder pairs), and is therefore derived from a separate prediction run to the virtual-screening metrics, which score binders against the full non-binder library. `slc_ranking_accuracy_82.csv` reports this quantity for all 82 SLC transporters meeting the pair-count threshold; the 23 transporters evaluated by virtual screening are a subset of these.

**Enrichment factor.**

```
EF@X% = (actives in top X% of ranked compounds) / (actives expected at random)
```

---

## Reproducing the reported statistics

```bash
pip install -r requirements.txt
python scripts/reproduce_statistics.py
python scripts/verify_druggability_statistics.py
```

`reproduce_statistics.py` recomputes every correlation, group comparison and effect size reported for the virtual-screening analysis from the CSV files in `data/`, and prints them alongside the values stated in the text so that any discrepancy is visible.

`verify_druggability_statistics.py` does the same for the experimental druggability analysis: hit recovery by protein family and by structural fold, the Fisher exact test comparing GPCRs with accessible-cavity transporters, and the chi-square test across fold classes.

Note that these scripts reproduce the *statistics*, not the underlying predictions. Regenerating the per-target metrics themselves requires the per-compound experimental data and model predictions, which are available under the Data Use Agreement described above.

---

## Calculation conditions

`CALC_CONDITIONS.md` records the conditions under which the per-target metrics were generated, including the prediction runs used, the sign convention, and analyses that were performed and subsequently found not to hold. It is included for transparency about how the reported values were arrived at.

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

---

## Contact

Correspondence: naoki.tarui@seedsupply.co.jp

Requests for access to per-compound data under a Data Use Agreement should be directed to the corresponding author.
