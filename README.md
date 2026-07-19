[README.md](https://github.com/user-attachments/files/30165400/README.md)
# membrane-druggability-vs-applicability

Aggregate data and analysis code for:

**Experimentally druggable membrane proteins can be invisible to AI virtual screening**
Naoki Tarui, Thuy Duong Nguyen, Masaharu Nakayama
SEEDSUPPLY INC., Japan

---

## What this repository contains

Target-level aggregate metrics supporting all figures and tables in the manuscript, together with the scripts that reproduce the reported statistics from those metrics.

**This repository does not contain per-compound data.** Individual compound structures and per-compound binding measurements from the Binder2030 platform are proprietary and available under a Data Use Agreement from the corresponding author for academic research; per-compound model predictions are available under the same agreement. A PubChem-overlap subset of the experimental data has been released separately (Tarui, N., Nakayama, M. & Nguyen, T. D. *SLAS Discov.* **39**, 100299, 2026).

---

## Data files

All files are one row per target.

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

Experimental druggability data (hit recovery for 296 GPCRs and 296 SLC transporters, with fold classification and structural descriptors) are provided as Supplementary Tables S1, S2 and S7 with the manuscript.

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

**Ranking accuracy.** Per-target ranking accuracy is computed over the experimental binder set only (minimum 10 binder pairs), and is therefore derived from a separate prediction run to the virtual-screening metrics, which score binders against the full non-binder library. `slc_ranking_accuracy_82.csv` reports this quantity for all 82 SLC transporters meeting the pair-count threshold; the 23 transporters evaluated by virtual screening are a subset of these.


**Enrichment factor.**

```
EF@X% = (actives in top X% of ranked compounds) / (actives expected at random)
```

---

## Reproducing the reported statistics

```bash
pip install pandas numpy scipy
python scripts/reproduce_statistics.py
```

The script recomputes every correlation, group comparison and effect size reported in the manuscript from the CSV files in `data/`, and prints them alongside the values stated in the text so that any discrepancy is visible.

Note that the script reproduces the *statistics*, not the underlying predictions. Regenerating the per-target metrics themselves requires the per-compound experimental data and model predictions, which are available under the Data Use Agreement described above.

---

## Calculation conditions

`CALC_CONDITIONS.md` records the conditions under which the per-target metrics were generated, including the prediction runs used, the sign convention, and analyses that were performed and subsequently found not to hold. It is included for transparency about how the reported values were arrived at.

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
