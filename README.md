# Binder2030 × Boltz-2 Benchmark: GPCR Virtual Screening

This repository contains analysis code and Boltz-2 prediction results for the following manuscript:

**"What Current AI Affinity Predictors Have and Have Not Learned: A Native-State Benchmark of Boltz-2 Virtual Screening across GPCRs"**

Naoki Tarui, Thuy Duong Nguyen, Masaharu Nakayama
SEEDSUPPLY INC.
Contact: naoki.tarui@seedsupply.co.jp
Website: www.seedsupply.co.jp

---

## Repository Contents

### data/

| File | Description |
|------|-------------|
| `gpcr_correlation_results.csv` | Per-target Pearson's r between Boltz-2 predicted and experimental pKd, GPCR (n=100, unified prediction batch) |
| `fpocket_gpcr_results.csv` | fpocket druggability score, pocket volume, and pocket count, GPCR (n=100) |
| `chai1_gpcr_comparison.csv` | Per-target Pearson's r for Chai-1 affinity prediction vs experimental pKd (GPCR, n=100), with Boltz-2 comparison |
| `vs_summary_all.csv` | Virtual screening performance metrics (20 GPCR targets): Pearson's r, separation index (Sep), AUROC, average precision, and enrichment factors |

### scripts/

| File | Description |
|------|-------------|
| `01_gpcr_correlation.py` | GPCR per-target correlation analysis |
| `02_fpocket_gpcr.py` | GPCR fpocket druggability analysis |
| `03_gpcr_improved_model.py` | Multi-descriptor model exploration |
| `04_residual_analysis.py` | Predicted-affinity spread analysis |
| `05_prepare_gpcr_inputs.py` | Boltz-2 input YAML preparation (requires proprietary compound data) |
| `06_run_fpocket_gpcr.py` | fpocket execution for GPCR |
| `07_vs_analysis.py` | Virtual screening enrichment analysis (separation index, AUROC, EF) |
| `08_vs_collect_results.py` | Virtual screening result collection and summary |

---

## Key Findings

### Per-target ranking accuracy (100 GPCRs, 1,270 pairs)

- Pearson's r ranged from −0.54 to +0.88 (mean 0.220, median 0.214)
- 21 targets (21%) showed strong positive correlation (r ≥ 0.5)
- 23 targets (23%) showed inverse correlation (r < 0)
- Ranking accuracy varies widely across targets, from near-perfect to inverse prediction

### Virtual screening (20 GPCR targets × 1,200 compounds)

- Screening success is governed by the separation between predicted binder and non-binder
  scores (separation index, Sep), not by rank-order correlation alone
- Sep tracks enrichment closely (Sep–EF@1% Pearson = 0.93; Sep–AUROC Pearson = 0.97),
  whereas ranking accuracy does not (r–EF@1% Pearson = 0.42; r–AUROC Pearson = 0.56)
- A high per-target r is necessary but not sufficient: several high-r targets
  (e.g., CCR2 r = 0.70, GPR85 r = 0.60, PTGER2 r = 0.58) failed to enrich true binders (EF@1% = 0)
- Conversely, some moderate-r targets screened well (e.g., ADRB1 r = 0.46, AUROC = 0.99;
  OX2R r = 0.38, AUROC = 0.92)

### Structural interpretation (post hoc)

- fpocket druggability score does not predict ranking accuracy (R² = 0.0009, p = 0.77);
  nearly all GPCRs have high druggability scores (ceiling effect)
- Binder contact consistency correlates with score separation (Pearson = 0.81 at 9 Å,
  positive-r subset), but requires experimental binder/non-binder labels and is therefore
  a retrospective interpretation, not a prospective predictor
- No label-free structural descriptor predicted screening success a priori

### Model dependence: Boltz-2 vs Chai-1 (100 GPCR targets)

- Boltz-2 and Chai-1 per-target correlations are only weakly related (Pearson's r = 0.34)
- Per-target predictability reflects model-specific behavior rather than an intrinsic
  property of the target

---

## Data Availability

The Binder2030 experimental Kd values and compound structures are proprietary to SEEDSUPPLY INC.
Aggregated and derived data supporting the manuscript are provided in this repository and the
Supplementary Information.

### Openly available here

- Boltz-2 predicted affinity values (affinity_pred_value)
- Per-target Pearson's r statistics (unified prediction batch)
- fpocket structural analysis results (druggability score, pocket volume)
- Virtual screening performance metrics (separation index, AUROC, EF)
- All Python analysis scripts

### Available upon reasonable request (subject to a data-sharing agreement)

- Experimental Kd values (pKd) measured by BST (Binder Selection Technology)
- For academic research purposes; contact naoki.tarui@seedsupply.co.jp

### Not available

- Compound structures (SMILES) and compound identifiers linked to experimental Kd data

---

## Requirements

```
Python 3.10
pandas >= 1.5
scipy >= 1.9
scikit-learn >= 1.1
numpy >= 1.23
fpocket (https://github.com/Discngine/fpocket)
Boltz-2 (https://github.com/jwohlwend/boltz)
Chai-1 (https://github.com/chaidiscovery/chai-lab)
```

---

## Citation

If you use this code or data, please cite the manuscript above. A link to the published
version will be added upon acceptance.
