# Binder2030 × Boltz-2 Benchmark: GPCR Virtual Screening

This repository contains analysis code and Boltz-2 prediction results for the following manuscript:

**"Structural determinants of AI-based affinity prediction define virtual screening success in GPCRs"**

Naoki Tarui, Thuy Duong Nguyen, Masaharu Nakayama  
SEEDSUPPLY INC.  
Contact: naoki.tarui@seedsupply.co.jp  
Website: www.seedsupply.co.jp

---

## Repository Contents

### data/

| File | Description |
|------|-------------|
| `gpcr_correlation_results.csv` | Per-target Pearson's r, GPCR (n=100) |
| `fpocket_gpcr_results.csv` | fpocket drug_score, GPCR (n=100) |
| `chai1_gpcr_comparison.csv` | Per-target Pearson's r for Chai-1 ipTM-based affinity prediction vs experimental pKd (GPCR, n=100); includes Boltz-2 comparison data |
| `vs_summary_all.csv` | Virtual screening performance metrics (11 GPCR targets) |

### scripts/

| File | Description |
|------|-------------|
| `01_gpcr_correlation.py` | GPCR per-target correlation analysis |
| `02_fpocket_gpcr.py` | GPCR fpocket analysis |
| `03_gpcr_improved_model.py` | GPCR improved prediction model |
| `04_residual_analysis.py` | Residual variance analysis |
| `05_prepare_gpcr_inputs.py` | Boltz-2 input YAML preparation (requires proprietary data) |
| `06_run_fpocket_gpcr.py` | fpocket execution for GPCR |
| `07_vs_analysis.py` | Virtual screening enrichment analysis |
| `08_vs_collect_results.py` | Virtual screening result collection and AUROC/EF calculation |

---

## Key Findings

### GPCR (100 targets, 1,270 pairs)

- Mean Pearson's r = 0.214
- CWxP motif in TM6 as structural predictor (P = 0.036*)
- CWxP subtype refinement improves prediction (P = 0.010**)
- Best combined model R² = 21.4%

### Virtual Screening Validation (11 GPCR targets × 1,200 compounds)

- High-accuracy targets (r ≥ 0.5): AUROC = 0.68–0.99, EF@1% = 8–58×
- Low-accuracy targets (r < 0): AUROC ≈ 0.5, EF@1% = 0
- CWxP(+) status predicts virtual screening success

### Exploratory Comparison: Boltz-2 vs Chai-1 (100 GPCR targets)

- Boltz-2 and Chai-1 show only weak concordance (Pearson's r = 0.317, P = 0.0013)
- CWxP(+) GPCRs: Boltz-2 mean r = 0.285 vs Chai-1 mean r = 0.106
- CWxP-associated structural principle appears model-specific to Boltz-2

---

## Data Availability

The Binder2030 experimental Kd values are proprietary to SEEDSUPPLY INC. and are not included in this repository.

### What is available here

- Boltz-2 predicted affinity values (affinity_pred_value)
- Per-target Pearson's r statistics
- fpocket structural analysis results (drug_score, pocket_volume)
- Virtual screening performance metrics (AUROC, EF)
- All Python analysis scripts

### What requires a Data Use Agreement

- Experimental Kd values (pKd) measured by BST (Binder Selection Technology)
- Available for non-commercial research upon request
- Contact: naoki.tarui@seedsupply.co.jp

### Not available

- Compound structures (SMILES)
- Compound identifiers linked to experimental Kd data

---

## Requirements

```
Python 3.10
pandas >= 1.5
scipy >= 1.9
scikit-learn >= 1.1
rdkit >= 2022.09
numpy >= 1.23
fpocket (https://github.com/Discngine/fpocket)
```

---

## Citation

If you use this code or data, please cite:

Tarui N, Nguyen TD, Nakayama M. Structural determinants of AI-based affinity prediction define virtual screening success in GPCRs. *Communications AI & Computing*, 2026. [DOI to be added upon publication]

---

## License

Code: MIT License  
Data: See Data Availability section above.

© 2026 SEEDSUPPLY INC. All rights reserved.
