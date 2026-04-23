"""
01_gpcr_correlation.py
======================
GPCR per-target Pearson's r analysis
Boltz-2 predicted vs experimental pKd

Input:  GPCR_merged_results2.csv (proprietary)
Output: gpcr_correlation_results.csv
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def calculate_per_target_correlation(df):
    """Calculate per-target Pearson's r between Boltz-2 predicted and experimental pKd."""
    results = []
    for target, group in df.groupby('Target'):
        if len(group) < 3:
            continue
        r, p = stats.pearsonr(group['affinity_pred_value_neg'], group['pKd'])
        results.append({
            'Target': target,
            'UniProt': group['UniProt'].iloc[0],
            'n': len(group),
            'Pearson_r': round(r, 3),
            'p_value': round(p, 4),
            'GPCR_class': group['GPCR_class'].iloc[0] if 'GPCR_class' in group.columns else 'Class A',
            'has_DRY': int(group['has_DRY'].iloc[0]) if 'has_DRY' in group.columns else None,
            'has_NPxxY': int(group['has_NPxxY'].iloc[0]) if 'has_NPxxY' in group.columns else None,
            'has_CWxP': int(group['has_CWxP'].iloc[0]) if 'has_CWxP' in group.columns else None,
            'cys_count': group['cys_count'].iloc[0] if 'cys_count' in group.columns else None,
            'trp_count': group['trp_count'].iloc[0] if 'trp_count' in group.columns else None,
            'asp_count': group['asp_count'].iloc[0] if 'asp_count' in group.columns else None,
        })
    return pd.DataFrame(results).sort_values('Pearson_r', ascending=False)

def main():
    # Load data (requires proprietary Binder2030 data)
    try:
        df = pd.read_csv('data/GPCR_merged_results2.csv')
    except FileNotFoundError:
        print("Error: GPCR_merged_results2.csv not found.")
        print("This file requires a Data Use Agreement with SEEDSUPPLY INC.")
        print("Contact: naoki.tarui@seedsupply.co.jp")
        return

    print(f"Loaded {len(df)} compound-target pairs across {df['Target'].nunique()} GPCRs")

    results = calculate_per_target_correlation(df)

    # Summary statistics
    print(f"\nSummary (n={len(results)} targets):")
    print(f"  Mean Pearson's r: {results['Pearson_r'].mean():.3f}")
    print(f"  Median Pearson's r: {results['Pearson_r'].median():.3f}")
    print(f"  Strong predictors (r >= 0.5): {(results['Pearson_r'] >= 0.5).sum()}")
    print(f"  Inverse predictors (r < 0): {(results['Pearson_r'] < 0).sum()}")

    # CWxP analysis
    if 'has_CWxP' in results.columns:
        cwxp_pos = results[results['has_CWxP'] == 1]['Pearson_r']
        cwxp_neg = results[results['has_CWxP'] == 0]['Pearson_r']
        stat, pval = stats.mannwhitneyu(cwxp_pos, cwxp_neg, alternative='greater')
        print(f"\nCWxP(+) vs CWxP(-): median {cwxp_pos.median():.3f} vs {cwxp_neg.median():.3f}, P={pval:.4f}")

    results.to_csv('data/gpcr_correlation_results.csv', index=False)
    print(f"\nSaved: data/gpcr_correlation_results.csv")

if __name__ == '__main__':
    main()
