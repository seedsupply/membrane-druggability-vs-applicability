"""
04_residual_analysis.py
=======================
Residual variance analysis for GPCR Boltz-2 predictions
Analyzes std_aff_value as post-hoc prediction reliability indicator

Input:  GPCR_merged_results2.csv (proprietary)
        gpcr_correlation_results.csv
Output: Console output with statistics
"""

import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def calculate_std_aff_value(df):
    """Calculate standard deviation of predicted affinity values per target."""
    return df.groupby('Target')['affinity_pred_value_neg'].std().reset_index(
        ).rename(columns={'affinity_pred_value_neg': 'std_aff_value'})

def main():
    # Load data
    try:
        df = pd.read_csv('data/GPCR_merged_results2.csv')
    except FileNotFoundError:
        print("Error: GPCR_merged_results2.csv not found.")
        print("This file requires a Data Use Agreement with SEEDSUPPLY INC.")
        return

    corr = pd.read_csv('data/gpcr_correlation_results.csv')

    # Calculate std_aff_value
    std_aff = calculate_std_aff_value(df)
    merged = corr.merge(std_aff, on='Target', how='left')

    print("=" * 60)
    print("std_aff_value Analysis (GPCR, n=100)")
    print("=" * 60)

    # Correlation with Pearson's r
    valid = merged.dropna(subset=['std_aff_value', 'Pearson_r'])
    r, p = stats.pearsonr(valid['std_aff_value'], valid['Pearson_r'])
    print(f"Correlation with Pearson's r: r={r:.3f}, P={p:.4f}")
    print(f"Variance explained: r²={r**2*100:.1f}%")

    # Threshold analysis
    print("\n--- Prediction reliability by std_aff_value threshold ---")
    for threshold in [0.5, 0.75, 1.0, 1.5]:
        high = merged[merged['std_aff_value'] >= threshold]['Pearson_r']
        low = merged[merged['std_aff_value'] < threshold]['Pearson_r']
        if len(high) > 0 and len(low) > 0:
            stat, pval = stats.mannwhitneyu(high, low, alternative='greater')
            print(f"threshold={threshold}: high(n={len(high)}) median={high.median():.3f} "
                  f"vs low(n={len(low)}) median={low.median():.3f}, P={pval:.4f}")

    # CWxP stratified analysis
    if 'has_CWxP' in merged.columns:
        print("\n--- std_aff_value by CWxP group ---")
        for cwxp in [1, 0]:
            group = merged[merged['has_CWxP'] == cwxp]['std_aff_value'].dropna()
            label = 'CWxP(+)' if cwxp == 1 else 'CWxP(-)'
            print(f"{label}: median={group.median():.3f}, mean={group.mean():.3f}")

    print(f"\nNote: std_aff_value is a post-hoc indicator —")
    print(f"computable only after Boltz-2 prediction execution.")

if __name__ == '__main__':
    main()
