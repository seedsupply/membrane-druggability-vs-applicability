"""
08_vs_collect_results.py
========================
Collect and summarize Boltz-2 virtual screening results
Generates figures and statistical comparisons

Input:  data/vs_summary_all.csv
        results/{TARGET}_vs_results.csv
Output: Summary statistics and plots
"""

import os
import pandas as pd
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

def load_summary():
    """Load virtual screening summary data."""
    try:
        df = pd.read_csv('data/vs_summary_all.csv')
        return df
    except FileNotFoundError:
        print("Error: Run 07_vs_analysis.py first.")
        return None

def group_analysis(df):
    print("=" * 60)
    print(f"Virtual Screening Performance Summary (n={len(df)} targets)")
    print("=" * 60)

    # By Pearson r group
    df['r_group'] = pd.cut(
        df['Pearson_r'],
        bins=[-1, 0, 0.5, 1],
        labels=['Low (r<0)', 'Mid (0≤r<0.5)', 'High (r≥0.5)']
    )

    print("\n--- By Pearson r group ---")
    for group in ['High (r≥0.5)', 'Mid (0≤r<0.5)', 'Low (r<0)']:
        subset = df[df['r_group'] == group]
        if len(subset) == 0:
            continue
        print(f"\n{group} (n={len(subset)}):")
        print(f"  AUROC: median={subset['AUROC'].median():.3f}, "
              f"range={subset['AUROC'].min():.3f}–{subset['AUROC'].max():.3f}")
        print(f"  EF@1%: median={subset['EF_1pct'].median():.1f}, "
              f"range={subset['EF_1pct'].min():.1f}–{subset['EF_1pct'].max():.1f}")

    # Statistical tests: ranking accuracy vs screening performance
    print("\n--- Ranking accuracy vs screening performance ---")
    r_auroc, p_auroc = stats.pearsonr(df['Pearson_r'], df['AUROC'])
    r_ef1, p_ef1 = stats.pearsonr(df['Pearson_r'], df['EF_1pct'])
    print(f"\nPearson r vs AUROC: r={r_auroc:.3f}, P={p_auroc:.4f}")
    print(f"Pearson r vs EF@1%: r={r_ef1:.3f}, P={p_ef1:.4f}")

def print_per_target_table(df):
    """Print per-target results table."""
    print("\n" + "=" * 80)
    print("Per-target Results (ordered by Pearson's r)")
    print("=" * 80)
    df_sorted = df.sort_values('Pearson_r', ascending=False)
    cols = ['Target', 'Pearson_r', 'Sep', 'AUROC', 'Average_Precision',
            'EF_0.5pct', 'EF_1pct', 'EF_5pct', 'EF_10pct']
    available_cols = [c for c in cols if c in df_sorted.columns]
    print(df_sorted[available_cols].to_string(index=False))

def main():
    df = load_summary()
    if df is None:
        return

    # Rename EF columns if needed
    ef_rename = {}
    for col in df.columns:
        if 'EF_' in col and 'pct' in col:
            new_col = col.replace('pct', '').replace('_', '@').replace('EF@', 'EF_') + 'pct'
            ef_rename[col] = col  # keep as is

    print_per_target_table(df)
    group_analysis(df)

    print("\n" + "=" * 60)
    print("Key Finding:")
    high = df[df['Pearson_r'] >= 0.5]
    low = df[df['Pearson_r'] < 0]
    print(f"High-accuracy targets (r≥0.5, n={len(high)}): "
          f"median AUROC={high['AUROC'].median():.3f}, "
          f"median EF@1%={high['EF_1pct'].median():.1f}×")
    print(f"Low-accuracy targets (r<0, n={len(low)}): "
          f"median AUROC={low['AUROC'].median():.3f}, "
          f"median EF@1%={low['EF_1pct'].median():.1f}×")

if __name__ == '__main__':
    main()
