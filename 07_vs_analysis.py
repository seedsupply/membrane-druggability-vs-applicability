"""
07_vs_analysis.py
=================
Virtual screening enrichment analysis
Calculates separation index (Sep), AUROC, and Enrichment Factor for GPCR targets

Input:  boltz2_outputs/{TARGET}/boltz_results_{TARGET}/predictions/
        data/GPCR_merged_results2.csv (proprietary, for binder labels)
Output: data/vs_summary_all.csv (per-target Sep, AUROC, AP, and EF@0.5/1/2/5/10%)
        results/{TARGET}_vs_results.csv
        The distributed data/vs_summary_all.csv reports the key columns used in
        the manuscript (Target, n_binders, Pearson_r, Sep, AUROC,
        Average_Precision, EF@0.5/1/5/10%); this script writes the full set of
        computed metrics, from which those columns are drawn.

Usage:
    python 07_vs_analysis.py --target CB1R
    python 07_vs_analysis.py --target all
"""

import os
import json
import argparse
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
import warnings
warnings.filterwarnings('ignore')

# Target metadata
TARGET_INFO = {
    'ADRB1': {'pearson_r': 0.464},
    'OX2R': {'pearson_r': 0.382},
    'CB1R': {'pearson_r': 0.830},
    'MTNR1B': {'pearson_r': 0.876},
    'DRD3': {'pearson_r': 0.571},
    'CXCR4': {'pearson_r': 0.449},
    'HTR6': {'pearson_r': -0.142},
    'HTR4': {'pearson_r': -0.116},
    'CCR3': {'pearson_r': -0.542},
    'CNR2': {'pearson_r': -0.539},
    'GPR37L1': {'pearson_r': -0.032},
    'SMO': {'pearson_r': -0.308},
    'OX1R': {'pearson_r': 0.659},
    'AGTR1': {'pearson_r': 0.570},
    'CHRM3': {'pearson_r': 0.814},
    'F2R': {'pearson_r': 0.537},
    'HTR7': {'pearson_r': 0.673},
    'PTGER2': {'pearson_r': 0.581},
    'CCR2': {'pearson_r': 0.699},
    'GPR85': {'pearson_r': 0.602},
}

EF_PERCENTILES = [0.5, 1, 2, 5, 10]

def collect_boltz2_results(target, output_base='boltz2_outputs'):
    """Collect Boltz-2 affinity predictions from output directory."""
    results = []
    pred_dir = os.path.join(output_base, target,
                            f'boltz_results_{target}', 'predictions')

    if not os.path.exists(pred_dir):
        print(f"  Warning: {pred_dir} not found")
        return pd.DataFrame()

    for compound_dir in os.listdir(pred_dir):
        json_file = os.path.join(pred_dir, compound_dir,
                                 f'affinity_{compound_dir}.json')
        if not os.path.exists(json_file):
            continue
        try:
            with open(json_file) as f:
                data = json.load(f)
            pred_value = data.get('affinity_pred_value')
            if pred_value is not None:
                results.append({
                    'Compound': compound_dir,
                    'pred_pKd': -float(pred_value)
                })
        except Exception:
            continue

    return pd.DataFrame(results)

def calculate_enrichment(df, target):
    """Calculate AUROC and EF at multiple percentiles."""
    df = df.sort_values('pred_pKd', ascending=False).reset_index(drop=True)
    n_total = len(df)
    n_binders = df['Is_Binder'].sum()

    if n_binders == 0:
        return None

    binder_rate = n_binders / n_total

    # AUROC and AP
    try:
        auroc = roc_auc_score(df['Is_Binder'], df['pred_pKd'])
        ap = average_precision_score(df['Is_Binder'], df['pred_pKd'])
    except Exception:
        auroc = ap = np.nan

    # Separation index (Sep):
    #   Sep = (mean_binder - mean_nonbinder) / (std_binder + std_nonbinder)
    # computed here in pred_pKd space (higher pred_pKd = stronger predicted
    # binding), so that a positive Sep indicates binders receive
    # systematically stronger predicted affinities than non-binders.
    binder_scores = df.loc[df['Is_Binder'] == 1, 'pred_pKd']
    nonbinder_scores = df.loc[df['Is_Binder'] == 0, 'pred_pKd']
    std_sum = binder_scores.std() + nonbinder_scores.std()
    if len(binder_scores) > 1 and len(nonbinder_scores) > 1 and std_sum > 0:
        sep = round((binder_scores.mean() - nonbinder_scores.mean()) / std_sum, 3)
    else:
        sep = np.nan

    # Enrichment Factors
    ef_dict = {}
    for pct in EF_PERCENTILES:
        n_top = max(1, int(n_total * pct / 100))
        hits = df.head(n_top)['Is_Binder'].sum()
        ef = (hits / n_top) / binder_rate if binder_rate > 0 else np.nan
        ef_dict[f'EF_{pct}pct'] = round(ef, 2)

    result = {
        'Target': target,
        'n_total': n_total,
        'n_binders': n_binders,
        'binder_rate_pct': round(binder_rate * 100, 2),
        'Sep': sep,
        'AUROC': round(auroc, 4),
        'Average_Precision': round(ap, 4),
        **ef_dict,
        'Pearson_r': TARGET_INFO.get(target, {}).get('pearson_r', np.nan),
    }
    return result, df

def main(target_arg):
    # Load binder data
    try:
        binders = pd.read_csv('data/GPCR_merged_results2.csv')
    except FileNotFoundError:
        print("Error: GPCR_merged_results2.csv not found.")
        print("This file requires a Data Use Agreement with SEEDSUPPLY INC.")
        return

    os.makedirs('results', exist_ok=True)

    targets = list(TARGET_INFO.keys()) if target_arg == 'all' else [target_arg]
    summary = []

    for target in targets:
        print(f"\n{'='*50}")
        print(f"Analyzing: {target}")

        # Collect predictions
        pred_df = collect_boltz2_results(target)
        if pred_df.empty:
            print(f"  No results found. Run Boltz-2 first.")
            continue

        # Merge with binder labels
        target_binders = binders[binders['Target'] == target][
            ['Compound', 'pKd']].copy()
        target_binders['Is_Binder'] = 1

        merged = pred_df.merge(target_binders, on='Compound', how='left')
        merged['Is_Binder'] = merged['Is_Binder'].fillna(0).astype(int)

        n_binders = merged['Is_Binder'].sum()
        print(f"  Compounds: {len(merged)}, Binders: {n_binders}")

        # Calculate enrichment
        result = calculate_enrichment(merged, target)
        if result is None:
            continue
        metrics, sorted_df = result

        summary.append(metrics)
        sorted_df.to_csv(f'results/{target}_vs_results.csv', index=False)

        print(f"  Sep: {metrics['Sep']}")
        print(f"  AUROC: {metrics['AUROC']:.4f}")
        for pct in EF_PERCENTILES:
            print(f"  EF@{pct}%: {metrics[f'EF_{pct}pct']:.2f}")

    if summary:
        summary_df = pd.DataFrame(summary)
        summary_df = summary_df.sort_values('Sep', ascending=False)
        summary_df.to_csv('data/vs_summary_all.csv', index=False)
        print(f"\nSaved: data/vs_summary_all.csv")
        print(summary_df[['Target','Pearson_r','Sep','AUROC',
                          'EF_1pct','EF_5pct']].to_string(index=False))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Virtual screening enrichment analysis'
    )
    parser.add_argument('--target', default='all',
                        help='Target name or "all"')
    args = parser.parse_args()
    main(args.target)
