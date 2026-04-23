"""
07_vs_analysis.py
=================
Virtual screening enrichment analysis
Calculates AUROC and Enrichment Factor for GPCR targets

Input:  boltz2_outputs/{TARGET}/boltz_results_{TARGET}/predictions/
        data/GPCR_merged_results2.csv (proprietary, for binder labels)
Output: data/vs_summary_all.csv
        results/{TARGET}_vs_results.csv

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
    'CB1R':    {'pearson_r': 0.838, 'has_CWxP': 1},
    'MTNR1B':  {'pearson_r': 0.836, 'has_CWxP': 1},
    'ADRB1':   {'pearson_r': 0.592, 'has_CWxP': 1},
    'CXCR4':   {'pearson_r': 0.451, 'has_CWxP': 1},
    'OX2R':    {'pearson_r': 0.326, 'has_CWxP': 0},
    'HTR6':    {'pearson_r': 0.009, 'has_CWxP': 0},
    'HTR4':    {'pearson_r': -0.040, 'has_CWxP': 1},
    'GPR37L1': {'pearson_r': -0.303, 'has_CWxP': 0},
    'CNR2':    {'pearson_r': -0.423, 'has_CWxP': 1},
    'SMO':     {'pearson_r': -0.557, 'has_CWxP': 0},
    'CCR3':    {'pearson_r': -0.614, 'has_CWxP': 0},
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
        'AUROC': round(auroc, 4),
        'Average_Precision': round(ap, 4),
        **ef_dict,
        'Pearson_r': TARGET_INFO.get(target, {}).get('pearson_r', np.nan),
        'has_CWxP': TARGET_INFO.get(target, {}).get('has_CWxP', np.nan),
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

        print(f"  AUROC: {metrics['AUROC']:.4f}")
        for pct in EF_PERCENTILES:
            print(f"  EF@{pct}%: {metrics[f'EF_{pct}pct']:.2f}")

    if summary:
        summary_df = pd.DataFrame(summary)
        summary_df.to_csv('data/vs_summary_all.csv', index=False)
        print(f"\nSaved: data/vs_summary_all.csv")
        print(summary_df[['Target','Pearson_r','AUROC',
                          'EF_1pct','EF_5pct']].to_string(index=False))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Virtual screening enrichment analysis'
    )
    parser.add_argument('--target', default='all',
                        help='Target name or "all"')
    args = parser.parse_args()
    main(args.target)
