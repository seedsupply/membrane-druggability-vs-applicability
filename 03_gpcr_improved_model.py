"""
03_gpcr_improved_model.py
=========================
GPCR improved prediction model
Variance explained by individual features and combined model

Input:  gpcr_correlation_results.csv
        fpocket_gpcr_results.csv
        GPCR_merged_results2.csv (proprietary, for std_aff_value)
Output: Supplementary Table 2 statistics
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings('ignore')

def calculate_r2(x, y):
    """Calculate r² and P value for single feature."""
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 5:
        return np.nan, np.nan
    r, p = stats.pearsonr(x[mask], y[mask])
    return round(r**2 * 100, 1), round(p, 4)

def combined_model(X, y, feature_names):
    """Multiple linear regression with StandardScaler."""
    mask = ~np.isnan(X).any(axis=1) & ~np.isnan(y)
    X_clean = X[mask]
    y_clean = y[mask]

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_clean)

    model = LinearRegression()
    model.fit(X_scaled, y_clean)
    y_pred = model.predict(X_scaled)
    R2 = r2_score(y_clean, y_pred)

    return round(R2 * 100, 1)

def main():
    # Load data
    corr = pd.read_csv('data/gpcr_correlation_results.csv')
    fpocket = pd.read_csv('data/fpocket_gpcr_results.csv')

    df = corr.merge(fpocket[['UniProt','druggability_score','pocket_volume','n_pockets']],
                    on='UniProt', how='left')

    y = df['Pearson_r'].values

    print("=" * 60)
    print("GPCR Feature Analysis (n=100 targets)")
    print("=" * 60)

    # 2. ipTM score (if available)
    if 'iptm_mean' in df.columns:
        r2, p = calculate_r2(df['iptm_mean'].values, y)
        print(f"ipTM score (mean):        r²={r2}%, P={p}")

    # 3. Compound logP
    if 'logP_mean' in df.columns:
        r2, p = calculate_r2(df['logP_mean'].values, y)
        print(f"Compound logP (mean):     r²={r2}%, P={p}")

    # 4. Compound MW diversity
    if 'MW_diversity' in df.columns:
        r2, p = calculate_r2(df['MW_diversity'].values, y)
        print(f"Compound MW diversity:    r²={r2}%, P={p}")

    # 5. fpocket druggability_score
    r2, p = calculate_r2(df['druggability_score'].values, y)
    print(f"fpocket druggability_score:       r²={r2}%, P={p}")

    # 6. std_aff_value (post-hoc)
    if 'std_aff_value' in df.columns:
        r2, p = calculate_r2(df['std_aff_value'].values, y)
        print(f"std_aff_value (post-hoc): r²={r2}%, P={p}")

    print("\n" + "=" * 60)
    print("Combined Model (GPCR)")
    print("=" * 60)

    # GPCR combined model
    feature_cols = []
    for col in ['logP_mean', 'iptm_mean', 'std_aff_value']:
        if col in df.columns:
            feature_cols.append(col)

    if len(feature_cols) >= 2:
        X = df[feature_cols].values
        R2 = combined_model(X, y, feature_cols)
        print(f"Features: {feature_cols}")
        print(f"Combined R² = {R2}%")

if __name__ == '__main__':
    main()
