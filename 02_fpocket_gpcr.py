"""
02_fpocket_gpcr.py
==================
fpocket druggability analysis for GPCR targets
Calculates drug_score and pocket_volume for GPCR structures

Requirements: fpocket (https://github.com/Discngine/fpocket)

Input:  AlphaFold2 PDB structures for GPCR targets
Output: fpocket_gpcr_results.csv
"""

import os
import subprocess
import pandas as pd
import numpy as np
import glob
import re
import warnings
warnings.filterwarnings('ignore')

def run_fpocket(pdb_file, output_dir):
    """Run fpocket on a single PDB file."""
    cmd = f"fpocket -f {pdb_file} -o {output_dir}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.returncode == 0

def parse_fpocket_results(fpocket_output_dir):
    """Parse fpocket output to extract drug_score and pocket_volume."""
    info_file = os.path.join(fpocket_output_dir, '*_info.txt')
    info_files = glob.glob(info_file)

    if not info_files:
        return None, None, 0

    drug_scores = []
    pocket_volumes = []

    for f in info_files:
        with open(f) as fh:
            content = fh.read()
        score = re.search(r'Druggability Score\s*:\s*([\d.]+)', content)
        volume = re.search(r'Volume\s*:\s*([\d.]+)', content)
        if score:
            drug_scores.append(float(score.group(1)))
        if volume:
            pocket_volumes.append(float(volume.group(1)))

    if not drug_scores:
        return None, None, 0

    # Return best pocket (highest drug_score)
    best_idx = drug_scores.index(max(drug_scores))
    return (round(max(drug_scores), 3),
            round(pocket_volumes[best_idx], 3) if pocket_volumes else None,
            len(drug_scores))

def main():
    # Load GPCR target list
    try:
        targets = pd.read_csv('data/gpcr_correlation_results.csv')[['Target', 'UniProt']]
    except FileNotFoundError:
        print("Error: Run 01_gpcr_correlation.py first.")
        return

    pdb_dir = 'structures/gpcr'
    output_dir = 'fpocket_results/gpcr'
    os.makedirs(output_dir, exist_ok=True)

    results = []
    for _, row in targets.iterrows():
        target = row['Target']
        uniprot = row['UniProt']
        pdb_file = os.path.join(pdb_dir, f"{uniprot}.pdb")

        if not os.path.exists(pdb_file):
            print(f"  Skipping {target}: PDB not found")
            continue

        target_out = os.path.join(output_dir, target)
        os.makedirs(target_out, exist_ok=True)

        print(f"  Running fpocket: {target}...")
        success = run_fpocket(pdb_file, target_out)

        if success:
            drug_score, pocket_volume, n_pockets = parse_fpocket_results(target_out)
            results.append({
                'UniProt': uniprot,
                'Target': target,
                'drug_score': drug_score,
                'pocket_volume': pocket_volume,
                'n_pockets': n_pockets
            })
            print(f"    drug_score={drug_score}, n_pockets={n_pockets}")

    df = pd.DataFrame(results)
    df.to_csv('data/fpocket_gpcr_results.csv', index=False)
    print(f"\nSaved: data/fpocket_gpcr_results.csv ({len(df)} targets)")
    print(f"Median drug_score: {df['drug_score'].median():.3f}")

if __name__ == '__main__':
    main()
