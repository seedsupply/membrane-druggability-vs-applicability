"""
06_run_fpocket_gpcr.py
======================
fpocket execution and result parsing for GPCR targets
Downloads AlphaFold2 structures and runs fpocket analysis

Requirements: fpocket (https://github.com/Discngine/fpocket)

Input:  gpcr_correlation_results.csv (UniProt IDs)
Output: fpocket_gpcr_results.csv
"""

import os
import subprocess
import pandas as pd
import numpy as np
import glob
import re
import requests
import warnings
warnings.filterwarnings('ignore')

def download_alphafold_structure(uniprot_id, output_dir):
    """Download AlphaFold2 predicted structure from EBI."""
    url = (f"https://alphafold.ebi.ac.uk/files/"
           f"AF-{uniprot_id}-F1-model_v4.pdb")
    pdb_path = os.path.join(output_dir, f"{uniprot_id}.pdb")

    if os.path.exists(pdb_path):
        return pdb_path

    response = requests.get(url, timeout=30)
    if response.status_code == 200:
        with open(pdb_path, 'w') as f:
            f.write(response.text)
        return pdb_path
    return None

def run_fpocket(pdb_file, output_dir):
    """Run fpocket on a single PDB file."""
    cmd = ['fpocket', '-f', pdb_file]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=output_dir)
    return result.returncode == 0

def parse_fpocket_info(fpocket_dir, target):
    """Parse fpocket output directory for best pocket statistics."""
    # fpocket creates {target}_out directory
    out_dir = os.path.join(fpocket_dir, f"{target}_out")
    if not os.path.exists(out_dir):
        return None, None, 0

    pocket_files = glob.glob(os.path.join(out_dir, 'pockets', 'pocket*_atm.pdb'))
    info_file = os.path.join(out_dir, f"{target}_info.txt")

    if not os.path.exists(info_file):
        return None, None, len(pocket_files)

    drug_scores = []
    pocket_volumes = []

    with open(info_file) as f:
        content = f.read()

    pockets = re.split(r'Pocket \d+', content)[1:]
    for pocket in pockets:
        score = re.search(r'Druggability Score\s*:\s*([\d.]+)', pocket)
        volume = re.search(r'Volume\s*:\s*([\d.]+)', pocket)
        if score:
            drug_scores.append(float(score.group(1)))
        if volume:
            pocket_volumes.append(float(volume.group(1)))

    if not drug_scores:
        return None, None, 0

    best_idx = int(np.argmax(drug_scores))
    return (round(max(drug_scores), 3),
            round(pocket_volumes[best_idx], 3) if pocket_volumes else None,
            len(drug_scores))

def main():
    # Load target list
    targets_df = pd.read_csv('data/gpcr_correlation_results.csv')[['Target','UniProt']]

    struct_dir = 'structures/gpcr'
    fpocket_dir = 'fpocket_results/gpcr'
    os.makedirs(struct_dir, exist_ok=True)
    os.makedirs(fpocket_dir, exist_ok=True)

    results = []
    for _, row in targets_df.iterrows():
        target = row['Target']
        uniprot = row['UniProt']

        print(f"Processing {target} ({uniprot})...")

        # Download structure
        pdb_path = download_alphafold_structure(uniprot, struct_dir)
        if not pdb_path:
            print(f"  Failed to download structure for {uniprot}")
            continue

        # Copy PDB with target name for fpocket
        import shutil
        target_pdb = os.path.join(fpocket_dir, f"{target}.pdb")
        shutil.copy(pdb_path, target_pdb)

        # Run fpocket
        success = run_fpocket(target_pdb, fpocket_dir)
        if not success:
            print(f"  fpocket failed for {target}")
            continue

        # Parse results
        drug_score, pocket_volume, n_pockets = parse_fpocket_info(
            fpocket_dir, target)

        results.append({
            'UniProt': uniprot,
            'Target': target,
            'drug_score': drug_score,
            'pocket_volume': pocket_volume,
            'n_pockets': n_pockets
        })
        print(f"  drug_score={drug_score}, pocket_volume={pocket_volume}, "
              f"n_pockets={n_pockets}")

    df = pd.DataFrame(results)
    df.to_csv('data/fpocket_gpcr_results.csv', index=False)
    print(f"\nSaved: data/fpocket_gpcr_results.csv ({len(df)} targets)")
    print(f"Median drug_score: {df['drug_score'].median():.3f}")
    print(f"Median pocket_volume: {df['pocket_volume'].median():.1f}")

if __name__ == '__main__':
    main()
