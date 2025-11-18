import allel
import numpy as np
import pandas as pd
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def pca_and_save_coords(vcf_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    print("\n[PHASE 1] Loading VCF...")
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    samples = callset['samples']

    print("-> Genotypes loaded:", genotypes.shape)

    print("-> Filtering out variants with missing data...")
    gn_alt = genotypes.to_n_alt()
    gn_alt = gn_alt[~np.any(np.isnan(gn_alt), axis=1)]
    print("-> Shape after filtering:", gn_alt.shape)

    print("-> Performing PCA...")
    gn_transposed = gn_alt.T
    scaler = StandardScaler()
    gn_scaled = scaler.fit_transform(gn_transposed)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(gn_scaled)
    explained = pca.explained_variance_ratio_ * 100
    print(f"-> Variance explained: {explained[:2]}")

    df = pd.DataFrame({
        'Sample': samples,
        'PC1': coords[:, 0],
        'PC2': coords[:, 1]
    })
    coords_csv = os.path.join(output_dir, 'pca_coords.csv')
    df.to_csv(coords_csv, index=False)
    print(f"PCA coordinates saved to {coords_csv}")

    # Guarda también la varianza explicada para usar después
    var_exp_path = os.path.join(output_dir, 'explained_variance.npy')
    np.save(var_exp_path, explained, allow_pickle=True)
    print(f"Explained variance saved to {var_exp_path}")

    return coords_csv, var_exp_path
