import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

GREEN = "\033[1;32m"
RESET = "\033[0m"

#default 8 clusters max
def infer_populations_automatically(vcf_path, output_dir, max_clusters=8):
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
    
    #Calculate optimal number of clusters using Elbow plot
    
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

    print("-> Determining optimal number of clusters...")
    inertia = []
    silhouettes = []
    K = range(2, max_clusters + 1)

    for k in K:
        km = KMeans(n_clusters=k, random_state=42)
        labels = km.fit_predict(coords)
        inertia.append(km.inertia_)
        silhouettes.append(silhouette_score(coords, labels))
    
    #plot Elbow
    
    plt.figure()
    plt.plot(K, inertia, marker='o')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Total inertia')
    plt.title('Elbow method')
    plt.grid(True)
    elbow_path = os.path.join(output_dir, 'elbow_clusters.png')
    plt.savefig(elbow_path)
    print(f"Elbow plot saved to {elbow_path}")

    best_k = K[np.argmax(silhouettes)]
    print(f"Optimal number of clusters (based on silhouette): {best_k}")

    kmeans = KMeans(n_clusters=best_k, random_state=42)
    cluster_labels = kmeans.fit_predict(coords)
    
    #generate a file with all the clusters specified
    
    df = pd.DataFrame({
        'Sample': samples,
        'PC1': coords[:, 0],
        'PC2': coords[:, 1],
        'Cluster': cluster_labels
    })

    output_csv = os.path.join(output_dir, 'auto_populations.csv')
    df.to_csv(output_csv, index=False)
    print(f"{GREEN} ✅ Inferred populations saved to {output_csv} {RESET}")
    
    #plot PCA
    
    plt.figure(figsize=(10, 6))
    for cluster in sorted(df['Cluster'].unique()):
        sub = df[df['Cluster'] == cluster]
        plt.scatter(sub['PC1'], sub['PC2'], label=f'Cluster {cluster}', alpha=0.7)
    plt.title(f'PCA + Automatic Clustering (k={best_k})')
    plt.xlabel(f'PC1 ({explained[0]:.1f}%)')
    plt.ylabel(f'PC2 ({explained[1]:.1f}%)')
    plt.legend()
    plt.tight_layout()
    pca_path = os.path.join(output_dir, 'pca_clusters.png')
    plt.savefig(pca_path)
    print(f"PCA plot saved to {pca_path}")

    return df
