import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from k_means_constrained import KMeansConstrained
from kneed import KneeLocator

GREEN = "\033[1;32m"
RESET = "\033[0m"

# Función para calcular distancias euclídeas manualmente
def euclidean_distances(X):
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            dist = np.linalg.norm(X[i] - X[j])
            D[i, j] = dist
            D[j, i] = dist
    return D

def simple_silhouette(X, labels):
    D = euclidean_distances(X)
    n = X.shape[0]
    S = np.zeros(n)
    for i in range(n):
        same_cluster = labels == labels[i]
        other_clusters = labels != labels[i]

        if np.sum(same_cluster) > 1:
            a = np.mean(D[i, same_cluster][D[i, same_cluster] != 0])
        else:
            a = 0

        b = np.min([
            np.mean(D[i, labels == lbl])
            for lbl in set(labels) if lbl != labels[i]
        ])
        S[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0
    return np.mean(S)

def cluster_and_plot(coords_csv, var_exp_path, output_dir, max_clusters=8, min_cluster_size=5):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(coords_csv)
    coords = df[['PC1', 'PC2']].values
    samples = df['Sample']

    explained = np.load(var_exp_path, allow_pickle=True)

    print("-> Determining optimal number of clusters...")
    inertia = []
    K = range(2, max_clusters + 1)

    for k in K:
        km = KMeansConstrained(n_clusters=k, size_min=min_cluster_size, random_state=42)
        km.fit(coords)
        inertia.append(km.inertia_)

    plt.figure()
    plt.plot(K, inertia, marker='o')
    plt.xlabel('Number of clusters (k)')
    plt.ylabel('Total inertia')
    plt.title('Elbow method')
    plt.grid(True)
    elbow_path = os.path.join(output_dir, 'elbow_clusters.png')
    plt.savefig(elbow_path)
    print(f"Elbow plot saved to {elbow_path}")

    # Usa KneeLocator para encontrar el codo automáticamente
    kneedle = KneeLocator(K, inertia, curve='convex', direction='decreasing')
    best_k = kneedle.knee
    print(f"Optimal number of clusters (based on elbow): {best_k}")

    kmeans = KMeansConstrained(n_clusters=best_k, size_min=min_cluster_size, random_state=42)
    cluster_labels = kmeans.fit_predict(coords)

    df['Cluster'] = cluster_labels

    output_csv = os.path.join(output_dir, 'auto_populations.csv')
    df.to_csv(output_csv, index=False)
    print(f"{GREEN}✅ Inferred populations saved to {output_csv}{RESET}")

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

