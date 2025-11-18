import os
import numpy as np
import pandas as pd
import allel
import matplotlib.pyplot as plt
from itertools import combinations

def compute_pairwise_fst(vcf_path, populations_csv, output_dir, window_size=50000, min_samples=3):
    os.makedirs(output_dir, exist_ok=True)

    print("\n[PHASE 4] Computing pairwise FST between populations...")

    # Load data
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    positions = callset['variants/POS']
    samples = callset['samples']

    df = pd.read_csv(populations_csv)
    df.columns = df.columns.str.lower()

    # Group samples by cluster
    cluster_dict = {
        str(cluster): df[df['cluster'] == cluster]['sample'].tolist()
        for cluster in df['cluster'].unique()
    }

    # Filter clusters with enough samples
    subpops = []
    labels = []
    for cluster, cluster_samples in cluster_dict.items():
        indices = [np.where(samples == s)[0][0] for s in cluster_samples if s in samples]
        if len(indices) >= min_samples:
            subpops.append(indices)
            labels.append(f"Pop{cluster}")
        else:
            print(f"⚠️  Skipping cluster {cluster}: only {len(indices)} sample(s)")

    if len(subpops) < 2:
        print("❌ Not enough valid populations for pairwise comparison.")
        return

    window_starts = np.arange(positions.min(), positions.max(), window_size)
    windows = np.column_stack((window_starts, window_starts + window_size))

    all_fst_results = []
    pair_labels = []

    for (i, j) in combinations(range(len(subpops)), 2):
        ac1 = genotypes[:, subpops[i]].count_alleles()
        ac2 = genotypes[:, subpops[j]].count_alleles()

        fst_vals = []
        for start, end in windows:
            mask = (positions >= start) & (positions < end)
            if np.sum(mask) < 2:
                fst_vals.append(np.nan)
                continue
            num, den = allel.hudson_fst(ac1[mask], ac2[mask])
            with np.errstate(invalid='ignore', divide='ignore'):
                fst_window = np.sum(num) / np.sum(den)
            fst_vals.append(fst_window)

        pair_label = f"{labels[i]} vs {labels[j]}"
        pair_labels.append(pair_label)
        all_fst_results.append(fst_vals)

        df_result = pd.DataFrame({
            "start": windows[:, 0],
            "end": windows[:, 1],
            "fst": fst_vals
        })
        df_result.to_csv(os.path.join(output_dir, f"{pair_label.replace(' ', '_')}_fst_windows.csv"), index=False)

        # Individual line plot per pair
        plt.figure(figsize=(10, 4))
        plt.plot(windows[:, 0], fst_vals, marker='o', linestyle='-', markersize=2)
        plt.title(f"Windowed FST: {pair_label}")
        plt.xlabel("Genomic position")
        plt.ylabel("FST (Hudson)")
        plt.ylim(0, 1)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{pair_label.replace(' ', '_')}_fst_line.png"))
        plt.close()

    # Global comparison plot - lines
    fig, axs = plt.subplots(len(all_fst_results), 1, figsize=(10, 3.5 * len(all_fst_results)), sharex=True)
    if len(all_fst_results) == 1:
        axs = [axs]
    for i, (label, values) in enumerate(zip(pair_labels, all_fst_results)):
        axs[i].plot(windows[:, 0], values, marker='o', linestyle='-', markersize=2)
        axs[i].set_title(label)
        axs[i].set_ylabel("FST")
        axs[i].set_ylim(0, 1)
    axs[-1].set_xlabel("Genomic position")
    plt.suptitle("Windowed FST comparison between population pairs", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fst_comparative_lines.png"), bbox_inches="tight")
    plt.close()

    # Global comparison - boxplots
    plt.figure(figsize=(10, 5))
    plt.boxplot([pd.Series(vals).dropna() for vals in all_fst_results], labels=pair_labels, patch_artist=True)
    plt.ylabel("FST")
    plt.title("Distribution of FST values between population pairs")
    plt.ylim(0, 1)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "fst_comparative_boxplots.png"))
    plt.close()

    # Summary table of average FST per comparison
    summary = []
    for label, values in zip(pair_labels, all_fst_results):
        mean_fst = pd.Series(values).mean(skipna=True)
        summary.append({'Comparison': label, 'Average_FST': round(mean_fst, 5)})

    df_summary = pd.DataFrame(summary).sort_values(by='Average_FST', ascending=False)
    df_summary.to_csv(os.path.join(output_dir, "fst_average_table.csv"), index=False)
    print("📊 Summary table saved: fst_average_table.csv")

    # FST heatmap
    try:
        import seaborn as sns

        unique_pops = sorted(set(sum([label.split(" vs ") for label in pair_labels], [])))
        matrix = pd.DataFrame(np.nan, index=unique_pops, columns=unique_pops)

        for label, values in zip(pair_labels, all_fst_results):
            pop1, pop2 = label.split(" vs ")
            mean_val = pd.Series(values).mean(skipna=True)
            matrix.loc[pop1, pop2] = mean_val
            matrix.loc[pop2, pop1] = mean_val  # symmetric

        plt.figure(figsize=(8, 6))
        sns.heatmap(matrix, annot=True, fmt=".3f", cmap="YlGnBu", vmin=0, vmax=1)
        plt.title("Average FST heatmap between populations")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "fst_average_heatmap.png"))
        plt.close()

        print("Average FST heatmap saved: fst_average_heatmap.png")
    except ImportError:
        print("⚠️  Could not generate FST heatmap (seaborn not installed)")

    print("✅ Windowed FST successfully calculated and visualized for all pairs.")

if __name__ == "__main__":
    compute_pairwise_fst("path/to/file.vcf", "path/to/populations.csv", "output/fst_by_pairs")
