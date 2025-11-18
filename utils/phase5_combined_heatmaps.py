import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_combined_heatmaps(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Dictionary with CSV file name and corresponding Z-score column
    csv_files = {
        'pi_zscore': ('pi.csv', 'pi_zscore'),
        'tajima_zscore': ('tajima_d.csv', 'tajima_d_zscore'),
        'ld_zscore': ('ld.csv', 'ld_r2_zscore'),
        'hd_zscore': ('hd_by_window.csv', 'hd_zscore')
    }

    statistics = list(csv_files.keys())
    populations = [d for d in os.listdir(input_dir)
                   if os.path.isdir(os.path.join(input_dir, d))]
    populations = sorted(populations)

    for pop in populations:
        data = []
        positions = None

        for stat, (filename, column) in csv_files.items():
            path = os.path.join(input_dir, pop, filename)
            if not os.path.exists(path):
                continue

            df = pd.read_csv(path)
            if 'start' not in df.columns or column not in df.columns:
                continue

            if positions is None:
                positions = df['start']

            data.append(df[column].values)

        if data and positions is not None:
            matrix = pd.DataFrame(data, index=statistics, columns=positions)
            plt.figure(figsize=(14, 3))
            sns.heatmap(matrix, cmap="coolwarm", center=0, cbar_kws={'label': 'Z-score'})
            plt.title(f"Combined Z-score heatmap - Population {pop}")
            plt.xlabel("Genomic position (start)")
            plt.ylabel("Statistic")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"{pop}_combined_heatmap.png"))
            plt.close()

if __name__ == "__main__":
    plot_combined_heatmaps("population_results", "combined_heatmaps")
