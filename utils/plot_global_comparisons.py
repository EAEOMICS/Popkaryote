import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_population_comparisons(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    statistics = {
        'pi': ('pi.csv', 'pi_zscore', 'Nucleotide Diversity (π) - Z-score'),
        'tajima': ('tajima_d.csv', 'tajima_d_zscore', "Tajima's D - Z-score"),
        'WT': ('Watterson_theta_global.csv', 'Watterson_theta_zscore', 'Watterson Theta - Z-score')
    }

    upper_percentile = 2.326  # Z-score for the 99th percentile

    populations = [d for d in os.listdir(input_dir)
                   if os.path.isdir(os.path.join(input_dir, d))]
    populations = sorted(populations)

    for key, (filename, colname, main_title) in statistics.items():
        fig, axs = plt.subplots(len(populations), 1, figsize=(10, 4 * len(populations)), sharex=True)

        if len(populations) == 1:
            axs = [axs]  # ensure list if only one population

        for i, pop in enumerate(populations):
            path = os.path.join(input_dir, pop, filename)
            if not os.path.exists(path):
                axs[i].set_title(f'{pop} (no data)')
                axs[i].axis('off')
                continue

            df = pd.read_csv(path)
            if 'start' not in df.columns or colname not in df.columns:
                axs[i].set_title(f'{pop} (missing columns)')
                axs[i].axis('off')
                continue

            axs[i].plot(df['start'], df[colname], label=key.upper())
            axs[i].axhline(y=upper_percentile, color='red', linestyle='--', label='Z = +2.326')
            axs[i].set_title(pop)
            axs[i].set_ylabel('Z-score')
            axs[i].legend()

        plt.tight_layout()
        plt.suptitle(main_title, fontsize=16, y=1.02)
        plt.savefig(os.path.join(output_dir, f'{key}_comparison.png'), bbox_inches='tight')
        plt.close()
