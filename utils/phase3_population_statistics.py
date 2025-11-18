import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pickle
from scipy.spatial.distance import squareform
import seaborn as sns

GREEN = "\033[1;32m"
RESET = "\033[0m"

def compute_population_statistics(vcf_path, populations_csv, output_dir, window_size=10000,step=None):

    os.makedirs(output_dir, exist_ok=True)

    df_pops = pd.read_csv(populations_csv)
    samples = df_pops['Sample'].tolist()
    clusters = df_pops['Cluster'].tolist()
    populations = sorted(df_pops['Cluster'].unique())
    print(f"Detected populations: {populations}")

    callset = allel.read_vcf(vcf_path, samples=samples)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    genotypes2 = allel.GenotypeArray(callset['calldata/GT'])
    
    haps=genotypes.to_haplotypes()
    haps_real=haps[:,:genotypes.shape[1]]
    genotypes=haps_real
    
    positions = callset['variants/POS']

    sorted_idx = np.argsort(positions)
    positions = positions[sorted_idx]
    genotypes = genotypes.take(sorted_idx, axis=0)

    sample_names = callset['samples']
    df_sample_cluster = pd.DataFrame({"sample": sample_names, "cluster": clusters})

    subpops = [df_sample_cluster[df_sample_cluster['cluster'] == pop].index.tolist() for pop in populations]
    with open(os.path.join(output_dir, 'subpops_indices.pkl'), 'wb') as f:
        pickle.dump(subpops, f)

    if step is None:
        window_starts = np.arange(positions.min(), positions.max(), window_size)
        windows = np.column_stack((window_starts, window_starts + window_size))
    else:
        window_starts = np.arange(positions.min(), positions.max()- window_size + step, step)
        windows = np.column_stack((window_starts, window_starts + window_size))

    for i, pop in enumerate(populations):
        if len(subpops[i]) < 2:
            print(f"⚠️  Skipping population {pop} (only populations with at least 2 samples))")
            continue
        pop_dir = os.path.join(output_dir, f"population_{pop}")
        os.makedirs(pop_dir, exist_ok=True)

        gen_pop = genotypes.take(subpops[i], axis=1)
        ac_pop = gen_pop.count_alleles()
                
        # π
        
        print(f"\n▶️ Calculating window-based π for population {pop}...")
        
        pi_total, _, _, n_sites = allel.windowed_diversity(pos=positions, ac=ac_pop, size=window_size, step=step)
        pi = pi_total / n_sites
        pi_z = (pi - np.nanmean(pi)) / np.nanstd(pi)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "pi": pi, "pi_zscore": pi_z}).to_csv(os.path.join(pop_dir, 'pi.csv'), index=False)
        
        
        mean=np.nanmean(pi_total)
        plt.figure()
        plt.plot(windows[:, 0], pi)
        plt.xlabel('Window start (bp)')
        plt.ylabel('π')
        plt.title(f'Nucleotide diversity - Population {pop}')
        plt.axhline(y=mean, color='red', linestyle='--', label='mean')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi_raw.png'))
        plt.close()
        
        #pi raw values
        mean=np.nanmean(pi)
        plt.figure()
        plt.plot(windows[:, 0], pi)
        plt.axhline(y=mean, color='red', linestyle='--', label='mean')
        plt.xlabel('Position (bp)')
        plt.ylabel('π value')
        plt.title(f'Population {pop} π per variant site')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi_raw_global.png'))
        
        #pi zscore
        plt.figure()
        plt.plot(windows[:, 0], pi_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Window start (bp)')
        plt.ylabel('π (z-score)')
        plt.title(f'Nucleotide diversity - Population {pop}')
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi.png'))
        plt.close()
        
        mean=np.nanmean(pi_total)
        
        # Barplot de pi_total
        plt.figure(figsize=(12, 6))
        plt.bar(windows[:, 0], pi_total, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
        plt.axhline(y=mean, color='red', linestyle='--', label='mean')
        plt.xlabel('Position (bp)')
        plt.ylabel('π value')
        plt.title(f'Population {pop} π (Barplot)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi_window_barplot.png'))
        
        mean=np.nanmean(pi)
        # Barplot de pi (raw)
        plt.figure(figsize=(12, 6))
        plt.bar(windows[:, 0], pi, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
        plt.axhline(y=mean, color='red', linestyle='--', label='mean')
        plt.xlabel('Position (bp)')
        plt.ylabel('π value')
        plt.title(f'Population {pop} π Raw (Barplot)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi_rawl_barplot.png'))
        
        # Barplot de pi_z (z-scores)
        plt.figure(figsize=(12, 6))
        plt.bar(windows[:, 0], pi_z, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Position (bp)')
        plt.ylabel('π (z-score)')
        plt.title(f'Population {pop} π Normalized (Barplot)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'pi_barplot.png'))
        
        print(f"{GREEN}✅ Population's {pop} π calculated{RESET}")
        
        # Tajima's D
        
        print(f"\n▶️ Calculating window-based Tajima's D for population {pop}...")
        
        #if len(pos_filtered) < 2:
        #    tajima_d = np.array([np.nan] * len(windows))
        #else:
        #    ac_filtered = gen_filtered.count_alleles()
        #    tajima_d, _, _ = allel.windowed_tajima_d(pos=pos_filtered, ac=ac_pop_filtered, windows=windows)
        
        tajima_d, _, _= allel.windowed_tajima_d(pos=positions, ac=ac_pop, size=window_size, step=step)
        tajima_z = (tajima_d - np.nanmean(tajima_d)) / np.nanstd(tajima_d)
        pd.DataFrame({"start": windows[:, 0], "end": windows[:, 1], "tajima_d": tajima_d, "tajima_d_zscore": tajima_z}).to_csv(os.path.join(pop_dir, 'tajima_d.csv'), index=False)
        
        mean=np.nanmean(tajima_d)
        plt.figure()
        plt.plot(windows[:, 0], tajima_d)
        plt.axhline(y=0, color='red', linestyle='--', label='D = 0')
        plt.axhline(y=mean, color='black', linestyle='--', label='mean')
        plt.xlabel('Position (bp)')
        plt.ylabel("Tajima's D")
        plt.title(f"Tajima's D - Population {pop}")
        max_abs = np.nanmax(np.abs(tajima_d))
        plt.ylim(-max_abs, max_abs)
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'tajima_raw.png'))
        plt.close()
        
        plt.figure()
        plt.plot(windows[:, 0], tajima_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Position (bp)')
        plt.ylabel("Tajima's D (z-score)")
        plt.title(f"Tajima's D - Population {pop}")
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'tajima.png'))
        plt.close()
        
        print(f"{GREEN}✅ Population's {pop} Tajima's D calculated{RESET}")

        # Hd 
        
        print(f"\n▶️ Calculating window-based Hd for population {pop}...")
        
        haps=gen_pop
        print("Number of SNPs without missing data for Hd:", haps.shape[0])
        with open(os.path.join(pop_dir, 'hd_global.txt'), 'w') as f:
            if haps.shape[0] == 0:
                print("⚠️  Global haplotype diversity (Hd) could not be calculated.")
                print("ℹ️  This happens when there are not enough polymorphic SNPs without missing data across ALL samples.")
                print("ℹ️  Consider checking per-population Hd results instead.")
                f.write("⚠️  Global haplotype diversity (Hd) could not be calculated.\n")
                f.write("ℹ️  Possibly due to insufficient polymorphic SNPs without missing data across all samples.\n")
                f.write("ℹ️  Check population-level results where haplotypic variation may exist.\n")
            else:
                num_haplotypes = np.unique(haps.T, axis=0).shape[0]
                print("Number of unique haplotypes:", num_haplotypes)
                hd_global = allel.haplotype_diversity(haps)
                f.write(f"Average haplotype diversity: {hd_global}\n")
                f.write(f"Number of unique haplotypes: {num_haplotypes}\n")

        # Window-based Hd
        print("\n▶️ Calculating window-based Hd...")
    
        sizes = [50, 100, 200]
        colors = ['blue', 'red', 'black']
        labels = ['size=50 / step 25', 'size=100 / step 50', 'size=200 / step 100']
        plt.figure(figsize=(12, 5))

        for Size, color, label in zip(sizes, colors, labels):
            stup = Size // 2  # step = size/2
            hd = allel.moving_haplotype_diversity(haps, size=Size, step=stup)
            starts = np.arange(0, len(positions) - Size + 1, stup)
            # Para cada ventana, calcula la posición media de los SNPs dentro
            window_pos = np.array([positions[s:s+Size].mean() for s in starts])
            plt.plot(window_pos, hd, color=color, label=label)

        plt.xlabel('Genomic position (bp)')
        plt.ylabel('Haplotype diversity')
        plt.title(f'Haplotype diversity (variant_size 50/100/200) - Population {pop}')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir,'Hd_50_100_200.png'))
        
        print(f"{GREEN}✅ Population's {pop} Hd calculated{RESET}")

        # LD
        
#        print(f"\n▶️ Calculating window-based Linkage Disequilibrium for population {pop}...")
#        
#        gts=genotypes2.take(subpops[i], axis=1)
#        mask_a2 = gts[:, :, 1] == -1
#        gts[:, :, 1][mask_a2] = gts[:, :, 0][mask_a2]
#        gn = gts.to_n_alt(fill=-1)
#        r = allel.rogers_huff_r(gn)
#        r2=squareform(r ** 2)
#        
#        mask = np.triu(np.ones_like(r2, dtype=bool))
#        plt.figure(figsize=(16, 12))
#        sns.heatmap(
#          r2,
#          mask=mask,
#          cmap='Reds',
#          square=True,
#          vmin=0, vmax=1,
#          cbar_kws={'label': '$r^2$'}
#        )
#        
#        plt.title(f'Linkage Disequilibrium (LD) Heatmap - Population {pop}')
#        plt.xlabel('Variant')
#        plt.ylabel('Variant')
#        plt.savefig(os.path.join(pop_dir, 'LD_heatmap.png'))
#        
#        print(f"{GREEN}✅ Population's {pop} Linkage Disequilibrium calculated{RESET}")
                
        #Watterson's Theta
        
        print(f"\n▶️ Calculating Watterson's Theta for population {pop}...")
        
        WT, _, _, _= allel.windowed_watterson_theta(pos=positions, ac=ac_pop, size=window_size, step=step)
        WT_z = (WT - np.nanmean(WT)) / np.nanstd(WT)
        pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "Watterson_theta": WT, "Watterson_theta_zscore": WT_z}).to_csv(os.path.join(pop_dir, 'Watterson_theta_global.csv'), index=False)
        
        mean=np.nanmean(WT)
        plt.figure()
        plt.plot(windows[:, 0], WT)
        plt.axhline(y=mean, color='black', linestyle='--', label='mean')
        plt.xlabel('Position (bp)')
        plt.ylabel("Watterson's Theta")
        plt.title(f"Global Watterson's Theta - Population {pop}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(pop_dir, 'WT_raw.png'))
        
        plt.figure()
        plt.plot(windows[:, 0], WT_z)
        plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
        plt.xlabel('Position (bp)')
        plt.ylabel("Watterson's Theta (z-score)")
        plt.title(f"Global Watterson's Theta (normalized) - Population {pop}")
        plt.tight_layout()
        plt.legend()
        plt.savefig(os.path.join(pop_dir, 'WT_normalized.png'))
        
        print(f"{GREEN}✅ Population's {pop} Watterson's Theta calculated{RESET}")

