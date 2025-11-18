import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from scipy.spatial.distance import squareform
import subprocess
import re

GREEN = "\033[1;32m"
RESET = "\033[0m"

#default window size=10,000

def compute_global_statistics(vcf_path, output_dir,reference, window_size=10000,step=None,threads=4,dates=None,month=True):

    os.makedirs(output_dir, exist_ok=True)

    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    #NEW for haploids, remove fake fased chromosome (-1 entries in variants)
    haps=genotypes.to_haplotypes()
    haps_real=haps[:,:genotypes.shape[1]]
    genotypes=haps_real
    
    positions = callset['variants/POS']

    sorted_indices = np.argsort(positions)
    positions = positions[sorted_indices]
    genotypes = genotypes.take(sorted_indices, axis=0)


    if step is None:
        window_starts = np.arange(positions.min(), positions.max(), window_size)
        windows = np.column_stack((window_starts, window_starts + window_size))
    else:
        window_starts = np.arange(positions.min(), positions.max()- window_size + step, step)
        windows = np.column_stack((window_starts, window_starts + window_size))
        
    ac = genotypes.count_alleles()

    # π
    print(f"\n▶️ Calculating window-based π...")
    
    pi_total, _, _, n_sites = allel.windowed_diversity(pos=positions, ac=ac, size=window_size, step=step)
    pi = pi_total / n_sites #Pi per variant site
    pi_z = (pi - np.nanmean(pi)) / np.nanstd(pi)
    pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "pi": pi, "pi_zscore": pi_z}).to_csv(os.path.join(output_dir, 'pi_global.csv'), index=False)
    mean=np.nanmean(pi_total)
    plt.figure()
    plt.plot(windows[:, 0], pi_total)
    plt.axhline(y=mean, color='red', linestyle='--', label='mean')
    plt.xlabel('Position (bp)')
    plt.ylabel('π value')
    plt.title('Global π')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_window_global.png'))
    
    mean=np.nanmean(pi)
    #pi raw values
    plt.figure()
    plt.plot(windows[:, 0], pi)
    plt.axhline(y=mean, color='red', linestyle='--', label='mean')
    plt.xlabel('Position (bp)')
    plt.ylabel('π value')
    plt.title('Global π per variant site')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_raw_global.png'))
    
    #pi zscore
    plt.figure()
    plt.plot(windows[:, 0], pi_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Position (bp)')
    plt.ylabel('π (z-score)')
    plt.title('Global π Normalized')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_global.png'))
    
    mean=np.nanmean(pi_total)
    # Barplot de pi_total
    plt.figure(figsize=(12, 6))
    plt.bar(windows[:, 0], pi_total, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
    plt.axhline(y=mean, color='red', linestyle='--', label='mean')
    plt.xlabel('Position (bp)')
    plt.ylabel('π value')
    plt.title('Global π (Barplot)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_window_global_barplot.png'))
    
    mean=np.nanmean(pi)
    # Barplot de pi (raw)
    plt.figure(figsize=(12, 6))
    plt.bar(windows[:, 0], pi, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
    plt.axhline(y=mean, color='red', linestyle='--', label='mean')
    plt.xlabel('Position (bp)')
    plt.ylabel('π value')
    plt.title('Global π Raw (Barplot)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_raw_global_barplot.png'))

    # Barplot de pi_z (z-scores)
    plt.figure(figsize=(12, 6))
    plt.bar(windows[:, 0], pi_z, width=windows[1, 0] - windows[0, 0], align='center', edgecolor='black')
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Position (bp)')
    plt.ylabel('π (z-score)')
    plt.title('Global π Normalized (Barplot)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pi_global_barplot.png'))
    
    print(f"{GREEN}✅ Global π calculated{RESET}")
    
    #Tajima's D
    print("\n▶️ Calculating window-based Tajima's D...")
    tajima_d, _, _= allel.windowed_tajima_d(pos=positions, ac=ac, size=window_size, step=step)
    tajima_z = (tajima_d - np.nanmean(tajima_d)) / np.nanstd(tajima_d)
    df=pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "tajima_d": tajima_d, "tajima_d_zscore": tajima_z})
    
    tajima_d_global = allel.tajima_d(ac)
    print("Tajima's D global:", tajima_d_global)
    df = pd.concat([
      df,
      pd.DataFrame({
        "start": [1],                    
        "end": [positions[-1]],         
        "tajima_d": [tajima_d_global],
        "tajima_d_zscore": [np.nan]      
      })
    ], ignore_index=True)
    
    df.to_csv(os.path.join(output_dir, 'tajima_d_global.csv'), index=False)
    
    plt.figure()
    plt.plot(windows[:, 0], tajima_d)
    plt.axhline(y=0, color='red', linestyle='--', label='D = 0')
    plt.axhline(y=tajima_d_global, color='black', linestyle='--', label='global')
    plt.xlabel('Position (bp)')
    plt.ylabel("Tajima's D")
    plt.title("Global Tajima's D")
    max_abs = np.nanmax(np.abs(tajima_d))
    plt.ylim(-max_abs, max_abs)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tajima_raw_global.png'))
    
    plt.figure()
    plt.plot(windows[:, 0], tajima_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Position (bp)')
    plt.ylabel("Tajima's D (z-score)")
    plt.title("Global Tajima's D (normalized)")
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'tajima_global.png'))
    
    print(f"{GREEN}✅ Global Tajima's D calculated{RESET}")
    # LD
#    gn_alt = genotypes.to_n_alt()
#    mask = ~np.any(np.isnan(gn_alt), axis=1)
#    gn_alt = gn_alt[mask]
#    pos_filtered = positions[mask]

#    ld_means = []
#    for start, end in windows:
#        mask_window = (pos_filtered >= start) & (pos_filtered < end)
#        gn_window = gn_alt[mask_window, :]
#        if gn_window.shape[0] < 2:
#            ld_means.append(np.nan)
#        else:
#            r = allel.rogers_huff_r(gn_window)
#            r2 = r ** 2
#            ld_means.append(np.nanmean(r2))

#    ld_z = (np.array(ld_means) - np.nanmean(ld_means)) / np.nanstd(ld_means)
#    pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "ld_r2": ld_means, "ld_r2_zscore": ld_z}).to_csv(os.path.join(output_dir, 'ld_global.csv'), index=False)
    
#    plt.figure()
#    plt.plot(windows[:, 0], ld_means)
#    plt.xlabel('Position (bp)')
#    plt.ylabel('LD (r²)')
#    plt.title('Global LD')
#    plt.tight_layout()
#    plt.savefig(os.path.join(output_dir, 'ld_raw_global.png'))
    
#    plt.figure()
#    plt.plot(windows[:, 0], ld_z)
#    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
#    plt.xlabel('Position (bp)')
#    plt.ylabel('LD (r² z-score)')
#    plt.title('Global LD (normalized)')
#    plt.tight_layout()
#    plt.savefig(os.path.join(output_dir, 'ld_global.png'))

#    print(f"\n▶️ Calculating window-based Linakge-Disequilibrium...")
#    
#    #Transform Missing in Diploid Genotype Array to All Homozygot Genotype Array
#    gts=allel.GenotypeArray(callset['calldata/GT'])
#    gts=gts.take(sorted_indices, axis=0)
#    mask_a2 = gts[:, :, 1] == -1
#    gts[:, :, 1][mask_a2] = gts[:, :, 0][mask_a2]
#    gn = gts.to_n_alt(fill=-1)
#    r = allel.rogers_huff_r(gn)
#    r2=squareform(r ** 2)
#    
#    mask = np.triu(np.ones_like(r2, dtype=bool))
#    plt.figure(figsize=(16, 12))
#    sns.heatmap(
#        r2,
#        mask=mask,
#        cmap='Reds',
#        square=True,
#        vmin=0, vmax=1,
#        cbar_kws={'label': '$r^2$'}
#    )
#    
#    plt.title('Linkage Disequilibrium (LD) Heatmap')
#    plt.xlabel('Variant')
#    plt.ylabel('Variant')
#    plt.savefig(os.path.join(output_dir, 'LD_heatmap.png'))
#
#    print(f"{GREEN}✅ Global Linkage Disequilibrium calculated{RESET}")
    
    # Hd (haplotype diversity)

    haps=genotypes
    
    print("Number of SNPs without missing data for Hd:", haps.shape[0])
    with open(os.path.join(output_dir, 'hd_global.txt'), 'w') as f:
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
    plt.title('Haplotype diversity (variant_size 50/100/200)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,'Hd_50_100_200.png'))
    
    print(f"{GREEN}✅ Global Window-based Hd calculated{RESET}")
    
    
    #Watterson's Theta
    print("\n▶️ Calculating window-based Watterson's Theta...")
    WT, _, _, _= allel.windowed_watterson_theta(pos=positions, ac=ac, size=window_size, step=step)
    WT_z = (WT - np.nanmean(WT)) / np.nanstd(WT)
    df_wt=pd.DataFrame({"start": windows[:,0], "end": windows[:,1], "Watterson_theta": WT, "Watterson_theta_zscore": WT_z})
    
    watterson_theta_global = allel.watterson_theta(pos=positions, ac=ac, start=1, stop=positions[-1])
    df_wt = pd.concat([
      df_wt,
      pd.DataFrame({
        "start": [1],                           
        "end": [positions[-1]],                 
        "Watterson_theta": [watterson_theta_global],
        "Watterson_theta_zscore": [np.nan]      
      })
    ], ignore_index=True)
    
    df_wt.to_csv(os.path.join(output_dir, 'Watterson_theta_global.csv'), index=False)
    
    plt.figure()
    plt.plot(windows[:, 0], WT)
    plt.axhline(y=watterson_theta_global, color='black', linestyle='--', label='global')
    plt.xlabel('Position (bp)')
    plt.ylabel("Watterson's Theta")
    plt.title("Global Watterson's Theta")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'WT_raw_global.png'))
    
    plt.figure()
    plt.plot(windows[:, 0], WT_z)
    plt.axhline(y=2.326, color='red', linestyle='--', label='Z = +2.326')
    plt.xlabel('Position (bp)')
    plt.ylabel("Watterson's Theta (z-score)")
    plt.title("Global Watterson's Theta (normalized)")
    plt.tight_layout()
    plt.legend()
    plt.savefig(os.path.join(output_dir, 'WT_global.png'))
    
    print(f"{GREEN}✅ Global Window-based Watterson's Theta calculated{RESET}")


    #Bactdating
    print("\n▶️ Dating the nodes of the bacterial phylogenetic tree...")
    output_dir=output_dir.rsplit("global", 1)[0]
    iqtree=f"iqtree -s {output_dir}core.full.aln -T {threads} -B 10000 --seqtype DNA --use-nn-model --prefix {output_dir}global/iqtree_output"
    subprocess.run(iqtree, shell=True, executable='/bin/bash')
    clonalF=f"ClonalFrameML {output_dir}global/iqtree_output.treefile {output_dir}core.full.aln {output_dir}global/cf_out"
    subprocess.run(clonalF,shell=True,executable='/bin/bash')
    clonalF2=f"Rscript utils/cfml_results.R {output_dir}global/cf_out"
    subprocess.run(clonalF2,shell=True,executable='/bin/bash')
    print(f"{GREEN}✅ Bacterial phylogenetic tree constructed {RESET}")
    if dates is not None:
        if  month == 'True':
            BactDating=f"Rscript utils/BACTDATING.R {dates} {output_dir}global True "
        else:
            BactDating=f"Rscript utils/BACTDATING.R {dates} {output_dir}global False "
        subprocess.run(BactDating,shell=True,executable='/bin/bash')
        print(f"{GREEN}✅ bacterial phylogenetic tree nodes dated{RESET}")
    

    print("\n▶️ Performing Linkage Disequilibrium study with LDWEAVER...")
    LDWEAVER=f"Rscript utils/LDWEAVER.R {reference} {output_dir}"
    subprocess.run(LDWEAVER,shell=True,executable='/bin/bash')
    print(f"{GREEN}✅ LDWEAVER Done...{RESET}")
    print(f"{GREEN}✅ Global statistics calculated and saved to:{RESET}", output_dir)