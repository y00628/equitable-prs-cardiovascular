import pandas as pd
import numpy as np
from pathlib import Path
import os
from functools import reduce
import regex as re
from venny4py.venny4py import *
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os

def str_to_float(s):
    try:
        f = float(s)
        return f
    except Exception as e:
        return np.nan


def get_all_dat(directory):
  """
  Gets a list of all .dat files in a directory and its subdirectories.

  Args:
    directory: The path to the directory.

  Returns:
    A list of file paths.
  """
  directory = Path(directory)

  # Search recursively for .dat files
  dat_files = list(directory.rglob("*.dat"))
  return dat_files

def gene_loci():
    assoc_test_fps = get_all_dat(os.path.join('..', 'data', 'gene_associations'))

    all_dfs = []
    for fp in assoc_test_fps:
        temp_df = pd.read_table(fp).drop(columns=['PANEL', 'FILE'])
        population = re.findall(r'chr\d+.(\w{3})', str(fp))[0]
        temp_df['TWAS.P'] = temp_df['TWAS.P'].apply(str_to_float)
        temp_df['TWAS.Z'] = temp_df['TWAS.Z'].apply(str_to_float)
        temp_df['POP'] = population
        temp_df = temp_df[(temp_df['TWAS.P'].notna()) & (temp_df['TWAS.Z'].notna())]
        temp_df = temp_df[temp_df['TWAS.P'] <= 0.05]
        all_dfs.append(temp_df)
    all_chrs_pops = pd.concat(all_dfs, ignore_index=True)


    eur_gene_loci = all_chrs_pops[(all_chrs_pops['POP'] == 'eur')]
    eas_gene_loci = all_chrs_pops[(all_chrs_pops['POP'] == 'eas')]
    amr_gene_loci = all_chrs_pops[(all_chrs_pops['POP'] == 'amr')]
    afr_gene_loci = all_chrs_pops[(all_chrs_pops['POP'] == 'afr')]

    populations = {
        'EUR': eur_gene_loci,
        'EAS': eas_gene_loci,
        'AMR': amr_gene_loci,
        'AFR': afr_gene_loci,
        
    }

    # Get chromosome boundaries from all data
    all_data = pd.concat(populations.values())
    chrom_boundaries = all_data.groupby("CHR").agg({"P0": "min", "P1": "max"}).reset_index()
    chrom_boundaries.rename(columns={"P0": "Min", "P1": "Max"}, inplace=True)

    # Calculate the min and max TWAS.P values across all populations
    min_pval = all_data["TWAS.P"].min()
    max_pval = all_data["TWAS.P"].max()

    # Create a figure with subplots for each chromosome horizontally
    n_chromosomes = len(chrom_boundaries)
    fig, axes = plt.subplots(1, n_chromosomes, figsize=(30, 4), sharey=True)

    # Plot for each population and chromosome
    for j, (chrom_idx, row) in enumerate(chrom_boundaries.iterrows()):
        ax = axes[j]  # Get the corresponding axis for the chromosome

        # Plot reference genome lines for all populations on this chromosome
        for i, (pop, df) in enumerate(populations.items()):
            ypos = i + 1  # The y-position is the row number for the population
            chrom_genes = df[df["CHR"] == row["CHR"]]
            
            for _, gene in chrom_genes.iterrows():
                # Normalize the TWAS.P value to get opacity (alpha) between 0 and 1
                p_value = gene["TWAS.P"]
                alpha = (p_value - min_pval) / (max_pval - min_pval)  # Normalized p-value between 0 and 1
                
                # Ensure alpha stays between 0 and 1
                alpha = min(max(alpha, 0), 1)

                # Plot the gene as an arrow with the corresponding opacity
                arrow = patches.FancyArrow(
                    gene["P0"], ypos, gene["P1"] - gene["P0"], 0,
                    width=0.05, head_width=0.15, head_length=20, 
                    color="blue", alpha=alpha,  # Apply opacity here
                    length_includes_head=True
                )
                ax.add_patch(arrow)

        # Set x-axis limits for this chromosome
        ax.set_xlim(row["Min"] - 100, row["Max"] + 100)  # Add some padding around the gene region
        ax.set_ylim(0, len(populations) + 1)  # Keep a consistent y-axis range

        # Label the chromosomes (set to be shown on all subplots)
        ax.set_xlabel(f"Chr {row['CHR']}", fontsize=10)
            
        # Remove x-axis ticks (hide the scale)
        ax.set_xticks([])

        # Remove y-axis ticks for clarity
        ax.set_yticks([])

        # Keep only the left and bottom spines visible, remove top and right spines
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        if j != 0:  # Remove left spine for all but the first subplot
            ax.spines['left'].set_visible(False)

    # Set y-axis labels only on the first subplot for populations with spacing
    for i, (pop, _) in enumerate(populations.items()):
        axes[0].text(-0.15, i + 1, pop, ha='right', va='center', fontsize=12)  # Adjust the x-position for spacing

    # Add custom legend for opacity
    # Create a color map and normalization for opacity based on p-values
    norm = mcolors.Normalize(vmin=min_pval, vmax=max_pval)
    cmap = mcolors.LinearSegmentedColormap.from_list("pval_opacity", ["blue", "white"])

    # Create a color bar (or legend) using the colormap
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Empty array for the color bar to work
    cax = fig.add_axes([0.87, 1, 0.1, 0.05])  # Adjust these values to fit your needs
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal', )
    cbar.set_label('TWAS P-value', rotation=0)
    cbar.set_ticks([0, 0.05])



    fig.suptitle('Multi-Population Gene Loci')
    #fig.supylabel('Population')
    fig.supxlabel('Chromosome')

    # Adjust layout and spacing between subplots
    plt.tight_layout(w_pad=0)
    print(f"Saved gene_loci plot to {os.path.join("..","plots","gene_loci.png")}")
    plt.savefig(os.path.join("..","plots","gene_loci.png"))



def generate_venn():
    assoc_test_fps = get_all_dat(os.path.join('..', 'data', 'gene_associations'))

    all_dfs = []
    for fp in assoc_test_fps:
        temp_df = pd.read_table(fp).drop(columns=['PANEL', 'FILE'])
        population = re.findall(r'chr\d+.(\w{3})', str(fp))[0]
        temp_df['TWAS.P'] = temp_df['TWAS.P'].apply(str_to_float)
        temp_df['TWAS.Z'] = temp_df['TWAS.Z'].apply(str_to_float)
        temp_df['POP'] = population
        temp_df = temp_df[(temp_df['TWAS.P'].notna()) & (temp_df['TWAS.Z'].notna())]
        temp_df = temp_df[temp_df['TWAS.P'] <= 0.05]
        all_dfs.append(temp_df)
    all_chrs_pops = pd.concat(all_dfs, ignore_index=True)
    eur_genes = set(all_chrs_pops[(all_chrs_pops['POP'] == 'eur')]['ID'].unique())
    eas_genes = set(all_chrs_pops[(all_chrs_pops['POP'] == 'eas')]['ID'].unique())
    amr_genes = set(all_chrs_pops[(all_chrs_pops['POP'] == 'amr')]['ID'].unique())
    afr_genes = set(all_chrs_pops[(all_chrs_pops['POP'] == 'afr')]['ID'].unique())

    gene_sets = {
        'EUR': eur_genes,
        'EAS': eas_genes,
        'AMR': amr_genes,
        'AFR': afr_genes,
        
    }
    out = os.path.join('..', 'plots', 'pop_gene_venn')
    venny4py(sets=gene_sets, out=out)
    print(f'Saved venn diagram to {out}')

def manhattan_helper(pop, p_threshold_significance=5e-8, subsample_cutoff=5e-1, subsample_rate=0.01):
    # Read the GWAS summary statistics file
    file_path = f'../data/raw/gwas/HF_Bothsex_{pop}_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz'
    df = pd.read_csv(file_path, compression='gzip', sep='\t', header=0)

    # Ensure required columns exist
    required_cols = {'#CHR', 'POS', 'inv_var_meta_p'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Missing required columns in input file: {required_cols - set(df.columns)}")

    # Rename columns for consistency
    df = df.rename(columns={'#CHR': 'CHR', 'inv_var_meta_p': 'P'})

    # Select relevant columns and drop missing values
    df = df[['CHR', 'POS', 'P']].dropna()
   
    # Convert P-values to -log10(P) for plotting
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    df = df.dropna()  # Drop any remaining NaNs
    df['-log10(P)'] = df['P'].apply(lambda x: -np.log10(x))

    # Convert chromosome numbers to integers for sorting
    df['CHR'] = df['CHR'].astype(int)

    # Apply filtering:
    df_significant = df[df['P'] <= subsample_cutoff]  # Keep all SNPs below threshold
    df_nonsignificant = df[df['P'] > subsample_cutoff].sample(frac=subsample_rate, random_state=42)  # Subsample 1% of non-significant SNPs

    # Combine both datasets
    df = pd.concat([df_significant, df_nonsignificant])

    # Sort by chromosome and position
    df = df.sort_values(by=['CHR', 'POS'])

    # Generate a Manhattan plot
    plt.figure(figsize=(20, 6))
    colors = ['blue', 'red']  # Alternating colors for chromosomes

    # Create chromosome index mapping
    df['Chrom_Position'] = df.groupby('CHR')['POS'].transform(lambda x: x - x.min())

    # Set up x-axis positions
    chrom_max = df.groupby('CHR')['Chrom_Position'].max()
    chrom_offsets = chrom_max.cumsum().shift(fill_value=0)
    df['Adjusted_POS'] = df['Chrom_Position'] + df['CHR'].map(chrom_offsets)

    # Define colors
    significant_color = 'red'
    gray_shades = ['darkgray', 'lightgray']

    # Plot each chromosome separately
    for i, chrom in enumerate(df['CHR'].unique()):
        chr_data = df[df['CHR'] == chrom]
    
        # Separate significant and non-significant SNPs
        chr_significant = chr_data[chr_data['P'] <= p_threshold_significance]
        chr_nonsignificant = chr_data[chr_data['P'] > p_threshold_significance]
    
        # Plot non-significant SNPs in alternating gray shades
        plt.scatter(chr_nonsignificant['Adjusted_POS'], chr_nonsignificant['-log10(P)'],
                    color=gray_shades[i % 2], s=5, alpha=0.6)
    
        # Plot significant SNPs in red
        plt.scatter(chr_significant['Adjusted_POS'], chr_significant['-log10(P)'],
                    color=significant_color, s=10)

    # Add significance threshold line
    plt.axhline(y=-np.log10(p_threshold_significance), color='black', linestyle='dashed', linewidth=1)

    # Format x-axis to show chromosome numbers
    tick_positions = chrom_offsets + chrom_max / 2
    plt.xticks(tick_positions, [str(chrom) for chrom in chrom_offsets.index])
    plt.ylim(0, 20)

    # Labels and title
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P-value)")
    plt.title(f"Manhattan Plot of GWAS Data (Population: {pop})")

    plt.savefig(f'../plots/manhattan/{pop}_manhattan_v4.png')

def manhattan():
    for pop in ['afr', 'eas', 'amr', 'eur']:
        print(f'Creating Manhattan plot for {pop}...')
        manhattan_helper(pop)
        print('Done!')
    
    
def data_prep(pop):
    """Reads GWAS summary statistics for all chromosomes and concatenates them."""
    all_dfs = []

    for chrom_num in range(1, 23):  # Loop through chromosomes 1-22
        file_path = f'../data/gene_associations/chr{chrom_num}.{pop}.assoc.dat'
        
        try:
            df = pd.read_csv(file_path, sep="\t")

            df = df.drop(columns=['PANEL', 'FILE'])

            # Select relevant columns
            df = df[['ID', 'CHR', 'P0', 'TWAS.P']].dropna()

            # Strip whitespace and convert to numeric, forcing errors to NaN
            df['TWAS.P'] = pd.to_numeric(df['TWAS.P'].str.strip(), errors='coerce')

            # Drop NaN values
            df = df.dropna()

            # Convert P-values to -log10(P) for plotting
            df['TWAS.P'] = df['TWAS.P'].astype(float)
            df['-log10(P)'] = df['TWAS.P'].apply(lambda x: -np.log10(x))    

            # Convert chromosome numbers to integers for sorting
            df['CHR'] = df['CHR'].astype(int)

            # Sort by chromosome and position
            df = df.sort_values(by=['CHR', 'P0'])

            all_dfs.append(df)  # Store processed df

        except FileNotFoundError:
            print(f"File not found: {file_path}, skipping chromosome {chrom_num}")

    # Concatenate all chromosome data
    final_df = pd.concat(all_dfs, ignore_index=True)

    return final_df


def miami_helper(pop):
    """Generates a Miami plot for chromosomes 1-22 comparing given population vs 'eur'."""
    
    # Get GWAS data for both populations across all chromosomes
    df1 = data_prep('eur')
    df2 = data_prep(pop)

    # Set up the plot
    plt.figure(figsize=(14, 6))
    sig_color = 'red'  # Color for significant SNPs
    grey_colors = ['lightgrey', 'darkgrey']  # Alternating grey for non-significant SNPs

    # Track chromosome positions for spacing
    chrom_offsets = {}
    cumulative_pos = 0

    for chrom in range(1, 23):
        chr_data1 = df1[df1['CHR'] == chrom]
        chr_data2 = df2[df2['CHR'] == chrom]

        if not chr_data1.empty and not chr_data2.empty:
            min_pos = min(chr_data1['P0'].min(), chr_data2['P0'].min())
            max_pos = max(chr_data1['P0'].max(), chr_data2['P0'].max())

            chr_data1 = chr_data1.copy()
            chr_data2 = chr_data2.copy()

            chr_data1.loc[:, 'plot_pos'] = chr_data1['P0'] - min_pos + cumulative_pos
            chr_data2.loc[:, 'plot_pos'] = chr_data2['P0'] - min_pos + cumulative_pos

            chrom_offsets[chrom] = cumulative_pos + (max_pos - min_pos) / 2
            cumulative_pos += max_pos - min_pos + 1e6  # Add buffer between chromosomes

            # Identify significant SNPs
            significant1 = chr_data1['-log10(P)'] >= -np.log10(5e-2)
            significant2 = chr_data2['-log10(P)'] >= -np.log10(5e-2)
            significant_annotate1 = chr_data1['-log10(P)'] >= -np.log10(1e-2)
            significant_annotate2 = chr_data2['-log10(P)'] >= -np.log10(5e-3)

            # Plot eur SNPs (bottom)
            plt.scatter(chr_data1['plot_pos'][significant1], -chr_data1['-log10(P)'][significant1],
                        color=sig_color, s=10, edgecolor='black')
            plt.scatter(chr_data1['plot_pos'][~significant1], -chr_data1['-log10(P)'][~significant1],
                        color=grey_colors[chrom % 2], s=5, alpha=0.6)
            
            # Plot pop SNPs (top)
            plt.scatter(chr_data2['plot_pos'][significant2], chr_data2['-log10(P)'][significant2],
                        color=sig_color, s=10, edgecolor='black')
            plt.scatter(chr_data2['plot_pos'][~significant2], chr_data2['-log10(P)'][~significant2],
                        color=grey_colors[chrom % 2], s=5, alpha=0.6)

            # Store all text annotations
            texts = []

            # Annotate significant SNPs with 'ID' (top)
            for i, row in chr_data2[significant_annotate2].iterrows():
                texts.append(plt.text(row['plot_pos'], row['-log10(P)'] + 0.2, row['ID'], 
                                    fontsize=7, ha='center', rotation=45))

            # Annotate significant SNPs with 'ID' (bottom)
            for i, row in chr_data1[significant_annotate1].iterrows():
                texts.append(plt.text(row['plot_pos'], -row['-log10(P)'] - 0.2, row['ID'], 
                                    fontsize=7, ha='center', rotation=45))

            # Adjust text positions to avoid overlap
            adjust_text(texts, 
            expand_text=(1.2, 1.5),  
            expand_points=(1.2, 1.5),  
            force_text=(0.3, 0.5),  
            force_points=(0.3, 0.5),  
            lim=200)  

    # Remove chromosome number blocks
    plt.xticks([])  

    # Adjust y-axis: Keep values negative but display them as positive
    max_y = max(df1['-log10(P)'].max(), df2['-log10(P)'].max()) + 1
    y_ticks = np.arange(0, np.ceil(max_y), 2)
    plt.yticks(list(-y_ticks) + list(y_ticks), labels=[str(int(abs(y))) for y in list(-y_ticks) + list(y_ticks)])
    plt.ylim(-max_y, max_y)  

    # Add genome-wide significance and 1e-3 threshold lines
    plt.axhline(y=0, color='black', linestyle='solid', linewidth=1)
    plt.axhline(y=-np.log10(5e-2), color='red', linestyle='dotted', linewidth=1)
    plt.axhline(y=np.log10(5e-2), color='red', linestyle='dotted', linewidth=1)

    # Add annotation for significance threshold
    plt.text(cumulative_pos * 1, -max_y + 0.25, "*Significance Threshold of 0.05", 
         color='red', fontsize=10, ha='right')

    # Labels and title
    plt.xlabel("Genomic Position")
    plt.ylabel("-log10(P-value)")
    plt.title(f"Miami Plot of GWAS Data\n{pop} (Top) vs eur (Bottom)")

    # Optionally save the plot
    plt.savefig(f'../plots/miami/{pop}.png', dpi=300)

def miami():
    for pop in ['afr', 'eas', 'amr', 'eur']:
        print(f'Creating Miami plot for {pop}...')
        miami_helper(pop)
        print('Done!')

def plot_go_analysis(populations):
    for pop in populations:
        print(f'Creating GO plot for {pop}...')
        file_path = os.path.join("..","data","GO_data",f"{pop}_GO_0.05.txt")
        
        df = pd.read_csv(file_path, sep='\t')
        df['log10_Pvalue'] = -np.log10(df['PValue'])
        df['term'] = df['Term'].str.split('~').str[-1]  # Extract GO term
        
        # Select top 15 terms based on P-value
        top_df = df.sort_values(by='log10_Pvalue', ascending=False).head(15)
        
        # Plot
        plt.figure(figsize=(10, 8))
        plt.barh(top_df['term'], top_df['log10_Pvalue'], color='skyblue')
        plt.xlabel('log10(Pvalue)')
        plt.ylabel('Term', fontsize=14)
        plt.title(f"Gene Ontology Analysis (pop: {pop})")
        plt.yticks(fontsize=12)
        
        # Formatting
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig(os.path.join("..","plots","GO_plots",f"{pop}_GO_0.05.png"), dpi=300)
        print('Done!')
def go_plots():
    plot_go_analysis(['afr', 'eur', 'amr', 'eas'])


