import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import multiprocessing as mp
from scipy import stats
import argparse
import seaborn as sns

def feat_len(dataset):
    return [row[2] - row[1] for _, row in dataset.iterrows()]

def count_genes_per_chr(chr_name, df, chr_col_index):
    return chr_name, df[df[chr_col_index] == chr_name].shape[0]

def parallel_gene_count(chr_lengths, df, chr_col_index=0):
    with mp.Pool(mp.cpu_count()) as pool:
        # Use pool.starmap to pass multiple arguments to the function
        results = pool.starmap(count_genes_per_chr, [(chr_name, df, chr_col_index) for chr_name in chr_lengths.keys()])
    
    # Convert the results into a dictionary
    return dict(results)

def occ_report(dataset1, intr_len, length_intergenic, intergenic):
    # dataset1 tiene el arreglo con todos los arreglos de los largos de los elementos moviles
    # intr_len es el largo de todos los intrones codificantes a proteina
    # length_intergenic es el largo de las zonas intergenicas
    # intergenic es el arreglo con arreglo de los tamaños de los elementos moviles en zonas intergenicas

    all_lens = 0
    all_intergenics = 0
    titles = ["SINE", "LINE", "LTR", "DNA", "RC", "Retroposon"]
    intergenics = []
    intronics = []

    
    for length1, title, length2 in zip(dataset1, titles, intergenic):
        all_lens += sum(length1)
        all_intergenics += sum(length2)

        occ_intr = sum(length1)/intr_len
        occ_intergenic = sum(length2)/length_intergenic

        #print(f"Para {title}")
        #print(f"Ocupancia en intrones : {occ_intr}")
        #print(f"Ocupancia intergenica: {occ_intergenic}")

        intronics.append(occ_intr)
        intergenics.append(occ_intergenic)

    #print(f"En total en intrones: {all_lens/intr_len}")
    #print(f"En total en intergenica: {all_intergenics/length_intergenic}")

    return(intronics, intergenics)

def t_test(data1, data2, titles):
    # Initialize lists to store results
    results = []

    # Iteratively perform the t-tests
    for length1, length2, title in zip(data1, data2, titles):
        t_stat, p_value = stats.ttest_ind(length1, length2)

        # Collect the data for this test
        result = {
            "Title": title,
            "Mean 1": np.mean(length1),
            "Mean 2": np.mean(length2),
            "T-Statistic": t_stat,
            "P-Value": p_value,
            "Significant": p_value < 0.05
        }
        results.append(result)

    # Convert results to a dataframe
    df_results = pd.DataFrame(results)
    return df_results

# A partir del em_per_x hacer un grafico del tamaño de todos los tipos de TE
def process_tes(dataset):
    # Extract lengths and types
    lengths = [(row[2] - row[1] + 1, row[6].split("/")[0]) for _, row in dataset.iterrows()]

    # Initialize lists for each TE type
    sine_lengths = []
    line_lengths = []
    ltr_lengths = []
    dna_lengths = []
    rc_lengths = []
    retroposon_lengths = []

    # Classify lengths by TE type
    for length, te_type in lengths:
        if te_type == "SINE":
            sine_lengths.append(length)
        elif te_type == "LINE":
            line_lengths.append(length)
        elif te_type == "LTR":
            ltr_lengths.append(length)
        elif te_type == "DNA":
            dna_lengths.append(length)
        elif te_type == "RC":
            rc_lengths.append(length)
        elif te_type == "Retroposon":
            retroposon_lengths.append(length)
    
    return [sine_lengths, line_lengths, ltr_lengths, dna_lengths, rc_lengths, retroposon_lengths]

def calculate_class_data(data):
    class_lengths = []
    class_freqs = []

    for lengths in data:
        n_elements = len(lengths)
        mean_length = sum(lengths) / n_elements if lengths else 0

        class_lengths.append(mean_length)
        class_freqs.append(n_elements)

    return class_lengths, class_freqs

# Tomar las ocupancias en intrones y zonas intergenicas, calcular la preferencia dividiendo
def pref_intron_inter(occs_intronic, occs_intergenic):
    prefs = []
    for intronic, intergenic in zip(occs_intronic, occs_intergenic):
        try:
            prefs.append(intronic/intergenic)
        except:
            prefs.append(None)
    print(prefs)
    return prefs

def plot_histograms_overlapping(data1, data2, supertitle, stats1, stats2, output_file):
    # Unpack statistics
    class_lengths1, class_freqs1 = stats1
    class_lengths2, class_freqs2 = stats2

    # Create a 2x3 grid of subplots
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()

    # Data and titles for subplots
    titles = ["SINE", "LINE", "LTR", "DNA", "RC", "Retroposon"]

    # Plot each histogram
    for ax, lengths1, lengths2, title, mean_length1, n_elements1, mean_length2, n_elements2 in zip(
        axs, data1, data2, titles, class_lengths1, class_freqs1, class_lengths2, class_freqs2
    ):
        ax.hist(lengths1, bins=30, alpha=0.5, edgecolor='black', label='Data1')  # Histogram for data1
        ax.hist(lengths2, bins=30, alpha=0.5, edgecolor='black', label='Data2')  # Histogram for data2

        ax.set_title(f'Lengths of {title} Elements')
        ax.set_xlabel('Length')
        ax.set_ylabel('Frequency')
        ax.legend([
            f'Intronic - Mean: {mean_length1:.2f}, N: {n_elements1}', 
            f'Intergenic - Mean: {mean_length2:.2f}, N: {n_elements2}'
        ])

    plt.suptitle(supertitle, fontsize=16)
    plt.tight_layout()
    
    # Save the plot to the specified file
    plt.savefig(output_file)

    # Optionally, close the plot to free up memory
    plt.close()

def occ_report_chr_by_chr(
    chr_lengths, introns_df, intergenic_df, fullintron_df, process_tes_func, feat_len_func, occ_report_func
):
    """
    Calculate occurrences and frequencies of intronic and intergenic features per chromosome.

    Parameters:
    - chr_lengths: dict, lengths of chromosomes {chr_name: length}.
    - introns_df: DataFrame, intron data.
    - intergenic_df: DataFrame, intergenic region data.
    - fullintron_df: DataFrame, full intron length data.
    - process_tes_func: callable, function to process transposable elements.
    - feat_len_func: callable, function to calculate feature lengths.
    - occ_report_func: callable, function to generate occurrence reports.

    Returns:
    - occs_per_chr_intronic: dict, occurrences per chromosome for intronic regions.
    - occs_per_chr_intergenic: dict, occurrences per chromosome for intergenic regions.
    - freqs_per_chr_intronic: dict, frequencies per chromosome for intronic regions.
    - freqs_per_chr_intergenic: dict, frequencies per chromosome for intergenic regions.
    """
    occs_per_chr_intronic = {}
    occs_per_chr_intergenic = {}
    freqs_per_chr_intronic = {}
    freqs_per_chr_intergenic = {}

    for chr_name, chr_length in chr_lengths.items():
        #print(f"Processing {chr_name}")

        df_chr_intron = process_tes_func(introns_df[introns_df[0] == chr_name])
        df_chr_intergenic = process_tes_func(intergenic_df[intergenic_df[0] == chr_name])
        df_chr_fullintron = sum(feat_len_func(fullintron_df[fullintron_df[0] == chr_name]))

        if df_chr_fullintron == 0:
            continue

        mini_intergenic_size = chr_length - df_chr_fullintron

        occs_per_chr_intronic[chr_name], occs_per_chr_intergenic[chr_name] = occ_report_func(
            df_chr_intron, df_chr_fullintron, mini_intergenic_size, df_chr_intergenic
        )

        freqs_per_chr_intronic[chr_name] = [len(el) for el in df_chr_intron]
        freqs_per_chr_intergenic[chr_name] = [len(el) for el in df_chr_intergenic]

        #print("\n")

    return occs_per_chr_intronic, occs_per_chr_intergenic, freqs_per_chr_intronic, freqs_per_chr_intergenic

def make_pref_df(df1, df2, key_column):
    """
    Divide numeric columns of two dataframes row-wise based on a key column.
    
    Parameters:
        df1 (pd.DataFrame): First dataframe with numeric values to be divided.
        df2 (pd.DataFrame): Second dataframe with numeric values to divide by.
        key_column (str): The column to align rows on (e.g., "Chromosome").
    
    Returns:
        pd.DataFrame: A new dataframe with the key column and the division results.
    """
    # Merge dataframes on the key column
    merged_df = pd.merge(df1, df2, on=key_column, suffixes=("_1", "_2"))

    # Identify numeric columns for both dataframes
    numeric_cols_1 = [col for col in df1.columns if col != key_column]
    numeric_cols_2 = [col for col in df2.columns if col != key_column]

    # Initialize result dataframe with the key column
    result_df = merged_df[[key_column]].copy()

    # Perform element-wise division for matching numeric columns
    for col1, col2 in zip(numeric_cols_1, numeric_cols_2):
        result_df[col1] = merged_df[f"{col1}_1"] / merged_df[f"{col2}_2"]

    return result_df

def pref_df_to_csv(df, filename, input_dir):
    # Reset index to make 'Chromosome' a column
    df.reset_index(inplace=True)
    df.rename(columns={"index": "Chromosome"}, inplace=True)
    # Save DataFrame to a CSV file    
    df.to_csv(f"{input_dir}/{filename}.csv", sep="\t", index=False)

def save_heatmap(df, key_column, filename="heatmap.png", title="Heatmap of Division Results", figsize=(10, 8)):
    """
    Create and save a heatmap from the resulting dataframe where each row corresponds to a chromosome.
    
    Parameters:
        df (pd.DataFrame): The dataframe containing the data.
        key_column (str): The column used as row labels (e.g., "Chromosome").
        filename (str): The file name where the heatmap will be saved.
        title (str): The title of the heatmap.
        figsize (tuple): The size of the figure.
    """
    # Set the index to the key column
    heatmap_data = df.set_index(key_column)

    # Plot the heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(
        heatmap_data,
        annot=True,  # Show values in each cell
        cmap="coolwarm",  # Choose a color palette
        fmt=".2f",  # Format for the numbers
        linewidths=0.5,  # Line width between cells
        cbar_kws={"label": "Value"},  # Label for the color bar
        center=1,  # Center the color palette at 1
    )

    # Add title and labels
    plt.title(title, fontsize=16)
    plt.xlabel("Features", fontsize=12)
    plt.ylabel("Chromosome", fontsize=12)

    # Save the plot
    plt.tight_layout()
    plt.savefig(filename, dpi=300)  # Save with high resolution
    plt.close()  # Close the plot to free up memory


def main(input_dir, organism_name):
    titles = ["SINE", "LINE", "LTR", "DNA", "RC", "Retroposon"]

    pcgenes_df = pd.read_csv(f"{input_dir}/em_per_gene.bed", sep="\t", header=None)

    introns_df = pd.read_csv(f"{input_dir}/em_per_intron.bed", sep="\t", header=None)
    intron_lengths = process_tes(introns_df)

    intergenic_df = pd.read_csv(f"{input_dir}/em_intergenic.bed", sep="\t", header=None)
    intergenic_lengths = process_tes(intergenic_df)

    only_intron_df = pd.read_csv(f"{input_dir}/introns.bed", sep="\t", header=None)
    only_intron_lengths = sum(feat_len(only_intron_df))

    only_gene_df = pd.read_csv(f"{input_dir}/genes.bed", sep="\t", header=None)
    only_gene_lengths = sum(feat_len(only_gene_df))

    chr_lengths_df = pd.read_csv(f"{input_dir}/chromosome_lengths.txt", sep="\t", header=None, names=["chr", "length"])
    chr_lengths = chr_lengths_df.set_index("chr").to_dict()["length"]

    genes_per_chr = parallel_gene_count(chr_lengths, pcgenes_df, chr_col_index=0)

    genome_size = sum(chr_lengths.values())
    genome_intergenic_size = genome_size - only_intron_lengths

    # Histogramas intronicos e intergenicos
    intron_lengths_graph, intron_freqs_graph = calculate_class_data(intron_lengths)
    intergenic_lengths_graph, intergenic_freqs_graph = calculate_class_data(intergenic_lengths)

    plot_histograms_overlapping(intron_lengths, intergenic_lengths, f"{organism_name}, intronic vs intergenic", (intron_lengths_graph, intron_freqs_graph), (intergenic_lengths_graph, intergenic_freqs_graph), f"{input_dir}/histograms_overlapping.png")

    # Length report
    length_df = pd.DataFrame(data=[intron_lengths_graph, intergenic_lengths_graph], index=["Intronic", "Intergenic"], columns=titles)
    preference = length_df.loc["Intronic"] / length_df.loc["Intergenic"]
    # Add the preference row to the DataFrame
    length_df.loc["Preference"] = preference
    length_df.to_csv(f"{input_dir}/length_report.tsv", sep="\t")

    # Frequency report
    freq_df = pd.DataFrame(data=[intron_freqs_graph, intergenic_freqs_graph], index=["Intronic", "Intergenic"], columns=titles)
    preference = freq_df.loc["Intronic"] / freq_df.loc["Intergenic"]
    # Add the preference row to the DataFrame
    freq_df.loc["Preference"] = preference
    freq_df.to_csv(f"{input_dir}/freq_report.tsv", sep="\t")

    # Occuapncy report
    # ocupancias intronicas de las clases de TE segun el orden de titles, y las ocupancias intergenicas de las clases de TE segun el orden de titles
    # La sumatoria de cada lista es el total de ocupancia para cada zona
    occ_intronic, occ_intergenic = occ_report(intron_lengths, only_intron_lengths, genome_intergenic_size, intergenic_lengths)
    occ_df = pd.DataFrame(data=[occ_intronic, occ_intergenic], index=["Intronic", "Intergenic"], columns=titles)
    preference = occ_df.loc["Intronic"] / occ_df.loc["Intergenic"]
    # Add the preference row to the DataFrame
    occ_df.loc["Preference"] = preference
    occ_df.to_csv(f"{input_dir}/occ_report.tsv", sep="\t")

    # t-test
    test_df = t_test(intron_lengths, intergenic_lengths, titles)
    test_df.to_csv(f"{input_dir}/t_test.tsv", sep="\t")

    # Chromosome by chromosome processes
    occs_per_chr_intronic, occs_per_chr_intergenic, freqs_per_chr_intronic, freqs_per_chr_intergenic = occ_report_chr_by_chr(chr_lengths, introns_df, intergenic_df, only_intron_df, process_tes, feat_len, occ_report)

    # Convert dictionary to DataFrame
    df_occs_per_chr_intronic = pd.DataFrame.from_dict(occs_per_chr_intronic, orient="index", columns=titles)
    pref_df_to_csv(df_occs_per_chr_intronic, "occ_per_chr_intronic", input_dir)

    df_occs_per_chr_intergenic = pd.DataFrame.from_dict(occs_per_chr_intergenic, orient="index", columns=titles)
    pref_df_to_csv(df_occs_per_chr_intergenic, "occ_per_chr_intergenic", input_dir)

    df_freqs_per_chr_intronic = pd.DataFrame.from_dict(freqs_per_chr_intronic, orient="index", columns=titles)
    pref_df_to_csv(df_freqs_per_chr_intronic, "freq_per_chr_intronic", input_dir)

    df_freqs_per_chr_intergenic = pd.DataFrame.from_dict(freqs_per_chr_intergenic, orient="index", columns=titles)
    pref_df_to_csv(df_freqs_per_chr_intergenic, "freq_per_chr_intergenic", input_dir)

    occs_pref_chr = make_pref_df(df_occs_per_chr_intronic, df_occs_per_chr_intergenic, "Chromosome")
    pref_df_to_csv(occs_pref_chr, "occ_per_chr_intronic_over_intergenic", input_dir)

    save_heatmap(occs_pref_chr, "Chromosome", filename=f"{input_dir}/heatmap_occs_per_chr_intronic_over_intergenic.png", title=f"Heatmap of Occupancy Preference by Chromosome in {organism_name}")

    freqs_pref_chr = make_pref_df(df_freqs_per_chr_intronic, df_freqs_per_chr_intergenic, "Chromosome")
    pref_df_to_csv(freqs_pref_chr, "freq_per_chr_intronic_over_intergenic", input_dir)

    save_heatmap(freqs_pref_chr, "Chromosome", filename=f"{input_dir}/heatmap_freqs_per_chr_intronic_over_intergenic.png", title=f"Heatmap of Frequency Preference by Chromosome in {organism_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process genomic data")
    parser.add_argument("input_dir", help="Base directory containing the input files")
    parser.add_argument("organism_name", help="Name of the organism")
    args = parser.parse_args()
    
    main(args.input_dir, args.organism_name)