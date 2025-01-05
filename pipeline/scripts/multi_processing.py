import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def main():

    occ_reports = {}
    length_reports = {}
    freq_reports = {}

    # iterate over the organism subdirectories
    for subdir in os.listdir(args.base_dir):
        if os.path.isdir(os.path.join(args.base_dir, subdir)):
            # iterate over the vectors files
            for vector in os.listdir(os.path.join(args.base_dir, subdir)):
                if vector == "occ_report.tsv":
                    # check the subdirectory name to name the dataframe
                    df = pd.read_csv(os.path.join(args.base_dir, subdir, vector), sep="\t")
                    occ_reports[subdir] = df
                elif vector == "length_report.tsv":
                    df = pd.read_csv(os.path.join(args.base_dir, subdir, vector), sep="\t")
                    length_reports[subdir] = df
                elif vector == "freq_report.tsv":
                    df = pd.read_csv(os.path.join(args.base_dir, subdir, vector), sep="\t")
                    freq_reports[subdir] = df

    organisms = []
    occ_preferences = []
    length_preferences = []
    freq_preferences = []

    for organism in occ_reports.keys():
        organisms.append(organism)
        occ_preferences.append(occ_reports[organism].loc[2])
        length_preferences.append(length_reports[organism].loc[2])
        freq_preferences.append(freq_reports[organism].loc[2])

    os.makedirs(f"{args.base_dir}/comparisons", exist_ok=True)

    df_occ = pd.DataFrame(occ_preferences, index=organisms)
    df_occ = df_occ.drop(columns=["Unnamed: 0"])
    df_occ.to_csv(f"{args.base_dir}/comparisons/occ_preferences.tsv", sep="\t")

    df_length = pd.DataFrame(length_preferences, index=organisms)
    df_length = df_length.drop(columns=["Unnamed: 0"])
    df_length.to_csv(f"{args.base_dir}/comparisons/length_preferences.tsv", sep="\t")

    df_freq = pd.DataFrame(freq_preferences, index=organisms)
    df_freq = df_freq.drop(columns=["Unnamed: 0"])
    df_freq.to_csv(f"{args.base_dir}/comparisons/freq_preferences.tsv", sep="\t")

    # plot heatmaps for each dataframe
    i = 1
    for df in [df_occ, df_length, df_freq]:
        plt.figure(figsize=(10,9))
        sns.heatmap(df, cmap="coolwarm", annot=True, fmt=".2f", center=1)
        plt.savefig(f"{args.base_dir}/comparisons/{i}.png")
        plt.close()
        i = i + 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparisons between multiple organisms")
    parser.add_argument("base_dir", help="Base directory containing organism subdirectories.")
    args = parser.parse_args()
    main()