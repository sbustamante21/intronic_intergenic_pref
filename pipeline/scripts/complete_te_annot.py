import pandas as pd
import argparse
import sys
import os


def main():
    # Receive the original dataframe
    parser = argparse.ArgumentParser(description="Clean the RepeatModeler's output file")
    parser.add_argument("input", help="Path to the input file - intersecion of RepeatModeler and RepeatMasker")
    parser.add_argument("input2", help="Path to the input file - RepeatModeler's unique features")
    parser.add_argument("input3", help="Path to the input file - RepeatMasker's unique features")
    parser.add_argument("output", help="Path to the output file")
    args = parser.parse_args()
    
    df = pd.read_csv(args.input, sep="\t", header=None)
    unique_rmodel = pd.read_csv(args.input2, sep="\t", header=None)
    unique_gcf = pd.read_csv(args.input3, sep="\t", header=None)

    # If column 7 is "Unknown", replace with column 14
    df[6] = df.apply(lambda row: row[13] if row[6] == "Unknown" else row[6], axis=1)

    # Drop column 8 to 14
    df.drop(df.columns[7:15], axis=1, inplace=True)
    df.drop_duplicates(inplace=True)

    # Delete rows 6 that have value "Unknown" from unique_rmodel
    unique_rmodel = unique_rmodel[unique_rmodel[6] != "Unknown"]

    # Add rows from uniques to df
    df = pd.concat([df, unique_rmodel, unique_gcf], ignore_index=True)
    df.drop_duplicates(inplace=True)

    # Save the cleaned dataframe to a new file
    df.to_csv(args.output, sep="\t", index=False, header=None)

if __name__ == "__main__":
    main()