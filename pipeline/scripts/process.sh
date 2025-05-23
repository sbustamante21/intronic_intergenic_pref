#!/bin/bash

# Ensure arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <GFF_file> <RepeatMasker BED file> <output_directory>  <Chromosome Lengths file>"
    exit 1
fi

# Input arguments
genome=$1
rm_file=$2
output_dir=$3
chr_lengths=$4

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Parsear la salida del repeatmasker
cp "$rm_file" "$output_dir/rmsk_parsed.bed"

# Sacar todos los genes
mawk '$3 == "gene"' "$genome" > "$output_dir/genes.gff"
mawk '!/^#/{print $1, $4, $5, ".", ".", $7, ".", $3, ".", $9}' "$output_dir/genes.gff" | sed 's/ /\t/g' > "$output_dir/genes.bed"

# Sacar todos los exones y pasarlos a bed
mawk '$3 == "exon"' "$genome" > "$output_dir/exons.gff"
mawk '{print $1, $4, $5, $3, ".", $7}' "$output_dir/exons.gff" | sed 's/ /\t/g' > "$output_dir/exons_pre1.bed"

# sort y merge a los exones
sort -k1,1 -k2,2n "$output_dir/exons_pre1.bed" > "$output_dir/exons_sorted.bed"
bedtools merge -s -i "$output_dir/exons_sorted.bed" > "$output_dir/exons.bed"

rm "$output_dir/exons.gff" "$output_dir/exons_pre1.bed" "$output_dir/exons_sorted.bed"

# Sacar intrones
bedtools subtract -a "$output_dir/genes.bed" -b "$output_dir/exons.bed" > "$output_dir/introns.bed"

# TE por intron
bedtools intersect -f 0.75 -wa -wb -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/introns.bed" > "$output_dir/em_per_intron.bed"

# También sacar TE por gen
bedtools intersect -f 1 -wa -wb -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/genes.bed" > "$output_dir/em_per_gene.bed"

# Saber cuantos elementos moviles hay en zonas intergenicas, nivel de genoma completo
bedtools intersect -v -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/genes.bed" > "$output_dir/em_intergenic.bed"

cp "$chr_lengths" "$output_dir/chromosome_lengths.txt"
