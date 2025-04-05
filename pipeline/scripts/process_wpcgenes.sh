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

# Sacar los genes en general
mawk '$3 == "gene"' "$genome" > "$output_dir/genes_all.gff"
mawk '!/^#/{print $1, $4, $5, ".", ".", $7, ".", $3, ".", $9}' "$output_dir/genes_all.gff" | sed 's/ /\t/g' > "$output_dir/genes_all.bed"

# Sacar todos los protein coding genes
mawk '$3 == "gene"' "$genome" | grep "protein_coding" > "$output_dir/genes.gff"
mawk '!/^#/{print $1, $4, $5, ".", ".", $7, ".", $3, ".", $9}' "$output_dir/genes.gff" | sed 's/ /\t/g' > "$output_dir/genes_pc.bed"

# Sacar todos los exones y pasarlos a bed
mawk '$3 == "exon"' "$genome" > "$output_dir/exons.gff"
mawk '{print $1, $4, $5, $3, ".", $7}' "$output_dir/exons.gff" | sed 's/ /\t/g' > "$output_dir/exons_pre1.bed"

# Intersectar exones con protein-coding genes (aqui se sacan exones de genes protein-coding)
bedtools intersect -a "$output_dir/exons_pre1.bed" -b "$output_dir/genes_pc.bed" > "$output_dir/exons_pre2.bed"

# sort y merge a los exones (de genes protein_coding)
sort -k1,1 -k2,2n "$output_dir/exons_pre2.bed" > "$output_dir/exons_sorted.bed"
bedtools merge -s -i "$output_dir/exons_sorted.bed" > "$output_dir/exons.bed"

# sort y merge a los genes (protein coding)
sort -k1,1 -k2,2n "$output_dir/genes_pc.bed" > "$output_dir/genes_pc_sorted.bed"
bedtools merge -s -i "$output_dir/genes_pc_sorted.bed" > "$output_dir/genes.bed"
rm "$output_dir/genes_pc.bed" "$output_dir/genes_pc_sorted.bed"

# merge a los genes no protein coding
sort -k1,1 -k2,2n "$output_dir/genes_all.bed" > "$output_dir/genes_all_sorted.bed"
bedtools merge -s -i "$output_dir/genes_all_sorted.bed" > "$output_dir/genes_all_merged.bed"

#rm "$output_dir/genes_all.bed" "$output_dir/genes_all_sorted.bed"

# Sacar intrones de genes en total
sort -k1,1 -k2,2n "$output_dir/exons_pre1.bed" > "$output_dir/exons_sorted_all.bed"
bedtools merge -s -i "$output_dir/exons_sorted_all.bed" > "$output_dir/exons_all.bed" # sort y merge a exones de genes en total
bedtools subtract -a "$output_dir/genes_all_merged.bed" -b "$output_dir/exons_all.bed" > "$output_dir/introns_all.bed"

rm "$output_dir/exons.gff" "$output_dir/exons_pre1.bed" "$output_dir/exons_pre2.bed" "$output_dir/exons_sorted.bed"

# Sacar intrones (de protein_coding genes)
bedtools subtract -a "$output_dir/genes.bed" -b "$output_dir/exons.bed" > "$output_dir/introns.bed"

# TE por intron (de protein_coding genes)
bedtools intersect -f 0.75 -wa -wb -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/introns.bed" > "$output_dir/em_per_intron.bed"

# También sacar TE por gen (de protein_coding gene)
bedtools intersect -f 1 -wa -wb -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/genes.bed" > "$output_dir/em_per_gene.bed"

# Saber cuantos elementos moviles hay en zonas intergenicas, nivel de genoma completo
bedtools intersect -v -a "$output_dir/rmsk_parsed.bed" -b "$output_dir/genes_all_merged.bed" > "$output_dir/em_intergenic.bed"

cp "$chr_lengths" "$output_dir/chromosome_lengths.txt"
