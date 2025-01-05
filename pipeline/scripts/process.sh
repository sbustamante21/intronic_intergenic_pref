#!/bin/bash

# Ensure arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <GFF_file> <RepeatMasker_file> <output_directory>"
    exit 1
fi

# Input arguments
genome=$1
rm_file=$2
output_dir=$3

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Parsear la salida del repeatmasker
tail -n +4 "$rm_file" | mawk '{print $5, $6, $7, $10, ".", $9, $11}' | sed 's/ /\t/g' | sed 's/\tC\t/\t-\t/g' > "$output_dir/rmsk_pre.bed"
mawk '{if ($7 ~ /(SINE|LINE|LTR|DNA|RC|Retroposon)/) print $0}' "$output_dir/rmsk_pre.bed" > "$output_dir/rmsk_onlyte.bed"
sed 's/\tDNA\t/\tDNA\/DNA_LTR_exNULL\t/g' "$output_dir/rmsk_onlyte.bed" > "$output_dir/rmsk_pre1.bed"
sed 's/\tLTR\t/\tLTR\/DNA_LTR_exNULL\t/g' "$output_dir/rmsk_pre1.bed" > "$output_dir/rmsk_pre2.bed"
grep -v "?" "$output_dir/rmsk_pre2.bed" > "$output_dir/rmsk_parsed.bed"

rm "$output_dir/rmsk_pre.bed" "$output_dir/rmsk_onlyte.bed" "$output_dir/rmsk_pre1.bed" "$output_dir/rmsk_pre2.bed"

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
