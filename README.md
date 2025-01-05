# Usage
Have a directory called genomes, with 1 subdirectory per organism.
Each organism must have a GFF, a RM OUT file and a chromosome lengths file. Two columns, one for chromosomes and one for lengths. No header.
Run scripts/process.sh
copy chromosome lengths to results/organism
run process_tes.py
run multi_processing to compare all organisms
