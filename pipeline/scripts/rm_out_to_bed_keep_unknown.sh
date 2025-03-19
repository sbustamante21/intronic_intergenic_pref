# Ensure arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <GFF_file> <output_file>"
    exit 1
fi

rm_file=$1
output_file=$2

# Parsear la salida del repeatmasker
tail -n +4 "$rm_file" | mawk '{print $5, $6, $7, $10, ".", $9, $11}' | sed 's/ /\t/g' | sed 's/\tC\t/\t-\t/g' > "rmsk_pre.bed"
mawk '{if ($7 ~ /(SINE|LINE|LTR|DNA|RC|Retroposon|Unknown)/) print $0}' "rmsk_pre.bed" > "rmsk_onlyte.bed"
sed 's/\tDNA\t/\tDNA\/DNA_LTR_exNULL\t/g' "rmsk_onlyte.bed" > "rmsk_pre1.bed"
sed 's/\tLTR\t/\tLTR\/DNA_LTR_exNULL\t/g' "rmsk_pre1.bed" > "rmsk_pre2.bed"
grep -v "?" "rmsk_pre2.bed" > "$output_file"

rm "rmsk_pre.bed" "rmsk_onlyte.bed" "rmsk_pre1.bed" "rmsk_pre2.bed"