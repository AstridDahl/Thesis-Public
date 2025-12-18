#!/bin/bash

# Input and output file names
input_file="signature_snv_filtered_annotated4.tsv"
output_file="kmer5_signature_snv_filtered_annotated4.tsv"

# Use awk to extract relevant columns and modify kmer_5 to include the alternative allele
awk 'BEGIN { FS=OFS="\t" }
     NR==1 { print "kmer_5", "alt_allele", "Donor_ID", "Tumor_Type"; next }  # Updated header with kmer_5
     { ref_allele=substr($79,6,1); alt_allele=$80; kmer_left=substr($79,2,2); kmer_right=substr($79,7,2); 
       print kmer_left "[" ref_allele ">" alt_allele "]" kmer_right, $80, $81, $82
       }' "$input_file" > "$output_file"

echo "Processing complete. Modified file saved as '$output_file'."
