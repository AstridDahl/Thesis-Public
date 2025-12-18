#!/bin/bash

# Input and output file names
input_file="signature_snv_filtered_annotated4.tsv"
output_file="signatures_only.tsv"

# Extract header
header=$(head -1 "$input_file")

# Get column numbers for signature columns and Donor_ID
cols=$(echo "$header" | tr '\t' '\n' | nl -v 1 | \
       grep -E "BI_COMPOSITE_SNV_SBS|Donor_ID|Tumor_Type" | \
       awk '{printf "%s,", $1}' | sed 's/,$//')

# Use cut to extract those columns (add 1 line of header + all data)
cut -f$cols "$input_file" > "$output_file"

echo "Extracted columns saved to '$output_file'."

# 50 signatures

# Chromosome      Start_position  End_position    Strand  Variant_Classification  
# Variant_Type    Reference_Allele        maja    Tumor_Seq_Allele2       snpOverlap      
# byFrequency     sample  normid  posAndType      ref_context     VAF  
# BI_COMPOSITE_SNV_SBS1_P BI_COMPOSITE_SNV_SBS2_P BI_COMPOSITE_SNV_SBS3_P BI_COMPO
# SITE_SNV_SBS4_P BI_COMPOSITE_SNV_SBS5_P BI_COMPOSITE_SNV_SBS6_S BI_COMPOSITE_SNV
# _SBS7a_S        BI_COMPOSITE_SNV_SBS7b_S        BI_COMPOSITE_SNV_SBS7c_S        
# BI_COMPOSITE_SNV_SBS8_P BI_COMPOSITE_SNV_SBS9_P BI_COMPOSITE_SNV_SBS10a_S       
# BI_COMPOSITE_SNV_SBS11_S        BI_COMPOSITE_SNV_SBS12_P        BI_COMPOSITE_SNV
# _SBS13_P        BI_COMPOSITE_SNV_SBS14_S        BI_COMPOSITE_SNV_SBS15_S        
# BI_COMPOSITE_SNV_SBS16_P        BI_COMPOSITE_SNV_SBS17a_P       BI_COMPOSITE_SNV
# _SBS17b_P       BI_COMPOSITE_SNV_SBS18_P        BI_COMPOSITE_SNV_SBS19_P        
# BI_COMPOSITE_SNV_SBS21_S        BI_COMPOSITE_SNV_SBS22_P        BI_COMPOSITE_SNV
# _SBS26_S        BI_COMPOSITE_SNV_SBS28_P        BI_COMPOSITE_SNV_SBS30_P        
# BI_COMPOSITE_SNV_SBS33_P        BI_COMPOSITE_SNV_SBS35_P        BI_COMPOSITE_SNV
# _SBS36_P        BI_COMPOSITE_SNV_SBS37_P        BI_COMPOSITE_SNV_SBS38_S        
# BI_COMPOSITE_SNV_SBS39_P        BI_COMPOSITE_SNV_SBS40_P        BI_COMPOSITE_SNV
# _SBS60_P        BI_COMPOSITE_SNV_SBS55_S        BI_COMPOSITE_SNV_SBS44_S        
# BI_COMPOSITE_SNV_SBS61_S        BI_COMPOSITE_SNV_SBS62_S        BI_COMPOSITE_SNV
# _SBS63_S        BI_COMPOSITE_SNV_SBS64_P        BI_COMPOSITE_SNV_SBS65_S        
# BI_COMPOSITE_SNV_SBS66_S        BI_COMPOSITE_SNV_SBS67_S        BI_COMPOSITE_SNV
# _SBS68_P        BI_COMPOSITE_SNV_SBS69_P        BI_COMPOSITE_SNV_SBS70_P        
# BI_COMPOSITE_SNV_SBS71_P        BI_COMPOSITE_SNV_SBS72_P        BI_COMPOSITE_SNV
# _SBS73_S        BI_COMPOSITE_SNV_SBS74_S        BI_COMPOSITE_SNV_SBS75_S        
# BI_COMPOSITE_SNV_SBS76_S        BI_COMPOSITE_SNV_SBS77_P        BI_COMPOSITE_SNV_SBS78_S       
# BI_COMPOSITE_SNV_SBS79_S        BI_COMPOSITE_SNV_SBS80_P        BI_COMPOSITE_SNV_SBS81_P        BI_COMPOSITE_SNV_SBS82_P        BI_COMPOSITE_SNV_SBS83_P     
# n_pc    eid     kmer_11 alt_allele      Donor_ID        Tumor_Type      colcode project_code.
