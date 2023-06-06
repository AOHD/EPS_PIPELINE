#!/bin/bash

# Input file name
input_file="/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/concatenatedRSEM_nofasta.gff"

# Output file name
output_file="/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/concatenatedRSEM_expanded.gff"

gene_counter=1

while IFS= read -r line; do
  # Original line
  name=";Name=$gene_counter"
  modified_line="${line/CDS/gene}$name"
  echo "$modified_line"

  # First expanded line
  IFS=';' read -ra fields <<< "$line"
  id=${fields[0]#*=}
  parent=";Parent=$id"
  expanded_line="${line/CDS/mRNA}"
  expanded_line="${expanded_line/ID=$id;/ID=$id\_mRNA;}"
  expanded_line="$expanded_line$parent"
  echo "$expanded_line"

  # Second expanded line
  IFS=';' read -ra fields <<< "$line"
  id_mrna=${fields[0]#*=}
  parent_mrna=";Parent=${id_mrna}_mRNA"
  expanded_line2="${line/CDS/exon}"
  expanded_line2="${expanded_line2/ID=$id;/ID=$id\_exon;}"
  expanded_line2="$expanded_line2$parent_mrna"
  echo "$expanded_line2"

  ((gene_counter++))

done < "$input_file" > "$output_file"
