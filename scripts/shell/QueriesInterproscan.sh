#!/usr/bin/env bash
# Working directionary, our group folder
WD=/user_data/ahd/EPS_PIPELINE
cd $WD
date=$(date '+%d_%b')
# Folder with all subsetted fasta files
FASTA=/user_data/ahd/EPS_PIPELINE/data/raw/queries_psiblast
# Results folder
mkdir $WD/interproscan_results/ips_queries/$date
RESULTS=$WD/interproscan_results/ips_queries/$date

# Loading interproscan
module load InterProScan/5.38-76.0-foss-2018a

for i in $FASTA/*EPS273*; do
i=`basename $i | rev | cut -c 5- | rev`
# Removing .faa from the name
echo $i
# Running interproscan
## Only Pfam and SUPERFAMILY database
## Only GFF3 file is written out
interproscan.sh -cpu 50 -appl Pfam -f GFF3 -o $RESULTS/$i -i $FASTA/$i.faa

# Remove trailing fasta sequence in gff3 file and leading three lines
F=$RESULTS/$i
sed '/^##FASTA$/,$d' ${F} | sed '1,/^##interproscan*/d' > $RESULTS/$i.gff3
done
module purge

#interproscan.sh -cpu 15 -appl Pfam -f GFF3 -o /user_data/ahd/EPS_PIPELINE/interproscan_results/ips_queries/levan -i /user_data/ahd/EPS_PIPELINE/data/raw/queries_psiblast/levan.faa