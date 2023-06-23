##All of these commands are to be run before the psiblast.sh script can be executed
##They ensure that the files needed to run PSI-BLAST (and gtdbtk de novo tree) are in the correct folders and are formated correctly.

#Moves all .faa files from their Prokka folders to the faa folder
for file in /user_data/ahd/EPS_PIPELINE/data/raw/prokka/*/*.faa; do cp $file /user_data/ahd/EPS_PIPELINE/data/raw/faa; echo $file; done

#Moves all .gff files from their Prokka folders to the gff folder
for file in /user_data/ahd/EPS_PIPELINE/data/raw/prokka/*/*.gff; do cp $file /user_data/ahd/EPS_PIPELINE/data/raw/gff; echo $file; done

ln -s /user_data/ahd/EPS_PIPELINE/data/raw/prokka/*/*.gff /user_data/ahd/EPS_PIPELINE/data/raw/gff

##Removes FASTA from gff for each numbered file and moves them to reduced_gff
##Run this in the gff folder
for file in *gff*; do sed '1,/##FASTA/!d' $file | sed 's/##FASTA.*//' > "/user_data/ahd/EPS_PIPELINE/data/raw/reduced_gff/$file"; echo $file; done

for file in *gff*; do sed '1,/##FASTA/!d' $file | sed 's/##FASTA.*//' > "/user_data/ahd/EPS_PIPELINE/isolate_genome/$file"; echo $file; done
##Now generate_gff.R can be run

##Makeblastdb
cat /user_data/ahd/EPS_PIPELINE/data/raw/faa/*.faa > /user_data/ahd/EPS_PIPELINE/databases/mydb.fas
module load BLAST+/2.12.0-gompi-2020b
makeblastdb -in /user_data/ahd/EPS_PIPELINE/databases/mydb.fas -dbtype prot -out /user_data/ahd/EPS_PIPELINE/databases/mydb

makeblastdb -in /user_data/ahd/EPS_PIPELINE/isolates/zooglea/prokka/MGIOFAJM.faa -dbtype prot -out /user_data/ahd/EPS_PIPELINE/isolates/zooglea/database/isolatedb


##Creating HQ-MAG symlinks from the symlink container
for filename in $(cat /user_data/ahd/EPS_PIPELINE/text.txt); do cp /projects/glomicave/HQMAG_database/symlinks/$filename /user_data/ahd/EPS_PIPELINE/data/raw/bins; echo $filename; done

#After psiblast is run, proximity filtration is run, and then IPS is run on the results of that.

## Metatranscriptomics ##

#######
#detectM
#######
for file in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/detectm/*/detectm_tpm_scores.tsv;
do cut -f 1,11,14  $file > /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/detectm_geneTPMannot/$(basename "$(dirname "$file")").tsv;
echo $(basename "$(dirname "$file")").tsv;
done


#######
#RSEM
#######

#Concatenates all .ffn (nucleotide) files from their Prokka folders to the metatranscriptomics folder
for file in /user_data/ahd/EPS_PIPELINE/data/raw/prokka/*/*.ffn; do cat $file >> /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/combined.ffn; echo $file; done

#Creates RSEM_TPM_summarised.tsv file
cut -f 1 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/results/LIB-Glomicave-0001_resequenced_mRNA.genes.results > \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv

#For loop that takes .genes.results files from /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM/results/, cuts the gene_id and TPM columns, renames the TPM column to the file name, and joins them all into one file by gene_id
for file in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/results/*.genes.results; \
do cut -f 1,6 $file | sed 's/TPM/'$(basename "$file")'/' | join - /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv > /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised_temp.tsv; \
mv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised_temp.tsv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv; \
echo $file; \
done


#Creates RSEM_TPM_summarised.tsv file
cut -f 1 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/results/LIB-Glomicave-0001_resequenced_mRNA.genes.results > \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.tsv

#For loop that takes .genes.results files from /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM/results/, cuts the gene_id and counts columns, renames the counts column to the file name, and joins them all into one file by gene_id
for file in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/results/*.genes.results; \
do cut -f 1,5 $file | sed 's/expected_count/'$(basename "$file")'/' | join - /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.tsv > /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised_temp.tsv; \
mv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised_temp.tsv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.tsv; \
echo $file; \
done


#Creates RSEM_TPM_MAG_summarised.tsv file
cut -f 1 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/LIB-Glomicave-0001_resequenced_mRNA_filtered.bam.genes.results.tsv > \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_MAG_summarised.tsv

#Removes "_filtered.bam.genes.results.tsv" from the end of the file names in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files
for file in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/*; \
do mv $file $(echo $file | sed 's/_filtered.bam.genes.results.tsv//'); \
echo $file; \
done

#For loop that takes .genes.results files from /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/, cuts the gene_id and MAG_TPM columns, renames the MAG_TPM column to the file name, and joins them all into one file by gene_id
for file in /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/*; \
do cut -f 1,10 $file | sed 's/TPM/'$(basename "$file")'/' | join - /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_MAG_summarised.tsv > /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_MAG_summarised_temp.tsv; \
mv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_MAG_summarised_temp.tsv /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_MAG_summarised.tsv; \
echo $file; \
done