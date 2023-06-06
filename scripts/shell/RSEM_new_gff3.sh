#Preparing the reference for RSEM
#I removed fasta sequences from the concatenated gff file /projects/glomicave/metatranscriptomics/04_dirseq/concatenated.gff
#And then i modified it using /user_data/ahd/EPS_PIPELINE/scripts/shell/gff3_modifier.sh to make it work with rsem-prepare-reference
rsem-prepare-reference --gff3 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/concatenatedRSEM_expanded.gff --bowtie2 /projects/glomicave/HQMAG_database/concatenated_assembly/HQ_sp_692_concat.fa /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/reference/combined.ORFs

#Running rsem-calculate-expression with 10 threads parallel
ls /projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/ | grep _R1.fq | sed 's/_R1.fq//' | \
parallel -j10 rsem-calculate-expression --paired-end --bowtie2 --bowtie2-sensitivity-level sensitive -p 10 \
/projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/{}_R1.fq \
/projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/{}_R2.fq \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/reference/combined.ORFs \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new_gff3/results/{}