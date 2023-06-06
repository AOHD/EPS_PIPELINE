#Preparing the reference for RSEM - without Zivile and Caitlin's extra rRNA filtering
rsem-prepare-reference --bowtie2 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/combined.ffn /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/combined.ORFs

#Running rsem-calculate-expression with 10 threads parallel
ls /projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/ | grep _R1.fq | sed 's/_R1.fq//' | \
parallel -j10 rsem-calculate-expression --paired-end --bowtie2 --bowtie2-sensitivity-level sensitive -p 10 \
/projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/{}_R1.fq \
/projects/glomicave/metatranscriptomics/01_mRNA_processed/processing_ribodetector/mRNA/{}_R2.fq \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/reference/combined.ORFs \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_old/results/{}