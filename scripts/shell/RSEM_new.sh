#Preparing the reference for RSEM - this is with Ziville and Caitlin's extra rRNA filtering
rsem-prepare-reference --bowtie2 /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/combined.ffn /user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/combined.ORFs

#Running rsem-calculate-expression with 10 threads parallel (5 threads)
ls /projects/glomicave/metatranscriptomics/03_coverm/reads_filtered/ | grep _R1.fastq | sed 's/_R1.fastq//' | \
parallel -j10 rsem-calculate-expression --paired-end --bowtie2 --bowtie2-sensitivity-level sensitive -p 5 \
/projects/glomicave/metatranscriptomics/03_coverm/reads_filtered/{}_R1.fastq \
/projects/glomicave/metatranscriptomics/03_coverm/reads_filtered/{}_R2.fastq \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/reference/combined.ORFs \
/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/{}