

#Write while loop that goes through every file in a directory and creates a symlink of every file with a number between 151 and 189 (and .bam at the end) in the name to a new directory
# Path: symlink_between_numbers.sh
source="/projects/glomicave/metatranscriptomics/03_coverm/bam_filtered/"
destination="/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/bam_files/"

counter=151

while [ $counter -le 189 ]; do
    ln -s $source*$counter*.bam $destination
    ((counter++))
done