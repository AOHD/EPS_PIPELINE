#!/usr/bin/env bash
#
module load GTDBTk/1.5.0-foss-2020b-Python-3.8.6
gtdbtk de_novo_wf --bacteria --genome_dir /user_data/ahd/EPS_PIPELINE/data/raw/bins --outgroup_taxon p__Patescibacteria --out_dir /user_data/ahd/EPS_PIPELINE/tree_results/gtdb_16112022 -x fa --cpus 8
module purge