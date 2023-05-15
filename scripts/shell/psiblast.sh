#!/usr/bin/env bash
data="/user_data/ahd/EPS_PIPELINE/databases/mydb"
query="/user_data/ahd/EPS_PIPELINE/data/raw/queries_psiblast"
WD="/user_data/ahd/EPS_PIPELINE"
cd $WD
mkdir $WD/"psiblast_results"
mkdir $WD/"psiblast_results/10_Nov"
##Query path definitions, non-MSA queries
operon_fasta=(
#acetan.faa
#alginate.faa
#amylovoran.faa
#cellulose.faa
#cellulose1.faa
#cellulose2.faa
#ColA.faa
#curdlan.faa
#diutan.faa
#gellan.faa
#gellan1.faa
#gellan2.faa
#HA_Pasteurella.faa
#HA_streptococcus.faa
#NulO_merged.faa
#pel_merged.faa
#pnag_eps.faa
#pnag_ica.faa
#pnag_pga.faa
#psl.faa
#rhizobium_eps.faa
#s88.faa
#salecan.faa
#stewartan.faa
#succinoglycan.faa
#vps.faa
#xanthan.faa
#burkholderia_eps.faa
#levan.faa
#synechan.faa
#methanolan.faa
#galactoglucan.faa
B_fragilis_PS_A.faa
B_fragilis_PS_B.faa
B_pseudomallei_EPS.faa
cepacian.faa
E_faecalis_PS.faa
emulsan.faa
EPS273.faa
GG.faa
glucorhamnan.faa
L_johnsonii_ATCC_33200_EPS_A.faa
L_johnsonii_ATCC_11506_EPS_B.faa
L_johnsonii_ATCC_2767_EPS_C.faa
L_lactis_EPS.faa
L_plantarum_HePS.faa
phosphonoglycan.faa
)
## PSI-BLAST of all operons in the GFF database from HQ-MAGs
module load BLAST+/2.12.0-gompi-2020b
for operon in ${operon_fasta[@]}; do
psiblast -query $query/$operon -db $data -out $WD/psiblast_results/$date/$operon -evalue 0.0001 -qcov_hsp_perc 20 -max_hsps 10 -max_target_seqs 100000 -outfmt 6 -num_iterations 20 -comp_based_stats 1 -num_threads 40
echo $operon BLASTed
done
module purge