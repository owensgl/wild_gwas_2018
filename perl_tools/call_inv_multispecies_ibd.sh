#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=2:00:00
#SBATCH --job-name=Inv_ibd
#SBATCH --array=3-11%50
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --exclusive=user
#SBATCH --output=/home/gowens/bin/logs/%A.%a.invibd_log

n="$SLURM_ARRAY_TASK_ID"
cd /scratch/gowens/wild_gwas/wild_gwas_2018/perl_tools/

git_species=$(cat /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt | sed '1d' | sed -n ${n}p | cut -f 1)
capital_species=$(cat /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt | sed '1d' | sed -n ${n}p | cut -f 2)
chr_raw=$(cat /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt | sed '1d' | sed -n ${n}p | cut -f 4)
chr_n=$(printf "%02d" $chr_raw)
chr="Ha412HOChr$chr_n"
mds_n=$(cat /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt | sed '1d' | sed -n ${n}p | cut -f 5)
mds_dir=$(cat /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.jan09.txt | sed '1d' |  sed -n ${n}p | cut -f 6)
mds="$mds_dir$mds_n"


all_species="texanus" #Only using texanus snp set for now.
species_abbreviation=""
if [ $git_species == "annuus" ]
then
#	all_species="all_annuus"
	species_abbreviation="Ann"
elif [ $git_species == "anomalus" ]
then
#	all_species="all_ano"
	species_abbreviation="Ano"
elif [ $git_species == "argophyllus" ]
then
#	all_species="all_arg"
	species_abbreviation="Arg"
elif [ $git_species == "petpet" ]
then
#	all_species="all_petiolaris"
	species_abbreviation="PetPet"
elif [ $git_species == "petfal" ]
then
#	all_species="all_petiolaris"
	species_abbreviation="PetFal"
elif [ $git_species == "petiolaris" ]
then
#	all_species="all_petiolaris"
	species_abbreviation="Pet"
elif [ $git_species == "niveus" ]
then
#	all_species="texanus"
	species_abbreviation="PetCan"
fi


#Shove variables into bash scripts so that we can call them using gnu parallel
directory="/scratch/gowens/wild_gwas/"
new_dir=$directory/${all_species}/$species_abbreviation.${chr}.${mds}.phylo
mkdir -p $new_dir
cat pull_out_phylo.sh | sed s/XXXchr/$chr/g | sed s/XXXmds/$mds/g | sed s/XXXgit_species/$git_species/g | sed s/XXXall_species/$all_species/g | sed s/XXXshort_species/$species_abbreviation/g > $new_dir/pull_samples.$git_species.$chr.$mds.sh

cat remap_phylo.sh | sed s/XXXchr/$chr/g | sed s/XXXmds/$mds/g | sed s/XXXcapital_species/$capital_species/g | sed s/XXXall_species/$all_species/g | sed s/XXXgit_species/$git_species/g > $new_dir/remap_call.$git_species.$chr.$mds.sh

#Run modified scripts in parallel
cd $new_dir
cat /scratch/gowens/wild_gwas/wild_gwas_2018/chrlist.txt | parallel -j 18 bash pull_samples.$git_species.$chr.$mds.sh
cat "../extra_samples_outgroups.txt" | cut -f 2 | sort | uniq | parallel -j 3 bash remap_call.$git_species.$chr.$mds.sh