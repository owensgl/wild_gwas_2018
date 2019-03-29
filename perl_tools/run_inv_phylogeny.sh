#Takes an inversion. Pulls out inversion regions in Ha412, gets the genes in that region, then pulls out orthologous genes in XRQ. Then we use GATK to pull out fasta sequences for a representattive set and calculate a phylogeny.
#Usage example: bash run_inv_phylogeny.sh 2

inv_n=$1 #The line number from the inversion sheet, starts with 2


#Parameters to update 
working_directory="/scratch/gowens/wild_gwas/temp" #Where all the output files are going to be produced
wild_gwas_git="/scratch/gowens/wild_gwas/wild_gwas_2018/"
inversion_version="v3" #Examples mar27, v3.
iqtree="/home/gowens/bin/iqtree-1.6.10-Linux/bin/iqtree" #IQ tree executable
outgroups="GRO_2043,DIV_1956"
xrq_ref="/scratch/gowens/ref/HanXRQr1.0-20151230.fa"
n_cores=30 #Number of cores for IQtree
make_trees="TRUE" #Set to TRUE if you want to make a phylogeny in IQ tree at the end
test_trees="TRUE" #Set to TRUE if you want to test inversion phylogenies at each gene.

cd $working_directory
inv_page="$wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.txt"
species=$(cat $inv_page | sed -n ${inv_n}p | cut -f 1)
chr_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 4)
chr_padded=$(printf %02d $chr_n)
chr="Ha412HOChr$chr_padded"
direction=$(cat $inv_page | sed -n ${inv_n}p | cut -f 6)
mds_n=$(cat $inv_page | sed -n ${inv_n}p | cut -f 5)
mds=${direction}$mds_n
Ha412_genes="$wild_gwas_git/resources/Han412-HO_gene_filt.gff3.gz"

#Samples from other species are pre selected and stored in these files
other_samples="$wild_gwas_git/resources/$species.physamples.txt"

#Make directories
mkdir -p $species
cd $species
mkdir -p "$chr.${mds}"
cd "$chr.${mds}"

species_abbreviation=""
#Make species abbreviations
if [ $species == "annuus" ]
then
	species_abbreviation="Ann"
elif [ $species == "argophyllus" ]
then
	species_abbreviation="Arg"
elif [ $species == "petpet" ]
then
	species_abbreviation="Pet"
elif [ $species == "petfal" ]
then
	species_abbreviation="Pet"
fi
if [ ! -f group_2.txt ]
then
#Get list of samples appropriate for each inversion test
cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 2) && ($4 > 0.9)' | shuf | head -n 5 | cut -f 1 > group_2.txt
cat $wild_gwas_git/MDS_outliers/Ha412HO/$species/Ha412HO_inv.$inversion_version.pcasites.$chr.$mds.genotypes.txt | awk '($2 != "NA") && ($6 == 0) && ($4 < 0.1)' | shuf | head -n 5 | cut -f 1 > group_0.txt
cat group_0.txt group_2.txt $other_samples > $chr.$mds.samplelist.txt

#Get HA412 genes in inversion region
cat $wild_gwas_git/MDS_outliers/Ha412HO/all/Ha412HO_inv.$inversion_version.inversions.regions.v1.txt | grep $species | grep $chr | grep $mds | cut -f 2-4 > inversion.regions.txt
n_regions=$(wc -l inversion.regions.txt | cut -f 1 -d " ")
rm inversion.ha412.genes.txt
for n in `seq $n_regions`
do
	region_chr=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 3)
	region_start=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 1)
	region_end=$(cat inversion.regions.txt | sed -n ${n}p | cut -f 2)
	zcat $Ha412_genes | perl $wild_gwas_git/perl_tools/pull_out_genes.pl $region_chr $region_start $region_end >> inversion.ha412.genes.txt
done	
#Get orthologous XRQ genes
cat $wild_gwas_git/resources/Ha412HO_HanXRQv1.1to1ortho.txt | grep -f inversion.ha412.genes.txt | cut -f 2 | tr -d ' '  > inversion.xrq.genes.txt
cat $wild_gwas_git/resources/HanXRQr1.0-20151230-EGN-r1.1.genelocations.txt | grep -f inversion.xrq.genes.txt  > inversion.xrq.genes.locations.txt

fi
if [ $test_trees == "TRUE" ]
then
#Make tested trees
	perl $wild_gwas_git/perl_tools/create_test_phylogenies.pl $species_abbreviation group_0.txt group_2.txt $other_samples $wild_gwas_git/resources/sample_info_file_all_samples_2018_12_07.tsv > $species.testtrees.treel
fi
#Make all the fasta files
bash $wild_gwas_git/perl_tools/make_all_fasta.sh $chr.$mds.samplelist.txt $species.testtrees.treel $wild_gwas_git/perl_tools $iqtree $outgroups $test_trees $xrq_ref $wild_gwas_git/resources/sample_info_file_all_samples_2018_12_07.tsv
ls | grep merged | grep fasta | perl $wild_gwas_git/perl_tools/concatenate_fasta.pl > $chr.$mds.fasta
if [ $make_trees == "TRUE" ]
then
	$iqtree -s $chr.$mds.fasta -bb 1000 -q partions.txt -o $outgroups -nt $n_cores -redo -pre $chr.$mds -alrt 1000
fi
