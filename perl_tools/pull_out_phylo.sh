i=$1
module load bcftools
module load tabix
chr="XXXchr"
mds="XXXmds"
git_species="XXXgit_species"
all_species="XXXall_species"
short_species="XXXshort_species"
samplefile="$chr.$mds.phylo.$git_species.txt"
full_samplefile="$chr.$mds.phylo.samplefile.txt"
directory="/scratch/gowens/wild_gwas"
outgroups="extra_samples_outgroups.txt"
genotyped_samplelist="$directory/$all_species/Texanus.tranche90.snp.gwas.90.bi.samplelist.txt"
if [ ! -f $samplefile ]
then
	perl /scratch/gowens/wild_gwas/wild_gwas_2018/perl_tools/pull_out_inv_homo_for_phylo.pl /scratch/gowens/wild_gwas/wild_gwas_2018/MDS_outliers/Ha412HO/$git_species/Ha412HO_inv.jan09.pcasites.$chr.$mds.genotypes.txt $genotyped_samplelist $short_species > $samplefile
fi
sleep 5s #Used so that the samplefile is present
cat $samplefile ../$outgroups > $full_samplefile.$i
bcftools view  -O v -a -S $full_samplefile.$i -c 1 ../$i.tranche90.snp.vcf.gz  > $i.tranche90.snp.$chr.$mds.all.vcf
for sample in `cat ../$outgroups`
do
	echo $sample > $chr.$mds.$sample.$i.samplelist.txt
	cat $samplefile >> $chr.$mds.$sample.$i.samplelist.txt
	bcftools view -S $chr.$mds.$sample.$i.samplelist.txt -c 1 -a -O v $i.tranche90.snp.$chr.$mds.all.vcf | bcftools view -c 1 -a -i 'INFO/AN>5' -O z > $i.tranche90.snp.$chr.$mds.$sample.vcf.gz

	tabix -p vcf $i.tranche90.snp.$chr.$mds.$sample.vcf.gz
	rm $chr.$mds.$sample.$i.samplelist.txt
done
rm $full_samplefile.$i
