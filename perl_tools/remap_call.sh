module load bcftools
module load tabix
module load bwa
module load samtools
module load bedtools
chr="XXXchr"
mds="XXXmds"
capital_species="XXXcapital_species"
git_species="XXXgit_species"
all_species="XXXall_species"
samplefile="$chr.$mds.sampleinfo.txt"
prefix="tranche90.snp.$chr.$mds"
pop=$1

#Merge all chromosomes including 00c
if [ ! -f $capital_species.$prefix.$pop.tmp.vcf.gz ]
then
	bcftools concat --threads 10 -O z -a Han*.$prefix.$pop.vcf.gz > $capital_species.$prefix.$pop.tmp.vcf.gz
fi

#Remap to Ha412HO
if [ ! -f $capital_species.$prefix.$pop.remappedHa412.vcf.gz ]
then
	perl /scratch/gowens/wild_gwas/wild_gwas_2018/perl_tools/xrqpos2ha412pos_bwa.pl $capital_species.$prefix.$pop.tmp.vcf.gz $capital_species.$prefix.$pop
	tabix -p vcf $capital_species.$prefix.$pop.remappedHa412.vcf.gz
fi

#Remove extra contigs

bcftools view --threads 10 -R /scratch/gowens/wild_gwas/wild_gwas_2018/Ha412HO.chrlist.bed -O v $capital_species.$prefix.$pop.remappedHa412.vcf.gz  | perl /scratch/gowens/pop_gen/vcf2ibd_dist.pl > /scratch/gowens/wild_gwas/wild_gwas_2018/Inv_ibd/Ha412HO/$git_species/$capital_species.$prefix.$pop.remappedtarget.ibd.txt



rm $capital_species.$prefix.$pop.Ha412conversionstats.txt
rm $capital_species.$prefix.$pop.unsorted.vcf
rm $capital_species.$prefix.$pop.sam
rm $capital_species.$prefix.$pop.fast*
rm $capital_species.$prefix.$pop.bed
rm $capital_species.$prefix.$pop.tmp.vcf.gz
rm Han*.$prefix.$pop.vcf.gz*


