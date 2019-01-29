i=$1
module load bcftools
module load tabix
chr="XXXchr"
mds="XXXmds"
git_species="XXXgit_species"
all_species="XXXall_species"
samplefile="$chr.$mds.sampleinfo.txt"
directory="/scratch/gowens/wild_gwas/"
cut -f 1 $samplefile > $chr.$mds.samples.txt
bcftools view  -O v -a -S $chr.$mds.samples.txt -c 1 ../$i.tranche90.snp.vcf.gz > $i.tranche90.snp.$chr.$mds.all.vcf
cut -f 2 $samplefile | sort | uniq > $chr.$mds.pops.txt

for pop in `cat $chr.$mds.pops.txt`
do
	grep $pop $samplefile | cut -f 1 > $chr.$mds.$pop.samples.txt

	cat $i.tranche90.snp.$chr.$mds.all.vcf | bcftools view -S $chr.$mds.$pop.samples.txt -c 1 -a -i 'INFO/AN>7' -O z > $i.tranche90.snp.$chr.$mds.$pop.vcf.gz
	tabix -p vcf $i.tranche90.snp.$chr.$mds.$pop.vcf.gz
done
rm $i.tranche90.snp.$chr.$mds.all.vcf
