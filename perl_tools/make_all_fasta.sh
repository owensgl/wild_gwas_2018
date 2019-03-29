inversion_genes="inversion.xrq.genes.locations.txt"
inversion_samplelist=$1
test_trees=$2
script_storage=$3 
iq_tree=$4
outgroups=$5
test_trees=$6
ref=$7
sampleinfo=$8
make_fasta_cores=8
iqtree_cores=30

module load gatk
n_genes=$(wc -l $inversion_genes | cut -f 1 -d " ")


for i in `seq $n_genes`
do
  if [ ! -s "merged.$i.fasta" ]
  then
    cat $inversion_genes | sed -n ${i}p | cut -f 2- > gene.bed
    cat $inversion_samplelist | parallel -j $make_fasta_cores bash $script_storage/make_fasta.sh $sampleinfo $ref $script_storage
    for f in *.fa; do (cat "${f}"; echo); done > merged.$i.fasta
    rm *.g.vcf*
    rm *.fa
    if [ $test_trees == "TRUE"]
    then
      $iq_tree -s merged.$i.fasta -nt $iqtree_cores -z $test_trees -n 0 -zb 10000 -au -o $outgroups -redo
      cat merged.$i.fasta.iqtree | perl $script_storage/parse_iqtree_test.pl > merged.$i.treetests.txt
      rm merged.$i.fasta.*
    fi
  fi
done

