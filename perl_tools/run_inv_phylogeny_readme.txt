This is the outline of what run_inv_phylogeny.sh does and how to run it. 



Usage: bash run_inv_phylogeny.sh #
-# is the line number of the inversion document (e.g. Ha412HO_inv.v3.txt), which is where it gets all the information about species and chromosome. Line # 1 is the header and won't work.

Basic outline:
-The script pulls out random samples that are with 0.15 of purity based on pcasites (e.g. 0-0.15, 0.85-1) and also called homozygous. It selects 5 samples per inv
ersion genotype. It also selects 2 individuals of all other species, and two outgroup samples.
-It finds XRQ genes in the inversion by finding HA412 genes in the inversion and then selecting their 1to1 orthologs
-For each gene, it produces a gvcf (with information of each site). That gvcf is converted to fasta while converting sites <2 reads and indels to "N". The fastas of each sample are combined.
-If requested, it tests possible inversion gene trees and outputs the support for each hypothesis
-At the end, all fastas are concatenated and a partition file is produced saying where genes start and end in the larger file.
-If requested, it creates a maximum likelihood tree of the full multigene fasta


Notes:
-It's designed for Cedar and uses the gatk4 module.
-The non-inversion samples are predetermined and in the resource directory (e.g. annuus.physamples.txt). You can change them if you want.
-The "make_trees" parameter will decide if it makes a tree of the concatenated fastas at the end. Set to "TRUE" if you want it to make a tree
-The "test_trees" parameter will decide if it tests inversion trees for each gene. 

