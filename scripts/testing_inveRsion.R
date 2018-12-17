library(inveRsion)
genofile <- "/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.beagle.remappedHa412.Ha412TMPChr13.txt"
gDat <-setUpGenoDatFile(file=genofile,sortMinor=TRUE,saveRes=FALSE)
haplo <- codeHaplo(gDat, 5, saveRes = TRUE, file = NULL, intSNP=FALSE,phasing="inversion&BP")
scanRes<-scanInv(haplo,window=0.5,saveRes=TRUE,geno=T)

load("scripts/scanRes.RData")


plot(scanRes,which="b",thBic=-Inf)
