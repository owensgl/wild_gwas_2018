library(inveRsion)
genofile <- "/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.env.90.bi.HanXRQChr05.txt"
gDat <-setUpGenoDatFile(file=genofile,sortMinor=TRUE,saveRes=FALSE)
haplo <- codeHaplo(gDat, 100, saveRes = F, file = NULL, intSNP=FALSE,phasing="inversion&BP")
scanRes<-scanInv(haplo,window=0.5,saveRes=TRUE,geno=T,)

