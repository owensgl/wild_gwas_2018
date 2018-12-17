library(bigsnpr)
library(tidyverse)

tmpfile <- tempfile()
snp_readBed("/media/owens/Copper/wild_gwas_2018/annuus/Annuus.tranche90.snp.gwas.90.bi.beagle.bed", backingfile = tmpfile)
obj.bigSNP <- snp_attach(paste0(tmpfile, ".rds"))
obj.bigSNP$map$chromosome <- as.numeric(gsub("HanXRQChr", "", obj.bigSNP$map$chromosome))

gwas_2PC <- read_csv("/media/owens/Copper/wild_gwas_2018/gwas/00ea2e84-5252-406e-8947-5c19168f7446.csv",comment="#") %>%
  rename(chr= `>CHR`,pos=Positions) %>% 
  select(chr,pos,TestStatistic,Beta0,Beta1, MAF) %>%
  mutate(chr= gsub("HanXRQChr","",chr),location=paste(chr,pos,sep="."), chr= as.numeric(chr))

gwas_20PC <- read_csv("/media/owens/Copper/wild_gwas_2018/gwas/bbb9c329-3dae-46d1-ad6f-d66a7645146c.csv",comment="#") %>%
  rename(chr= `>CHR`,pos=Positions) %>% 
  select(chr,pos,TestStatistic,Beta0,Beta1, MAF) %>%
  mutate(chr= gsub("HanXRQChr","",chr),location=paste(chr,pos,sep="."), chr= as.numeric(chr))

#Comparing GWAS


inner_join(gwas_2PC %>% select(chr,pos,Beta1) %>% rename(Beta_2PC=Beta1),
           gwas_20PC %>% select(chr,pos,Beta1) %>% rename(Beta_20PC=Beta1)) %>%
  mutate(beta_dif = Beta_2PC - Beta_20PC) %>%
  ggplot(.,aes(x=beta_dif)) + geom_histogram() 


snp.maf <- snp_MAF(obj.bigSNP$genotypes)

chr <- obj.bigSNP$map$chromosome
pos <- obj.bigSNP$map$physical.pos

hist(colSums(obj.bigSNP$genotypes[,1:length(chr)]) / 1228)
snp.maf[1]
original_data <- tibble(chr = chr, pos = pos, snp.maf = snp.maf)

full_join(original_data, gwas_2PC) -> original_data

original_data %>%
  filter(is.na(TestStatistic)) %>%
  filter(snp.maf >= 0.05) %>% View()
  ggplot(aes(snp.maf)) + geom_histogram()

subset_data <- subset(obj.bigSNP, ind.col = !is.na(original_data$TestStatistic), backingfile = subset_file)
obj.smallSNP <- snp_attach(subset_data)
G <- obj.smallSNP$genotypes

gwas %>% filter(pos == 167559702)

ind.keep2 <- snp_clumping(G, infos.chr = obj.smallSNP$map$chromosome, thr.r2 = 0.1, S=gwas_2PC$TestStatistic,
                          infos.pos = obj.bigSNP$map$physical.pos,
                          ncores=10)

sum(pull(gwas_2PC[ind.keep2,5]) * round(jitter(obj.smallSNP$genotypes[5,ind.keep2]/2,amount=0.01)))
sum(abs(obj.smallSNP$genotypes[40,ind.keep2] - 2) * pull(gwas_2PC[ind.keep2,5]))

abs(obj.smallSNP$genotypes[2,ind.keep2] - 2)


test <- snp_attachExtdata()
G <- big_copy(test$genotypes, ind.col = 1:1000)
CHR <- test$map$chromosome[1:1000]
POS <- test$map$physical.pos[1:1000]
y01 <- test$fam$affection - 1

# PCA -> covariables
obj.svd <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
#> Phase of clumping at r2 > 0.2.. keep 952 SNPs.
#> 
#> Iteration 1:
#> Computing SVD..
#> 
#> Converged!

# train and test set
ind.train <- sort(sample(nrow(G), 400))
ind.test <- setdiff(rows_along(G), ind.train) # 117

# GWAS
gwas.train <- big_univLogReg(G, y01.train = y01[ind.train],
                             ind.train = ind.train,
                             covar.train = obj.svd$u[ind.train, ])
# clumping
ind.keep <- snp_clumping(G, infos.chr = CHR,
                         ind.row = ind.train,
                         S = abs(gwas.train$score))
# -log10(p-values) and thresolding
summary(lpS.keep <- -predict(gwas.train)[ind.keep])
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.000565 0.132700 0.322159 0.482021 0.668394 7.511245 
thrs <- seq(0, 4, by = 0.5)
nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))

# PRS
prs <- snp_PRS(G, betas.keep = gwas.train$estim[ind.keep],
               ind.test = ind.test,
               ind.keep = ind.keep,
               lpS.keep = lpS.keep,
               thr.list = thrs)

# AUC as a function of the number of predictors
aucs <- apply(prs, 2, AUC, target = y01[ind.test])
library(ggplot2)
bigstatsr:::MY_THEME(qplot(nb.pred, aucs)) +
  geom_line() +
  scale_x_log10(breaks = nb.pred) +
  labs(x = "Number of predictors", y = "AUC")


cmsa.logit <- big_spLogReg(X = G, y01.train = y01[ind.train], 
                           ind.train = ind.train, 
                           covar.train = obj.svd$u[ind.train, ],
                           alphas = c(1, 0.5, 0.05, 0.001),
                           ncores = 2)

preds <- predict(cmsa.logit, X = G, ind.row = ind.test, 
                 covar.row = obj.svd$u[ind.test, ])



