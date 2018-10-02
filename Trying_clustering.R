library(kmodR)
library(mvoutlier)
library(cclust)
win.regions <- readRDS(paste(data_name,".",window_size,".windows.rds",sep=""))

win.regions %>% select(mds01:mds20)

tmp <- win.regions
k5 <-kmeans(win.regions %>% select(mds01:mds20),5 )
k1 <-kmod(win.regions %>% select(mds01:mds20),1, l=200)

aq <- aq.plot(win.regions %>% select(mds01:mds20),alpha=0.00001)

tmp$cluster <- aq$outliers


tmp %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-2)]) %>% 
  ggplot(aes(x=mid,y=value,color=cluster)) + geom_point(alpha=0.5) +
  facet_grid(mds~chrom,scales="free") + theme_bw() +
  scale_color_brewer(palette = "Set1")


for (i in 2:20){
  cc <- cclust(as.matrix(tmp %>% filter(cluster == "TRUE") %>% select(mds01:mds20)),i,method="neuralgas")
  print(i)
  print(
  clustIndex(cc,as.matrix(tmp %>% filter(cluster == "TRUE") %>% select(mds01:mds20)), index="all")[1]
  )
}
cc <- cclust(as.matrix(tmp %>% filter(cluster == "TRUE") %>% select(mds01:mds20)),8,method="neuralgas")


pdf("example_cluster_arg.pdf",height=16,width=50)
full_join(tmp, tmp %>% filter(cluster == "TRUE") %>% cbind(cc$cluster) %>%
            rename(kmeans_cluster = 'cc$cluster')) %>%
  gather(., mds, value, colnames(win.regions)[5:(ncol(win.regions)-1)]) %>%
  ggplot(aes(x=mid,y=value,color=as.factor(kmeans_cluster))) + geom_point(alpha=1) +
  facet_grid(mds~chrom,scales="free") + theme_bw() 
dev.off() 

full_join(tmp, tmp %>% filter(cluster == "TRUE") %>% cbind(cc$cluster) %>%
            rename(kmeans_cluster = 'cc$cluster')) %>%
  filter(cluster == "TRUE") %>% 
  ggplot(.,aes(chrom)) + geom_histogram(aes(fill=chrom),stat="count") +
  facet_wrap(~kmeans_cluster,scales="free_y")


data(X)
data(Y)
data(dat)
res <- locoutNeighbor(dat,X,Y,variant="knn",usemax=1,chisqqu=0.975,indices=c(1,11,24,36),
                      propneighb=0.1,npoints=100)


library(flexmix)
data("Nclus")

mymclust <- function ( formula = .~. , diagonal = TRUE )
{
  retval <- new ("FLXMC" , weighted = TRUE ,
                 formula = formula , dist = "mvnorm",
                 name = " my model - based clustering ")
  retval@defineComponent <- function ( para ) {
    logLik <- function (x , y ) {
      mvtnorm::dmvnorm (y , mean = para$center ,
                          sigma = para$cov , log = TRUE )
    }
    predict <- function ( x) {
      matrix ( para$center , nrow = nrow(x ),
               ncol = length ( para$center ), byrow = TRUE )
    }
    new ("FLXcomponent", parameters =
           list ( center = para$center , cov = para$cov ),
         df = para$df , logLik = logLik , predict = predict )
  }
  retval@fit <- function (x , y , w , ...) {
    para <- cov.wt (y , wt = w )[ c ("center", " cov ")]
    df <- (3 * ncol (y) + ncol (y )^2)/2
    if ( diagonal ) {
      para$cov <- diag ( diag ( para$cov ))
      df <- 2 * ncol ( y)
    }
    retval@defineComponent (c ( para , df = df ))
  }
  retval
}

ggplot(data=as.tibble(Nclus)) + geom_point(aes(x=V1,y=V2))
m1 <- flexmix(Nclus ~ 1, k = 4, model = FLXmclust)


library(ggrepel)

windows <- win.regions %>%
  mutate_(the_mds = "mds10-pos" ) %>% 
  filter(the_mds < high_cutoff & chrom == cluster_chr) %>% pull(n)

dist(win.fn.snp(windows))

out <- cov_pca(win.fn.snp(windows),k=2)
matrix.out <- t(matrix(out[4:length(out)],ncol=nrow(samples),byrow=T))
out <- matrix(out[4:length(out)],ncol=nrow(samples),byrow=T) %>% as.tibble() 
colnames(out) <- pull(samples)

out <- as_tibble(cbind(nms = names(out), t(out))) %>% 
  rename(name=nms,PC1=V1,PC2=V2) %>% 
  mutate(PC1 = as.double(PC1), PC2 = as.double(PC2))

out %>%
  mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)) %>%
  ggplot(.,aes(x=PC1,y=PC2,label=n)) + geom_point() + theme_bw() +
  geom_text_repel()

tmp <- hclust(dist(out[2:3]), method="complete")

out %>% pull(PC1) %>%  as.numeric(.) -> tmp_PC1
centers_matrix <- matrix(data=c(0,0.1,0.15,-0.1,0.3,-0.2),nrow=3)
kmeans_cluster <-kmeans(tmp_PC1, 3, centers=))
kmeans_cluster <-kmeans(matrix.out, 3,centers=centers_matrix)


try <- flexmix(PC1 ~ PC2, data = out %>%
                 mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)), k = 3)

var_value <- c(FALSE,F)
flexmix(.~., k = 3, model = FLXMRlmm(PC1 ~ PC2, random = ~ 1,varFix=c(Random = T, Residual = T)),
        data = out %>%
          mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)), 
        cluster=kmeans_cluster$cluster)
out$cluster <- try@cluster
out$cluster <-kmeans_cluster$cluster
out %>%
  mutate(PC1 = as.double(PC1), PC2 = as.double(PC2)) %>%
  ggplot(.,aes(x=PC1, y=PC2,color=as.factor(cluster))) + geom_point()

Nclus

example <- as.tibble(Nclus)
example$cluster <- ex1@cluster

example %>%
  ggplot(.,aes(x=V1, y=V2,color=as.factor(cluster))) + geom_point()


id <- rep(1:50, each = 10)
x <- rep(1:10, 50)
sample <- data.frame(y = rep(rnorm(unique(id)/2, 0, c(5, 2)), each = 10) +
                       rnorm(length(id), rep(c(3, 8), each = 10)) +
                       rep(c(0, 3), each = 10) * x,
                     x = x,
                     id = factor(id))
fitted <- flexmix(.~., k = 2, model = FLXMRlmm(y ~ x, random = ~ 1),
                  data = sample, control = list(tolerance = 10^-3),
                  cluster = rep(rep(1:2, each = 10), 25))
1.58
rotation.ss <- tibble(rotation = numeric(),betweenss = numeric())
for (i in seq(0.02, 3.14, 0.02)){
  rotated.matrix <- Rotation(matrix.out,i)
  
  rotated.kmeans <- kmeans(rotated.matrix[,1], 3, centers=c(min(rotated.matrix[,1]),(min(rotated.matrix[,1])+max(rotated.matrix[,1]))/2,max(rotated.matrix[,1])))
  rotated.tibble <- tibble(PC1 = as.numeric(rotated.matrix[,1]),
                            PC2 = as.numeric(rotated.matrix[,2]),
                            cluster = rotated.kmeans$cluster)

 # ggplot(rotated.tibble,aes(x=PC1,y=PC2,color=as.factor(cluster))) + 
  #  geom_point() + scale_color_brewer(palette = "Set1") + ggtitle(paste("BetweenSS =",round(rotated.kmeans$betweenss,3)))
  tmp.tibble <- tibble(rotation = as.numeric(i),betweenss = as.numeric(rotated.kmeans$betweenss))
  rotation.ss <- rbind(rotation.ss, tmp.tibble)
}
plot(rotation.ss)
