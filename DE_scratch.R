#' Differential Expression from scratch with Single Cell RNA Sequencing data
#' Notes from Kamil (4/12/19)
library(ggplot)
library(ggrepel)

# Load UMI count data into dataframe (samples, groups, counts)
AllCounts <- readRDS("compbio/fatfly_zumi_test_run.dgecounts.rds")
umi <- as.matrix(AllCounts$umicount$exon$all) 
samples <- c("sta-gfp-2R", "fed-neg-2L", "sta-neg-2L", "fed-gfp-1R",
             "fed-neg-1R", "fed-gfp-1L", "sta-gfp-1L", "sta-gfp-2L", 
             "sta-neg-1L", "fed-gfp-2R", "fed-neg-2R", "fed-gfp-2L", 
             "fed-neg-1L", "sta-neg-1R", "sta-neg-2R", "sta-gfp-1R")
colnames(umi) <- samples
groups <- stringr::str_split_fixed(colnames(log2cpm), '-', 3)[,1:2]
counts <- colSums(umi)
annot <- data.frame(samples, groups, counts)

# log transform cpm for every gene (log base 2) ~ tends to fit values toward a normal distribution
hist(umi, main="Raw UMI counts by gene", xlab="counts") # raw
hist((umi+1)/1e6, main="CPM by gene") # cpm 
log2cpm <- t(apply(umi, 1, function(x) log2(((x+1)/sum(x)) * 1e6))) #log2(cpm)
hist(log2cpm, main = "Log2(CPM) by gene") 

# principal components analysis of ALL samples
pca <- prcomp(t(log2cpm))
pca <- summary(pca)
pca_df <- cbind(annot, pca$x)

## scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca_names <- factor(colnames(pca$x), levels = colnames(pca$x))
df <- data.frame(pca.var.per, pca_names)
ggplot(data=df, aes(x=pca_names, y=pca.var.per)) + geom_bar(stat="identity")

## plot PC2 vs. PC1 (color by sta/fed)
ggplot(data=pca_df, aes(x=PC1, y=PC2, label=samples, color=X1)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1 vs. PC2 (fed vs. starved)")

## plot PC2 vs. PC1 (color by gfp/neg)
ggplot(data=pca_df, aes(x=PC1, y=PC2, label=samples, color=X2)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1 vs. PC2 (gfp+ vs. gfp-)")

# model correlation of each meta variable with PC1/PC2
lmPC1 <- lm(PC1 ~ X2, pca_df)
lmPC2 <- lm(PC2 ~ X1+X2+counts, pca_df)
summary(lmPC1)
summary(lmPC2)

## identify genes with highest loading values for PC1 (gfp+ vs gfp-)
diffExp <- sort(pca$rotation[,1])
diffExp[(15245-50):15245] 

## look for genes with a given name
#names(diffExp)[grepl("NPF", names(diffExp))]

## plot gene expression (cpm/log2cpm) across samples
gene.df <- data.frame(samples=colnames(log2cpm),
                      log2cpm=log2cpm["NPFR",])
ggplot(data=gene.df, aes(x=samples, y=log2cpm)) + geom_bar(stat="identity")

reversed <- t(apply(log2cpm, 1, function(x) 2^x))
gene.df <- data.frame(samples = colnames(reversed),
                      cpm=reversed["NPFR",])
ggplot(data=gene.df, aes(x=samples, y=cpm)) + geom_bar(stat="identity")


# perform PCA (GFP ONLY)
samples.gfp <- grepl("gfp",colnames(log2cpm))
log2cpm.gfp <- log2cpm[,samples.gfp]
annot.gfp <- data.frame(samples = annot$samples[samples.gfp],
                        groups = annot$X1[samples.gfp],
                        count = annot$counts[samples.gfp])
pca.gfp <- prcomp(t(log2cpm.gfp))
pca.gfp <- summary(pca.gfp)
pca.gfp_df <- cbind(annot.gfp, pca.gfp$x)

## scree plot (GFP ONLY)
pca.gfp.var <- pca.gfp$sdev^2
pca.gfp.var.per <- round(pca.gfp.var/sum(pca.gfp.var)*100, 1)
pca_names.gfp <- factor(colnames(pca.gfp$x), levels = colnames(pca.gfp$x))
df.gfp <- data.frame(pca.gfp.var.per, pca_names.gfp)
ggplot(data=df.gfp, aes(x=pca_names.gfp, y=pca.gfp.var.per)) + geom_bar(stat="identity")

## plot PC2 vs. PC2, GFP ONLY (fed vs. starved)
ggplot(data=pca.gfp_df, aes(x=PC1, y=PC2, label=samples, color=groups)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.gfp.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.gfp.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("PC1 vs. PC2, GFP ONLY (fed vs. starved)")

## model correlation of various features with PC1/PC2
lmGFP.PC1 <- lm(PC1 ~ count, pca.gfp_df)
lmGFP.PC2 <- lm(PC2 ~ groups, pca.gfp_df)
summary(lmGFP.PC1)
summary(lmGFP.PC2)

## identify genes with highest loading values for PC2 (fed vs. starved)
diffExp <- sort(pca.gfp$rotation[,2])
diffExp[(15245-100):15245] # up-regulated in fed
diffExp[1:40] # up-regulated in starved

## plot log2(cpm) gene expression across samples
grouped.log2cpm.gfp <- log2cpm.gfp[,c("fed-gfp-1R", "fed-gfp-1L", "fed-gfp-2R", "fed-gfp-2L",
                                      "sta-gfp-1R","sta-gfp-1L", "sta-gfp-2R", "sta-gfp-2L")]
gene.gpf.df <- data.frame(samples=colnames(grouped.log2cpm.gfp),
                          cpm=grouped.log2cpm.gfp["euc",]) #Obp44a
ggplot(data=gene.gpf.df, aes(x=samples, y=cpm)) + geom_bar(stat="identity")

## plot cpm gene expression across samples
reversed <- t(apply(grouped.log2cpm.gfp, 1, function(x) 2^x))
gene.gpf.df <- data.frame(samples = colnames(reversed),
                          cpm=reversed["bnb",])
ggplot(data=gene.gpf.df, aes(x=samples, y=cpm)) + geom_bar(stat="identity")

## plot two cells against each other
library(GGally)
df <- data.frame( c1 = log2cpm[,"fed-gfp-1R"], c2 = log2cpm[,"fed-neg-1R"],
                  c3 = log2cpm[,"fed-gfp-1L"], c4 = log2cpm[,"fed-neg-1L"],
                  c5 = log2cpm[,"fed-gfp-2R"], c6 = log2cpm[,"fed-neg-2R"],
                  c7 = log2cpm[,"fed-gfp-2L"], c8 = log2cpm[,"fed-neg-2L"],
                  c9 = log2cpm[,"sta-gfp-1R"],c10 = log2cpm[,"sta-neg-2R"],
                 c11 = log2cpm[,"sta-gfp-1L"], c12 = log2cpm[,"sta-neg-2L"],
                 c13 = log2cpm[,"sta-gfp-2R"], c14 = log2cpm[,"sta-neg-2R"],
                 c15 = log2cpm[,"sta-gfp-2L"], c16 = log2cpm[,"sta-neg-2L"])

df <- data.frame( fed_1R = log2cpm[,"fed-gfp-1R"], 
                  fed_1L = log2cpm[,"fed-gfp-1L"], 
                  fed_2R = log2cpm[,"fed-gfp-2R"], 
                  fed_2L = log2cpm[,"fed-gfp-2L"], 
                  sta_1R = log2cpm[,"sta-gfp-1R"],
                  sta_1L = log2cpm[,"sta-gfp-1L"],
                  sta_2R = log2cpm[,"sta-gfp-2R"],
                  sta_2L = log2cpm[,"sta-gfp-2L"])

plotList <- list(
  ggplot(data=df, aes(x=fed_1R,y=sta_1R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01), 
  ggplot(data=df, aes(x=fed_1R,y=sta_1L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_1R,y=sta_2R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_1R,y=sta_2L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_1L,y=sta_1R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01), 
  ggplot(data=df, aes(x=fed_1L,y=sta_1L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_1L,y=sta_2R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_1L,y=sta_2L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2R,y=sta_1R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01), 
  ggplot(data=df, aes(x=fed_2R,y=sta_1L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2R,y=sta_2R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2R,y=sta_2L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2L,y=sta_1R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01), 
  ggplot(data=df, aes(x=fed_2L,y=sta_1L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2L,y=sta_2R)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01),
  ggplot(data=df, aes(x=fed_2L,y=sta_2L)) + geom_point(alpha=0.05) + geom_rug(alpha = 0.01)
)


pm <- ggmatrix(plotList,
               4, 4, 
               yAxisLabels = c("fed 1R", "fed 1L", "fed 2R", "fed 2L"), 
               xAxisLabels = c("sta 1R", "sta 1L", "sta 2R", "sta 2L"), 
               byrow=TRUE)
pm

ggplot(data=df, aes(x=fed_1R,y=sta_1L)) + geom_point(alpha=0.05) 
ggplot(data=df, aes(x=fed_1R,y=sta_1L)) + geom_point(alpha=0.05) 
ggplot(data=df, aes(x=fed_1R,y=sta_1L)) + geom_point(alpha=0.05) 

pm <- ggmatrix(
  plotList,
  2, 3,
  c("A", "B", "C"),
  c("D", "E"),
  byrow = TRUE
)
pm

#' for 4 cells:
#' 1 2
#' 1 3
#' 1 4
#' 
library(ggplot2)

x <- data.frame(letters[1:10],abs(rnorm(10)),abs(rnorm(10)),type="x")
y <- data.frame(letters[1:10],abs(rnorm(10)),abs(rnorm(10)),type="y")
# in reality the number of row could be larger than 10 for each x and y

all <- rbind(x,y)
colnames(all) <- c("name","val1","val2","type")
p <- ggplot(all, aes(val1, val2))+ geom_smooth(method = "lm")  + geom_point() +
     facet_grid(~type) 