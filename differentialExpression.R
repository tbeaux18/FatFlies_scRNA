#' edgeR script
library(edgeR)

# Load UMI count data into dataframe (samples, groups, counts)
AllCounts <- readRDS("compbio/fatfly_zumi_test_run.dgecounts.rds")
umi <- as.matrix(AllCounts$umicount$exon$all) 

# Load sample & design data
samples.df <- read.csv("samples.csv")
design.df <- read.csv("design.csv")
test_group <- c(as.integer(stringr::str_split_fixed(design.df$test, ",", 2)))
ctrl_group <- c(as.integer(stringr::str_split_fixed(design.df$ctrl, ",", 2)))

# Re-order count matrix using sample sheet
umi <- umi[,samples.df$barcode_sequence]
colnames(umi) <- samples.df$sample_name

# Extract group names from names 
group <- stringr::str_split_fixed(colnames(umi), '_', 3)[,1]

# Build DGEList object for edgeR
y <- DGEList(counts=umi, group=group, samples=samples.df)

# remove sta.neg.1L ???

## Compute cpm, lcpm across samples (edgeR)
cpm <- cpm(umi)
lcpm <- cpm(umi, log=TRUE)



## Run and plot all sample PCA w/ scree plot (GLIMMA)
lcpm.pca <- prcomp(t(cpm))
lcpm.pca <- summary(lcpm.pca)
lcpm.pca.df <- cbind(samples.df, lcpm.pca$x)

## plot PCA variance explained across all samples (scree plot)
lcpm.pca.var <- lcpm.pca$sdev^2 
lcpm.pca.var.per <- round(lcpm.pca.var/sum(lcpm.pca.var)*100, 1)
lcpm.pc_names <- factor(colnames(lcpm.pca$x), levels = colnames(lcpm.pca$x))
df <- data.frame(lcpm.pca.var.per, lcpm.pc_names)
df["text"] <- paste("Var exp: ", df$lcpm.pca.var.per, sep="")
lcpm.pca.all.g <- ggplot(data=df, aes(x=lcpm.pc_names, y=lcpm.pca.var.per, text=text)) + geom_bar(stat="identity")
lcpm.pca.all.scree <- ggplotly(lcpm.pca.all.g, tooltip = c("text"))
lcpm.pca.all.scree

## plot PC1 vs. PC2 across all samples
## plot PC2 vs. PC1 (color by sta/fed)
lcpm.pca.all<-ggplot(data=lcpm.pca.df, aes(x=PC1, y=PC2, label=sample_name, color=experiment_group)) +
                      geom_text() +
                      xlab(paste("PC1 - ", lcpm.pca.var.per[1], "%", sep="")) +
                      ylab(paste("PC2 - ", lcpm.pca.var.per[2], "%", sep="")) +
                      ggtitle("PC1 vs. PC2 (fed vs. starved)")
ggplotly(lcpm.pca.all)

## Find genes expressed in 90% of all samples, and plot percent of common genes detected (ggplot)

## Run and plot test/ctrl PCA w/ scree plot (GLIMMA)

## Plot cpm & log2cpm experimental vs. control (interactive faceted graph) ***

## Compute mean cpm/lcpm across experimental groups and plot experimental vs control (same as above)

## Add ENTREZ IDs (org.Dm.eg.db)

## Filter genes of interest and visualize filtered genes (group-wise lcpm before vs. after)

## Normalize genes using ??? visualize normalization with PCA & RLE plots (PCA of lcpm before vs. after)

## Fit error model of gene expression (edgeR)

## Perform differential expression (edgeR)

## Plot logFC with volcano plot (GLIMMA)

## Plot heatmap of top DE genes (heatmap)

## Future goal: Perform gene-set & pathway enrichment analyses

## Future goal: Build in multiple experiment tests


# Get restricted sample names (put somewhere else?)
rsample_names = as.vector(samples.df$sample_name[samples.df$experiment_group == test_group[1] | 
                                                   samples.df$experiment_group == ctrl_group[1]])


## Run and plot test/ctrl group PCA w/ scree plot
#####
# Filter genes without >1 log2cpm in more than one cell
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# Grab cpm and log(cpm) (also mean and median)
cpm <- cpm(umi)
lcpm <- cpm(umi, log=TRUE)

# Perform normalization
y <- calcNormFactors(y)

# Model dispersions of each gene
y <- estimateDisp(y)

# Perform exact test to identify DE genes
et <- exactTest(y)
topTags(et)

# 
(x/((sum(x)*1e6))) + 2
log2((x/((sum(x)*1e6))) + 28)
log2cpm <- t(apply(umi,1,function(x) ""))
log2cpm <- t(apply(umi,1,function(x) log2((x/((sum(x)*1e6))) + 28)))
log2cpm <- t(apply(umi, 1, function(x) log2(((x)/((sum(x)*1e6)))))); hist(log2cpm)


# Filter genes
