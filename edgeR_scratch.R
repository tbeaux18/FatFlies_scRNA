#' edgeR script

library(edgeR)

# Load UMI count data into dataframe (samples, groups, counts)
AllCounts <- readRDS("compbio/fatfly_zumi_test_run.dgecounts.rds")
umi_all <- as.matrix(AllCounts$umicount$exon$all) 

# Name columns of UMI table
samples <- c("sta-gfp-2R", "fed-neg-2L", "sta-neg-2L", "fed-gfp-1R",
             "fed-neg-1R", "fed-gfp-1L", "sta-gfp-1L", "sta-gfp-2L", 
             "sta-neg-1L", "fed-gfp-2R", "fed-neg-2R", "fed-gfp-2L", 
             "fed-neg-1L", "sta-neg-1R", "sta-neg-2R", "sta-gfp-1R")
colnames(umi_all) <- samples

# Restrict to samples of interest
umi <- umi_all[,grepl("gfp",colnames(umi_all))] 
umi <- umi_all[,sort(colnames(umi))]

# Extract group names from names
group <- stringr::str_split_fixed(colnames(umi), '-', 3)[,1]

# Build DGEList object for edgeR
y <- DGEList(counts=umi, group=group)

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
