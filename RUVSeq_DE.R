library(RUVSeq)
library(zebrafishRNASeq)
library(EDASeq)
library(edgeR)

## load and expect zfGenes df
data(zfGenes)
head(zfGenes)
tail(zfGenes)

## filter out unexpressed genes ~ obviously need to modify for celseq2 UMI
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered <- zfGenes[filter,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

## store data in S4 class
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet(as.matrix(filtered), 
                           phenoData = data.frame(x, row.names=colnames(filtered)))

## upper-quantile normaliation (EDAseq)
set <- betweenLaneNormalization(set, which="upper")

## estimate factors of unwanted variation 
set1 <- RUVg(set, spikes, k=1)
pData(set1) 

## DE analysis using negative binomial GLM approach (edgeR)
## TODO: Iterate through list of test & control groups and re-run differential expression
design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
