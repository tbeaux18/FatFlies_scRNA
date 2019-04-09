#' Compare expression between a test group and a control group
#' 
#' TODO: FUNCTIONALITY / VALIDITY CHECKS
#' automate umicount.rds file import (w) 
#' automate experiment group import
#' automate test_group/ctrl_group
#' convert this data into test/ctrl SingleCellSequencing object
#' check that column names match cell bar codes
#' check that experiment groups match cell bar codes
#' check that test groups match cell bar codes
#'      
#' library(zebrafishRNASeq)
#' data(zfGenes)
#' head(zfGenes)
#' tail(zfGenes)
#' Go.db & org.Dm.eg.db packages needed for GO analysis

library(RUVSeq)
library(EDASeq)
library(edgeR)


## Load experiment group data (fed-gfp=1,fed-neg=2, sta-gfp=3, sta-neg=4)
sample_names <- c("sta-gfp-2R", "fed-neg-2L", "sta-neg-2L", "fed-gfp-1R",
                  "fed-neg-1R", "fed-gfp-1L", "sta-gfp-1L", "sta-gfp-2L", 
                  "sta-neg-1L", "fed-gfp-2R", "fed-neg-2R", "fed-gfp-2L", 
                  "fed-neg-1L", "sta-neg-1R", "sta-neg-2R", "sta-gfp-1R")
experiment_groups <- c(3,2,4,1,2,1,3,3,4,1,2,1,2,4,4,3)
experiment <- data.frame(experiment_groups, sample_names)
experiment

## Load comparison data
test_group <- c(1,1,3,3,4) # fed-gfp+ // fed-gfp+ // sta-gfp+ // sta-gfp+ // sta-gfp- 
ctrl_group <- c(3,2,4,1,2) # sta-gfp+ // fed-gfp- // sta-gfp- // fed-gfp+ // fed-gfp-
comparisons <- data.frame(test_group, ctrl_group)
comparisons

## Load UMI count data into dataframe
AllCounts <- readRDS("compbio/fatfly_zumi_test_run.dgecounts.rds")
data.inex <- as.matrix(AllCounts$umicount$inex$all) 
colnames(data.inex) <- sample_names

## filter out lowly expressed genes (must have more than 5 reads in at least two cells)
filter <- apply(data.inex,1,function(x) length(x[x>5])>=2)
filtered <- data.inex[filter,] ## num of genes eliminated:nrow(data.inex) - nrow(filtered)

## get gene names & spike names
spikes <- filtered[grepl("^ERCC",rownames(filtered)),0]; spikes <- c(row.names(spikes)); length(spikes)
genes <- filtered[!grepl("^ERCC",rownames(filtered)),0]; genes <- c(row.names(genes)); length(genes)

## automate construction of group vector
## re-organize dataframe by group
column_order <- c("fed-gfp-1R","fed-gfp-1L","fed-gfp-2R","fed-gfp-2L",
                  "fed-neg-1R","fed-neg-1L","fed-neg-2R","fed-neg-2L",
                  "sta-gfp-1R","sta-gfp-1L","sta-gfp-2R","sta-gfp-2L",
                  "sta-neg-1R","sta-neg-1L","sta-neg-2R","sta-neg-2L")
filtered <- filtered[,column_order]

## remove samples that are not being compared
samples <- c("fed-gfp-1R","fed-gfp-1L","fed-gfp-2R","fed-gfp-2L",
              "sta-gfp-1R","sta-gfp-1L","sta-gfp-2R","sta-gfp-2L")
filtered <- filtered[,samples]

## create EDASeq object using only *compared samples*
x <- as.factor(rep(c("fed-gfp", "sta-gfp"), each=4))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))

## plot relative log-expression across libraries
#pdf("sample_plots.pdf")
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(set, col=colors[x], cex=1.2)
#dev.off()

## upper-quantile normaliation (EDAseq)
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(set, col=colors[x], cex=1.2)

## estimate factors of unwanted variation 
set1 <- RUVg(set, spikes, k=1)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
pData(set1) 

## DE analysis using negative binomial GLM approach (edgeR)
design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
tbl <- lrt$table
tbl <- tbl[order(tbl$PValue),]


################################# RUVSeq/edgeR Vignette Code #################################

## Gene Ontology
require(org.Dm.eg.db)
alias <- as.list(org.Dm.egALIAS2EG)
gene_ids <- c(rownames(tbl)); length(gene_ids)
ls <- alias[gene_ids]

## remove any rows that don't have aliases
tbl <- tbl[rownames(tbl) %in% names(alias), ]; length(tbl)
entrez_ids <- unlist(alias[rownames(tbl)], use.names=F); length(entrez_ids)
rownames(lrt$table) <- entrez_ids

# Remove pathway identifiers that do not map to any entrez gene id

## Gene Ontology Analysis
go <- goana(lrt, species="Dm")
topGO(go, sort="up")


################################# END OF SCRIPT #################################


diffExp <- function(x="um?"){
  ## do something
  
  ## return differentially expressed genes
  genes <- c("gene1","gene2","gene3")
  log_effect <- c(0.08,0.07,0.06)
  diffExpGenes <- data.frame(genes, log_effect)
  return(diffExpGenes)
}

##//diffExp()