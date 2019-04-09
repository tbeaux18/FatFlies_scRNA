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
#' Go.db & org.Dm.eg.db packages needed for GO analysis

library(RUVSeq)
library(EDASeq)
library(edgeR)

###### LOAD UMI COUNT DATA (needs to be automated) ######

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

############################################

###### Filtering lowly expressed genes  ######

library(ggplot2)
pdf("plot3.pdf")
count <- data.frame(data.inex)
cpm_matr <- data.frame(cpm(count))

df <- rbind(data.frame(count=count$fed.gfp.1R, sample="fed.gfp.1R", cpm=cpm_matr$fed.gfp.1R), 
            data.frame(count=count$fed.gfp.1L, sample="fed.gfp.1L", cpm=cpm_matr$fed.gfp.1L),
            data.frame(count=count$fed.gfp.2R, sample="fed.gfp.2R", cpm=cpm_matr$fed.gfp.2R),
            data.frame(count=count$fed.gfp.2L, sample="fed.gfp.2L", cpm=cpm_matr$fed.gfp.2L),
            
            data.frame(count=count$fed.neg.1R, sample="fed.neg.1R", cpm=cpm_matr$fed.neg.1R),
            data.frame(count=count$fed.neg.1L, sample="fed.neg.1L", cpm=cpm_matr$fed.neg.1L),
            data.frame(count=count$fed.neg.2R, sample="fed.neg.2R", cpm=cpm_matr$fed.neg.2R),
            data.frame(count=count$fed.neg.2L, sample="fed.neg.2L", cpm=cpm_matr$fed.neg.2L),
            
            data.frame(count=count$sta.gfp.1R, sample="sta-gfp-2L", cpm=cpm_matr$sta.gfp.1R),
            data.frame(count=count$sta.gfp.1L, sample="sta.gfp.1L", cpm=cpm_matr$sta.gfp.1L),
            data.frame(count=count$sta.gfp.2R, sample="sta.gfp.2R", cpm=cpm_matr$sta.gfp.2R),
            data.frame(count=count$sta.gfp.2L, sample="sta.gfp.2L", cpm=cpm_matr$sta.gfp.2L),
            
            data.frame(count=count$sta.neg.1R, sample="sta.neg.1R", cpm=cpm_matr$sta.neg.1R),
            data.frame(count=count$sta.neg.1L, sample="sta.neg.1L", cpm=cpm_matr$sta.neg.1L),
            data.frame(count=count$sta.neg.2R, sample="sta.neg.2R", cpm=cpm_matr$sta.neg.2R),
            data.frame(count=count$sta.neg.2L, sample="sta.neg.2L", cpm=cpm_matr$sta.neg.2R))

ggplot(df, aes(x=df$count)) + geom_histogram(binwidth=5)   + xlim(c(-2,500)) +  ylim(c(-2,3000)) + facet_wrap(~sample)
ggplot(df, aes(x=df$cpm)) + geom_histogram(binwidth=100) + xlim(c(-2,2500)) + ylim(c(-2,2500)) + facet_wrap(~sample)
ggplot(df, aes(x=df$count, y=df$cpm)) + geom_point() + xlim(c(-2,4200)) + ylim(c(-2,18000)) + facet_wrap(~sample) + geom_hline(yintercept = 2100)
ggplot(df, aes(x=df$count, y=df$cpm)) + geom_point(alpha=0.2) + geom_hline(yintercept=9000, color="steelblue", alpha=0.7) + geom_vline(xintercept=2100, color="red", alpha=0.7) + facet_wrap(~sample) 

## a gene must have > 5 reads >2 cells)
filter <- apply(data.inex,1,function(x) length(x[x>5])>=2)
filtered <- data.inex[filter,] ## num of genes eliminated:nrow(data.inex) - nrow(filtered)

## get gene names & spike names
spikes <- filtered[grepl("^ERCC",rownames(filtered)),0]; spikes <- c(row.names(spikes))
genes <- filtered[!grepl("^ERCC",rownames(filtered)),0]; genes <- c(row.names(genes))

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

###### UMI Normalization  ######

## plot relative log-expression across libraries
#pdf("plots.pdf")
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

## upper-quantile normaliation (EDAseq)
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

## estimate factors of unwanted variation 
set1 <- RUVg(set, spikes, k=2)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
pLotPCA(set1, col=colors[x], cex=1.2)
pData(set1) 

## DE analysis ~ model w/ negative binomial GLM  (edgeR)
design <- model.matrix(~x +W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
plotMD(lrt)
topTags(lrt)

results <- lrt$table
results <- results[order(lrt$table$PValue),]

################################# RUVSeq/edgeR Vignette Code #################################

results <- lrt$table
results <- results[order(lrt$table$PValue),]

###### Gene Ontology ####### 

## download entrez_ids (ignore rows w/ no compatible EID for now)
require(org.Dm.eg.db)
alias2eg <- as.list(org.Dm.egALIAS2EG)
results <- results[rownames(results) %in% names(alias2eg), ]

## add entrez_ids to table
named_entrez_ids <- alias2eg[rownames(results)]; 
rm(entrez_ids); entrez_ids = character()
for (i in 1:length(named_entrez_ids)){
  entrez_ids <- c(entrez_ids, named_entrez_ids[[i]][1])
  }; results <- cbind(results, entrez_ids)

head(results)

## Gene Ontology Analysis
#go <- goana(lrt, species="Dm")
#topGO(go, sort="up")


###########################


dev.off()
################################# END OF SCRIPT #################################


###### Backburner ideas ###### 

#diffExp <- function(x="um?"){
  ## do something
  
  ## return differentially expressed genes
#  genes <- c("gene1","gene2","gene3")
#  log_effect <- c(0.08,0.07,0.06)
#  diffExpGenes <- data.frame(genes, log_effect)
#  return(diffExpGenes)
#}

##//diffExp()