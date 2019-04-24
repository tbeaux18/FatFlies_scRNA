#!/usr/bin/env Rscript
library("edgeR")

# arg1 = count data
# arg2 = cell data
# arg3 = test groups
args = commandArgs(trailingOnly=TRUE)

# All count data from zUMIs pipeline
AllCounts <- readRDS(args[1])

# supplied from SampleSheet.csv at beginning of pipeline
design.data <- read.csv(args[2], header=TRUE, sep=",")

# sorting the design matrix lexographically by cell barcode sequence 
sorted.design <- design.data[order(design.data[,"barcode_sequence"]),]

# may not need this column
sorted.design$cat_group <- paste(sorted.design$gfp_state,
                                 sorted.design$condition,
                                 sep=".")

sorted.design$col_name <- paste(sorted.design$gfp_state,
                                sorted.design$condition,
                                sorted.design$cell_id,
                                sep=".")



# UMI count matrix columns are in lexographic order of cell barcodes
umi.count.inex <- as.matrix(AllCounts$umicount$inex$all)
colnames(umi.count.inex) <- sorted.design[["col_name"]]



# grabbing the columns that are solely gfp positive to be compared
# need to figure out how to grab this data without hardcoded column names
gfp.umi.count <- umi.count.inex[,c("G.fed.1", "G.fed.5", "G.fed.2", "G.fed.6", "G.starved.13", "G.starved.10", "G.starved.14", "G.starved.9")]
# create grouping based on design data
umi.group <- substr(colnames(gfp.umi.count), 3, 3)




run_edgeR_diff_exp <- function(count_matrix, group, replicate_num) {
  
  # instantiates DGEList object
  y <- DGEList(counts=count_matrix, group=group)
  
  # filtering low expressed genes across conditions
  # keep genes that have at least 1 cpm in at least greater than or equal to rep number per condition
  # unsupervised approach
  keep <- rowSums(cpm(y) > 1) >= replicate_num
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalizing scaling factors using modified TMM to handle zero inflation
  y <- calcNormFactors(y, method="TMMwzp")
  
  # creating design matrix based on group
  # includes an intercept
  design <- model.matrix(~group)
  
  # estimating dispersion parameters based on experiment design
  y <- estimateDisp(y, design)
  
  # fitting the model
  fit <- glmFit(y, design)
  
  # performs the differential expression coef 2 vs coef 1
  lrt <- glmLRT(fit, coef=2)
  
  # generating a table of all results
  glm_results <- topTags(lrt, n = nrow(count_matrix), sort.by = "none")
  write.csv(glm_results, file='all_diff_exp_results.csv')
  
  # filtering results based on FDR significance
  sig_results <- glm_results[glm_results$table$FDR < .05,]
  write.csv(sig_results, file='sig_diff_exp_results.csv')
  
  # create smear plot of DEGs at 0.05 threshold
  deGenes <- decideTestsDGE(lrt, p=0.05)
  deGenes <- rownames(lrt)[as.logical(deGenes)]
  jpeg(file="deg_smear.jpeg")
  plotSmear(lrt, de.tags=deGenes)
  abline(h=c(-1, 1), col=2)
}



run_edgeR_diff_exp(gfp.umi.count, umi.group, 4)
dev.off()