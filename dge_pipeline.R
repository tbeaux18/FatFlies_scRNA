#!/usr/bin/env Rscript
library("edgeR")

# arg1 = count data
# arg2 = design data
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


# creating DGEList object and filtering rows that do not meet cpm threshold
umi.y <- DGEList(counts=gfp.umi.count, group=umi.group)
keep.umi <- rowSums(cpm(umi.y)>1) >= 2
umi.y <- umi.y[keep.umi, , keep.lib.sizes=FALSE]

# calculating the normalizing factors for each cell
umi.y <- calcNormFactors(umi.y, method="TMM")

# creating experiment design, includes y intercept
design <- model.matrix(~umi.group)

# estimates dispersion factors
umi.y <- estimateDisp(umi.y, design)

# fitting the model based on design
fit <- glmFit(umi.y, design)

# performing the differential expression
lrt <- glmLRT(fit, coef=2)

# isolating top differentially expressed genes Coef 2 vs Coef 1
#temp.lrt <- topTags(lrt, nrow(gfp.umi.count)) %>% as.data.frame

umi.glm_results <- topTags(lrt, n = nrow(gfp.umi.count), sort.by = "none")

# all results
write.csv(umi.glm_results, file='all_diff_exp_results.csv')

# sum(umi.glm_results$table$FDR < .05)

umi.sig <- umi.glm_results[umi.glm_results$table$FDR < .05,]

write.csv(umi.sig, file='sig_diff_exp_results.csv')
