#!/usr/bin/env Rscript
library("edgeR")
library("ggplot2")
library("plotly")
library("crosstalk")
library("Glimma")

# arg1 = count data // compbio/fatfly_zumi_test_run.dgecounts.rds
# arg2 = cell data // cell_data.csv
# arg3 = test groups // design.csv
# arg4 = subgroup name // gfp+
args = commandArgs(trailingOnly=TRUE)

# for testing only
args = c('compbio/fatfly_zumi_test_run.dgecounts.rds',
         'cell_data.csv',
         'design.csv',
         'gfp+')

# All count data from zUMIs pipeline () 
AllCounts <- readRDS(args[1])

# cell data supplied from SampleSheet.csv at beginning of pipeline  
cell_data <- read.csv(args[2], header=TRUE, sep=",")

# load design data (designates which samples should be compared)
design <- read.csv(args[3])
test_group <- c(as.integer(stringr::str_split_fixed(design$test, ",", 2)))
ctrl_group <- c(as.integer(stringr::str_split_fixed(design$control, ",", 2)))

# subgroup name for loglog-cpm plot
subgroup <- args[4]

# create all-sample UMI count matrix using cell_data sample_names as column names           // TODO: warn if column names repeat - create a table?
umi.count <- as.matrix(AllCounts$umicount$inex$all)
umi.count <- umi.count[,cell_data$barcode_sequence]
colnames(umi.count) <- cell_data$sample_name

# Remove AllCounts to conserve memory ~ it's a huge file!
rm("AllCounts")

# create UMI matrix that only has data for [~1st] comparison                                // TODO: convert to function and loop through
comparison = c(test_group[1], ctrl_group[1])
included <- cell_data$experiment_group %in% comparison
incl.umi.colnames <- cell_data$sample_name[included]
incl.umi.count <- umi.count[,which(colnames(umi.count) %in% incl.umi.colnames)]

# create 'group' character vector for [~1st] comparison
# each level is a condition
incl.umi.condition <- as.factor(cell_data$condition[included])
conditions <- levels(incl.umi.condition)

# calculate number of replicates per condition (rpc) for [~1st] comparison                // if ncol(incl.umi.count) %% length(comparison) == 0, otherwise throw error
rpc <- (ncol(incl.umi.count) / length(comparison))

# filtering low expressed genes across conditions for [~1st] comparison
# keep genes that have at least 1 cpm in at least greater than or equal to rep number per condition
# unsupervised approach
keep <- rowSums(cpm(incl.umi.count)>2) >= rpc
incl.umi.count <- incl.umi.count[keep, ]

# compute logcpm across entire library for [~1st] comparison
lcpm.incl <- cpm(incl.umi.count, log=TRUE)

# find average/median logcpm across conditions for each gene for [~1st] comparison
lcpm.test <- apply(lcpm.incl[ ,which(incl.umi.condition %in% conditions[1])],1,median)
lcpm.ctrl <- apply(lcpm.incl[ ,which(incl.umi.condition %in% conditions[2])],1,median)

# create interactive plot for median log2cpm values for condition vs. condition [~1st] comparison // TODO: make this a function! ~ paste("| ",conditions[1]," - ",conditions[2]," |",sep="")
plot.df <- data.frame( test = lcpm.test, ctrl = lcpm.ctrl)
plot.df["difference"] <- abs(plot.df$test-plot.df$ctrl)
plot.df <- cbind(plot.df, lcpm.incl)
plot.df["hover"] <- paste("Gene:", rownames(lcpm.incl), sep="")
plot.names <- c(conditions[1], conditions[2], "|difference|", colnames(plot.df[,4:length(colnames(plot.df))]))
plot.df <- SharedData$new(plot.df)

g <- ggplot(data=plot.df, aes(x=ctrl,y=test,text=hover)) + geom_point()
g <- g + labs(title = paste("Median log<sub>2</sub>(cpm) of ",conditions[1]," vs. ", conditions[2]," cells (",subgroup,")",sep=""),
           x = paste("log<sub>2</sub>(cpm) in ",conditions[2]," cells",sep=""),
           y = paste("log<sub>2</sub>(cpm) in ",conditions[1]," cells",sep=""))
gg <- highlight(ggplotly(g, tooltip=c("text")), 
                color="red",
                on = "plotly_selected")
gg <- toWebGL(gg)
ct <- crosstalk::bscols(widths=c(12,12), gg, DT::datatable(plot.df, style="bootstrap", width="100%", colnames=plot.names))
htmltools::save_html(ct, "test.html")

# differential expression

run_edgeR_diff_exp <- function(count_matrix, group, replicate_num) {
  
  # instantiates DGEList object
  y <- DGEList(counts=count_matrix, group=group)
  
  # filtering low expressed genes across conditions
  # keep genes that have at least 1 cpm in at least greater than or equal to rep number per condition
  # unsupervised approach
  keep <- rowSums(cpm(y) > 1) >= replicate_num
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # normalizing scaling factors using modified TMM to handle zero inflation
  y <- calcNormFactors(y, method="TMM")
  
  # creating design matrix based on group
  # includes an intercept
  design <- model.matrix(~group)
  
  # estimating dispersion parameters based on experiment design
  y <- estimateDisp(y, design)
  
  # fitting the model
  fit <- glmFit(y, design)
  
  # performs the differential expression coef 2 vs coef 1
  lrt <- glmLRT(fit, coef=2)
  
  # Glimma interactive visualization of results
  glMDPlot(lrt, counts=count_matrix, groups=group, transform=TRUE)
  
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



run_edgeR_diff_exp(incl.umi.count, incl.umi.condition, rpc)
dev.off()
