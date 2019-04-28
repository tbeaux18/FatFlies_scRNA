#!/usr/bin/env Rscript
library("edgeR")
library("plyr")
library("dplyr")
library("ggplot2")
library("plotly")
library("crosstalk")
library("Glimma")
library("DT")
library("stringr")
library("htmltools")

# arg1 = count data
# arg2 = cell data
# arg3 = design groups to test

args = commandArgs(trailingOnly=TRUE)

run_edgeR_diff_exp <- function(count_matrix, group, replicate_num, basename, qc=FALSE, fdr_thres, output_path) {

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
  # making significance at 0.08 FDR for now, will need to change
  glm_results <- topTags(lrt, n = nrow(count_matrix), sort.by = "none")
  sig_results <- glm_results[glm_results$table$FDR < fdr_thres,]

  if (qc==TRUE){
    all.filename <- paste(basename, "qc_all_diff_exp_results.csv", sep="_")
    all.full.out.path <- paste(output_path, all.filename, sep="/")
    write.csv(glm_results, file=all.full.out.path)
    
    sig.filename <- paste(basename, "qc_sig_diff_exp_results.csv", sep="_")
    sig.full.out.path <- paste(output_path, sig.filename, sep="/")
    write.csv(sig_results, file=sig.full.out.path)
   
    deGenes <- decideTestsDGE(lrt, p=fdr_thres)
    deGenes <- rownames(lrt)[as.logical(deGenes)]
    smear.filename <- paste(basename, "qc_deg_smear.jpeg", sep="_")
    smear.full.out.path <- paste(output_path, smear.filename, sep="/")
    jpeg(file=smear.full.out.path)
    plotSmear(lrt, de.tags=deGenes)
    abline(h=c(-1, 1), col=2)
  } else {
    all.exp.filename <- paste(basename, "exp_all_diff_exp_results.csv", sep="_")
    all.exp.full.out.path <- paste(output_path, all.exp.filename, sep="/")
    
    sig.exp.filename <- paste(basename, "exp_sig_diff_exp_results.csv", sep="_")
    sig.exp.full.out.path <- paste(output_path, sig.exp.filename, sep="/")
    write.csv(glm_results, file=all.exp.full.out.path)
    write.csv(sig_results, file=sig.exp.full.out.path)

    # create smear plot of DEGs at fdr_thres
    deGenes <- decideTestsDGE(lrt, p=fdr_thres)
    deGenes <- rownames(lrt)[as.logical(deGenes)]
    smear.filename <- paste(basename, "exp_deg_smear.jpeg", sep="_")
    smear.full.out.path <- paste(output_path, smear.filename, sep="/")
    jpeg(file=smear.full.out.path)
    plotSmear(lrt, de.tags=deGenes)
    abline(h=c(-1, 1), col=2)
  }


}

qc_diff_exp <- function(umi.count.matrix, cell.data.df, design.data.df, fdr_num, output_path) {

  # iterates down grabbing values from both columns
  for (i in 1:length(design.data.df)) {

    # grabs the comparisons to be made and creates a char vector for of column names to be subsetted
    control.comparison <- as.vector(design.data.df[i,])

    # subsets the cell.data to only the rows in the given control.comparison group
    sub.cell.data <- cell.data.df %>%
                      group_by(experiment_group) %>%
                      filter(experiment_group %in% control.comparison)

    # grabs max replicate count for unsupervised filtering during the DE steps
    replicate.count <- sub.cell.data %>% group_by(gfp_state) %>% tally()
    rep.num <- max(replicate.count$n)

    # makes a vector of col_name to subset the entire umi.count.matrix
    cell.name <- sub.cell.data %>% pull(col_name)

    # convert from tibble to df
    sub.cell.data <- as.data.frame(sub.cell.data)

    # resets the row names to col_name
    rownames(sub.cell.data) <- sub.cell.data$col_name

    # creates a subsetted matrix with new variable name
    qc.name <- levels(factor(sub.cell.data[cell.name,]$condition))[1]

    # subsets the matrix to only columns being compared
    umi.count.matrix.sub <- umi.count.matrix[,cell.name]

    # builds the interactive plots
    build_interactive_plots(umi.count.matrix.sub, cell.data.df, control.comparison, rep.num, "gfp", output_path)

    # for QC, comparing gfp pos vs neg
    comp.group <- factor(sub.cell.data[cell.name,]$gfp_state)
    run_edgeR_diff_exp(umi.count.matrix.sub, comp.group, rep.num, qc.name, qc=TRUE, fdr_thres=fdr_num, output_path)

  }
}



experimental_diff_exp <- function(umi.count.matrix, cell.data.df, testing.groups, fdr.threshold, basefile, output_path) {

  # subsetting data all the way down testing group
  sub.cell.data <- cell.data.df %>%
                    group_by(experiment_group) %>%
                    filter(experiment_group %in% testing.groups)

  # grabbing max replicate count between condition column, may need to fix in the future
  replicate.count <- sub.cell.data %>% group_by(condition) %>% tally()
  rep.num <- max(replicate.count$n)

  # vector of column names to subset the matrix
  exp.names <- sub.cell.data$col_name

  # should use experiment groups, but will keep as condition for now
  exp.group <- sub.cell.data$condition

  # subsetting the matrix
  exp.umi.count <- umi.count.matrix[,exp.names]

  # builds the interactive plots
  build_interactive_plots(exp.umi.count, cell.data.df, testing.groups, rep.num, "gfp+", output_path)

  # runs the differential expression between testing groups
  # only 2 pairwise comparisons
  run_edgeR_diff_exp(exp.umi.count, exp.group, rep.num, basefile, qc=FALSE, fdr.threshold, output_path)
}




build_interactive_plots <- function(sub.umi.count.matrix, cell.data.df, comparison.groups, rep_num, subgroup, output_path) {

  # grabbing column names based on comparison
  cell.col.data <- cell.data.df %>%
                    group_by(experiment_group) %>%
                    filter(experiment_group %in% comparison.groups)


  # factoring combined state column condition.gfp_state
  cell.condition.combined <- levels(factor(cell.col.data$comb_condition))

  # filtering the subset matrix such that at least 4 cells have 1 cpm in an unsupervised method
  # it does not recalculate library sizes due to not being a DGE object
  keep <- rowSums(cpm(sub.umi.count.matrix) > 1) >= rep_num
  sub.umi.count.matrix <- sub.umi.count.matrix[keep, ]

  # compute logcpm across entire library of subset
  lcpm <- cpm(sub.umi.count.matrix, log=TRUE)

  # find average/median logcpm across condition
  lcpm.group1 <- apply(lcpm[,cell.col.data %>%
                              filter(comb_condition==cell.condition.combined[1]) %>%
                              pull(col_name)], 1, median)
  lcpm.group2 <- apply(lcpm[,cell.col.data %>%
                              filter(comb_condition==cell.condition.combined[2]) %>%
                              pull(col_name)], 1, median)

  # create interactive plot for median log2cpm values for condition vs. condition [~1st] comparison // TODO: make this a function! ~ paste("| ",conditions[1]," - ",conditions[2]," |",sep="")
  plot.df <- data.frame( group1 = lcpm.group1, group2 = lcpm.group2)
  plot.df["difference"] <- abs(plot.df$test-plot.df$ctrl)
  plot.df <- cbind(plot.df, lcpm)
  plot.df["hover"] <- paste("Gene:", rownames(lcpm), sep="")
  plot.names <- c(cell.condition.combined[1], cell.condition.combined[2], "|difference|", colnames(plot.df[,4:length(colnames(plot.df))]))

  plot.df <- SharedData$new(plot.df)

  g <- ggplot(data=plot.df, aes(x=group2,y=group1,text=hover)) + geom_point()
  g <- g + labs(title = paste("Median log<sub>2</sub>(cpm) of ",cell.condition.combined[1]," vs. ", cell.condition.combined[2]," cells (",subgroup,")",sep=""),
                x = paste("log<sub>2</sub>(cpm) in ",cell.condition.combined[2]," cells",sep=""),
                y = paste("log<sub>2</sub>(cpm) in ",cell.condition.combined[1]," cells",sep=""))
  gg <- highlight(ggplotly(g, tooltip=c("text")),
                  color="red",
                  on = "plotly_selected")
  gg <- toWebGL(gg)
  ct <- crosstalk::bscols(widths=c(12,12), gg, DT::datatable(plot.df, style="bootstrap", width="100%", colnames=plot.names))

  html.filename <- paste(cell.condition.combined[1], "vs", cell.condition.combined[2], "graph.html", sep="_")
  html.full.out.path <- paste(output_path, html.filename, sep="/")
  htmltools::save_html(ct, html.full.out.path)
}


# All count data from zUMIs pipeline
all.counts <- readRDS(args[1])

# cell data supplied from SampleSheet.csv at beginning of pipeline
cell.data <- read.csv(args[2], header=TRUE, sep=",")
cell.data$col_name <- paste(cell.data$gfp_state, cell.data$condition, cell.data$cell_id, sep=".")
rownames(cell.data) <- cell.data$col_name
cell.data$comb_condition <- paste(cell.data$condition, cell.data$gfp_state, sep=".")
condition.group <- factor(cell.data$condition)
condition.levels <- levels(condition.group)


# load design data (designates which samples should be compared)
design <- read.csv(args[3])
test.group <- c(as.integer(stringr::str_split_fixed(design$test, ",", 2)))
ctrl.group <- c(as.integer(stringr::str_split_fixed(design$control, ",", 2)))
design.df <- data.frame(test.group, ctrl.group)

# create all-sample UMI count matrix using cell_data sample_names as column names           // TODO: warn if column names repeat - create a table?
# umi counts include intron counts (best to use per zUMIs due to accurate UMI collapsing)
umi.count <- as.matrix(all.counts$umicount$inex$all)

# critical that umi.count columns and cell.data rows are in same order
umi.count <- umi.count[,cell.data$barcode_sequence]
colnames(umi.count) <- cell.data[["col_name"]]

# remove initial matrix
rm(all.counts)

# creating working dir and create r_result dir path
current.dir <- getwd()
dir.create("r_results")
full.result.path <- paste(current.dir, "r_results", sep="/")

# run diff exp on control groups
qc_diff_exp(umi.count, cell.data, design.df, 0.08, full.result.path)

# run diff exp on testing groups
experimental_diff_exp(umi.count, cell.data, test.group, 0.08, "SvsF", full.result.path)


dev.off()
