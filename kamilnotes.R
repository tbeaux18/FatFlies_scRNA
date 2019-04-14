


## commands run with kamil
#' hist(x[,1])
#' hist(x/1e6)
#' hist((x+1)/1e6)
#' hist(log2((x+1)/1e6))
#' log2cpm <- apply(data.inex, 1, function(x) log2((x+1)/1e6))
#' dim(log2cpm)
#' log2cpm <- t(apply(data.inex, 1, function(x) log2((x+1)/1e6)))
#' dim(log2cpm)
#' plotPCA(log2cpm)
#' log2cpm(1:5,1:5)
#' log2cpm <- apply(data.inex, 1, function(x) log2(((x+1)/sum(x))*1e6))
#' hist(log2cpm[,1])
#' log2cpm <- t(apply(data.inex, 1, function(x) log2(((x+1)/sum(x))*1e6))
#' colSums(data.inex)
#' prcomp(t(log2cpm))
#' pca <- prcomp(t(log2cpm))
#' names(pca)
#' pca = summary(pca)
#' names(pca)
#' pca$importance
#' pca$x
#' stringr::str_split_fixed(colnames(log2cpm), '-', 3)
#' meta <- data.frame(stringr::str_split_fixed(colnames(log2cpm), '-', 3))
#' meta <- cbind(meta,pca$x)
#' head(meta)
#' mod <- lm(PC1 ~ X1+X2+X3, meta)
#' summary(mod)
#' mod <- lm(PC2 ~ X1+X2+X3, meta)
#' summary(mod)
#' pca$rotation[1:5,]
#' sorted_PC1 <- sort(pca$rotation[,1])
#' head(sorted_PC1)
#' tail(sorted_PC1)
#' plot(log2cpm["erm",])
#' plot(log2cpm["erm",], meta$PC1)