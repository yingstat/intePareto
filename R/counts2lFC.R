#' @title Calculate log2FoldChange
#'
#' @aliases counts2lFC
#'
#' @description
#' \code{counts2lFC} calculates log2FoldChange from counts.
#'
#' @param countData
#' a matrix of counts.
#'
#' @param colData
#' a data.frame with at least a single column. Rows of colData correspond to
#' columns of countData.
#'
#' @param condition
#' the formula expresses how the counts for each gene depend on the variables in
#' colData. The comparisons will be based on the alphabetical order of the
#' levels by default. You can also specify the reference level by ref parameter
#'
#' @param ref
#' specifying the reference level
#'
#' @param type
#' shrinkage estimator, default is "apeglm", the adaptive t prior shrinkage
#' estimator from the 'apeglm' package.
#'
#' @param apeAdapt
#' logical, should apeglm use the MLE estimates of LFC to adapt the prior, or
#' use default.
#'
#' @param ...
#' refer to DESeq2::lfcShrink() for more detailed parameters.
#'
#' @return resLFC
#' a dataframe contains log2FoldChange.
#'
#' # Please note this is a downsampling of the original data.
#' @import DESeq2
#' @importFrom stats relevel
#'
counts2lFC <- function(countData,
                       colData,
                       condition,
                       ref,
                       type = "apeglm",
                       apeAdapt = FALSE,...){
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countData,
                                        colData = colData,
                                        design = ~ condition)
  dds$condition <- stats::relevel(dds$condition, ref=ref)
  dds <- DESeq2::DESeq(dds, fitType="mean")

  resLFC <- DESeq2::lfcShrink(dds,
                              coef = DESeq2::resultsNames(dds)[2],
                              type = type,
                              apeAdapt = apeAdapt,...)
  resLFC <- as.data.frame(resLFC)
  return(resLFC)
}
