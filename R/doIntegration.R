#' @title Do integrative analysis
#'
#' @aliases doIntegration
#'
#' @description
#' \code{doIntegration} calculate log2FoldChange of RNA-Seq and ChIP-Seq and
#' then calculate Z scores for each marker.
#'
#' @param res
#' a list result from \code{doMatch} function.
#'
#' @param ref
#' specifying the reference level
#'
#' @param type
#' shrinkage estimator, default is "apeglm", the adaptive t prior shrinkage
#' estimator from the 'apeglm' package.
#'
#' @param apeAdapt
#' logical, should apeglm use the MLE estimates of LFC to adapt the prior,
#' or use default.
#'
#' @return df_final
#' a dataframe contains log2FoldChange of RNA-Seq and ChIP-Seq and Z scores
#' for each marker.
#'
#' @import DESeq2
#' @importFrom stats na.omit
#' @examples
#' data(res)
#' \donttest{doIntegration(res = res,ref="wild.type")
#' }
#'
#' @export
#'
doIntegration <- function(res,
                          ref,
                          type = "apeglm",
                          apeAdapt = FALSE){
  checkInputParam <- function(res, ref, type, apeAdapt) {

    # simple input param checks
    if(missing(res) || is.null(res) ||
       missing(ref) || is.null(ref) ||
       missing(type) || is.null(type) ||
       missing(apeAdapt) || is.null(apeAdapt)) {
      stop("arguments res, ref, type, apeAdapt must be specified")
    }


    # simple input param checks
    if(is.null(res)) {
      stop("res not specified")
    }

    if(is.list(res) == FALSE) {
      stop("res must be a list")
    }
    if(is.null(ref)) {
      stop("ref not specified")
    }
    if(is.character(ref) == FALSE) {
      stop("ref must be a character indicating the reference level of the
           condition")
    }
    if(!c("res.rna")%in%names(res)){
      stop('res must contain a dataframe named res.rna')
    }

    if(!c("res.chip")%in%names(res)){
      stop('res must contain a dataframe named res.chip')
    }

    if(!c("matched.data")%in%names(res)){
      stop('res must contain a dataframe named matched.data')
    }
  }

  checkInputParam(res = res,
                  ref = ref,
                  type = type,
                  apeAdapt = apeAdapt)

  message('use DESeq2 to calculate log2FoldChange of RNA-Seq')
  rnaseq.counts <- res$res.rna[,-1]
  rownames(rnaseq.counts) <- res$res.rna$external_gene_name
  condition <- factor(do.call(rbind, strsplit(x = colnames(rnaseq.counts),
                                              split = '_REP'))[,1])

  colData <- data.frame(condition = condition,
                        type = rep("RNA",length(condition)))
  rownames(colData) <- colnames(rnaseq.counts)

  rnalFC <- counts2lFC(countData = rnaseq.counts,
                       colData = colData,
                       condition = condition,
                       ref = ref,
                       type = type,
                       apeAdapt = apeAdapt)
  colnames(rnalFC) <- paste0("RNAseq.", colnames(rnalFC))

  message('use DESeq2 to calculate log2FoldChange of ChIP-Seq')

  chipseq.counts <- res$res.chip[,-1]
  rownames(chipseq.counts) <- res$res.chip$external_gene_name
  markers <- unique(do.call(rbind, strsplit(x = colnames(chipseq.counts),
                                            split = '_HM_'))[,1])
  chiplFC <- list()
  for (i in markers){
    id <- grep(i, colnames(chipseq.counts))
    c.names <- grep(i, colnames(chipseq.counts), value = TRUE)
    c.replicates <- do.call(rbind, strsplit(x = c.names, split = '_HM_'))[,2]
    c.cond <- do.call(rbind, strsplit(x = c.replicates, split = '_REP'))[,1]
    counts.data <- chipseq.counts[,id]
    condition <- factor(c.cond)
    colData <- data.frame(condition = condition,
                          type= rep("chip",length(condition)))
    rownames(colData) <- colnames(counts.data)
    chiplFC[[i]] <- counts2lFC(countData = counts.data,
                               colData = colData,
                               condition= condition,
                               ref = ref,
                               type = type,
                               apeAdapt = apeAdapt)

  }

  chiplFC.df <- do.call("cbind",chiplFC)

  rnaseq.chipseq.res <- merge(rnalFC, chiplFC.df, by=0)

  message("Z score calculation")

  id <- grep(".log2FoldChange", colnames(rnaseq.chipseq.res))
  rownames(rnaseq.chipseq.res) <- rnaseq.chipseq.res$Row.names
  rnaseq.chipseq.res <- rnaseq.chipseq.res[,id]
  rnaseq.chipseq.res <- stats::na.omit(rnaseq.chipseq.res)
  df_final <- rnaseq.chipseq.res

  for (i in markers){
    df_final[[paste0("z.",i)]] <-
    (df_final$RNAseq.log2FoldChange/stats::sd(
      df_final$RNAseq.log2FoldChange))*
      (df_final[[paste0(i,".log2FoldChange")]]/
         stats::sd(df_final[[paste0(i,".log2FoldChange")]]))
  }

  return(df_final)

}
