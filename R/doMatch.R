#' @title Match the RNA-Seq and ChIP-Seq data on the gene level
#'
#' @aliases doMatch
#'
#' @description
#' \code{doMatch} computes the number of reads (counts) fall into specific
#' genomic region such as promoter or genebody for ChIP-Seq, and calculate the
#' gene expression in counts, and then match the RNA-Seq and ChIP-Seq data on
#' the gene level with the method of "weighetd mean" or "highest".
#'
#' @param rnaMeta
#' metadata for RNA-Seq include column named "condition" indicates the
#' experiment condition or cell type, and column named "files" indicates the
#' paths of cprresponing abundance.tsv file that is returned from Kallisto.
#'
#' @param chipMeta
#' metadata for ChIP-Seq include column of "mark" column indicates the markers
#' of histone modifications, column of "condition" indicates the experiment
#' condition or cell type, and "files" column indicates the paths and the file
#' names of the aligned bam files.
#'
#' @param region
#' region has to be specified as "promoter" or "genebody".
#'
#' @param method
#' method has to be specified as "weighted.mean" or"highest" if region is set as
#' "promoter".
#'
#' @param ensemblDataset
#' Ensembl Dataset you want to use. To see the different datasets available
#' within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by
#' listDatasets(mart).
#'
#' @param host
#' specify the archived versions of Ensembl. To see the available archived
#' versions do: biomaRt::listEnsemblArchives()
#'
#' @param fragLength
#' extend reads toward the 3'-end to the average DNA fragment size obtained
#' after DNA size selection.
#'
#' @param promoter.length
#' the length of the promoter region.
#'
#' @return A list with the following three items.
#'
#' @return res.rna
#' a data frame contains RNA-Seq counts
#'
#' @return res.chip
#' a data frame contains ChIP-Seq counts
#'
#' @return matched.data
#' a dataframe contains matched RNA-Seq counts and ChIP-Seq counts.
#'
#' @importFrom GenomicRanges GRanges resize seqnames
#' @importFrom stats aggregate
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments readGAlignments
#' @import utils
#' @importFrom biomaRt useMart getBM
#' @examples
#' data(test_rna_meta)
#' data(test_chip_meta)
#' \donttest{
#' for(i in test_rna_meta$SRR){
#' test_rna_meta$files <- system.file("extdata",paste0(i,".tsv"),
#' package = "intePareto")
#' }
#' for(i in test_chip_meta$SRR){
#' test_chip_meta$files <- system.file("extdata", paste0(i,".bam"),
#' package = "intePareto")
#' }
#' doMatch(rnaMeta = test_rna_meta,
#' chipMeta = test_chip_meta,
#' region = "promoter",
#' method = "weighted.mean",
#' host = "http://aug2017.archive.ensembl.org",
#' ensemblDataset = "mmusculus_gene_ensembl")
#' }
#' @export
#'
doMatch <- function(rnaMeta,
                    chipMeta,
                    region,
                    method,
                    ensemblDataset,
                    host,
                    fragLength = 180,
                    promoter.length=5000){
  ### check input
  checkInputParam <- function(rnaMeta,
                              chipMeta,
                              region,
                              method,
                              ensemblDataset,
                              host) {

    # simple input param checks
    if(missing(rnaMeta) || is.null(rnaMeta) ||
       missing(chipMeta) || is.null(chipMeta) ||
       missing(region) || is.null(region) ||
       missing(method) || is.null(method) ||
       missing(host) || is.null(host) ||
       missing(ensemblDataset) || is.null(ensemblDataset)) {
      stop("arguments rnaMeta, chipMeta, region,
           method, host, ensemblDataset must be specified")
    }

    # check rnaMeat
    if(class(rnaMeta)!="data.frame"){
      stop("rnaMeta must be a data.frame")
    }

    if(!c("condition")%in%colnames(rnaMeta)){
      stop('rnaMeta has to contain "condition" column that indicates
           the condition of experiment.')
    }

    if(!c("files")%in%colnames(rnaMeta)){
      stop('rnaMeta has to contain "files" column indicates the paths of
           cprresponing abundance.tsv file calculated from Kallisto')
    }

    if(nrow(rnaMeta) <= 1) {
      stop('rnaMeta has to contain at least 2 rows
           (one for each biological condition)')
    }


    # check chipMeta
    if(class(chipMeta)!="data.frame"){
      stop("chipMeta must be a data.frame")
    }
    # check dataset
    mart = biomaRt::useMart("ensembl")
    martlist <- biomaRt::listDatasets(mart)
    if(!ensemblDataset%in%martlist$dataset){
      stop("ensemblDataset must be a dataset available within a biomaRt
           you can e.g. do: mart = useMart('ensembl'),
           followed by listDatasets(mart).")
    }

    if(!c("mark")%in%colnames(chipMeta)){
      stop('chipMeta has to contain "mark" column indicates the markers of
           histone modifications')
    }

    if(!c("condition")%in%colnames(chipMeta)){
      stop('chipMeta has to contain "condition" column that indicates the
           condition of experiment')
    }

    if(!c("files")%in%colnames(chipMeta)){
      stop('chipMeta has to contain "files" column indicates the paths of
           aligned bam files')
    }

    if(!is.character(chipMeta$files)){
      stop('the class of the paths of the aligned bam files should be
           character')
    }

    if(nrow(chipMeta) <= 1) {
      stop('chipMeta has to contain at least 2 rows
           (one for each biological condition)')
    }


    # check conditions in chipMeta and rnaMeta
    if(length(unique(chipMeta$condition)) != 2 |
       length(unique(rnaMeta$condition)) != 2) {
      stop('two distinct biological conditions
           expected in chipMeta and rnaMeta')
    }
    if(length(unique(c(chipMeta$condition, rnaMeta$condition))) != 2) {
      stop('identical biological conditions
           expected in chipMeta and rnaMeta')
    }


    # check region and method
    if(length(region) != 1) {
      stop('region has to be specified as "promoter" or "genebody"')
    }

    if(!is.character(region)) {
      stop('region has to be specified as "promoter" or "genebody"')
    }

    if(!region%in%c("promoter", "genebody")){
      stop('region has to be specified as "promoter" or "genebody"')
    }

    if(region=="promoter"){

      if(length(method) != 1) {
        stop('method has to be specified as "weighted.mean" or
             "highest" if region is set as "promoter"')
      }

      if(!is.character(method)) {
        stop('method has to be specified as "weighted.mean" or
             "highest" if region is set as "promoter"')
      }

      if(!method%in%c("weighted.mean","highest")){
        stop('method has to be specified as "weighted.mean" or
             "highest" if region is set as "promoter"')
      }
    }
  }

  checkInputParam(rnaMeta = rnaMeta, chipMeta = chipMeta,
                  region = region, method = method, host = host,
                  ensemblDataset = ensemblDataset)

  message("get RNA counts")
  ### get RNA counts
  est_counts <- function(x) {
    ex <- utils::read.table(x,sep = "\t",header = TRUE,row.names = "target_id")
    ex[,"est_counts",drop=FALSE]
  }

  files <- NULL
  for (i in unique(rnaMeta$condition)){
    a <- rnaMeta[rnaMeta$condition==i,]$files
    names(a) <- paste(i,"_REP",1:length(a), sep = "")
    files <- c(a, files)
  }
  fileNames <- names(files)

  exprsList <- list()
  for (f in fileNames) {
    x <- est_counts(files[f])
    names(x) <- f
    exprsList[[f]] <- x
  }
  exprs <- do.call("cbind",exprsList)

  # annotation form biomaRt package
  ensembl <- biomaRt::useMart("ensembl", dataset = ensemblDataset, host = host)

  N <- data.frame(biomaRt::getBM(attributes=c("ensembl_transcript_id_version",
                                              "ensembl_gene_id",
                                              "external_gene_name",
                                              "chromosome_name",
                                              "start_position",
                                              "end_position",
                                              "transcript_start",
                                              "transcript_end",
                                              "transcription_start_site",
                                              "transcript_length",
                                              "strand"),
                                 filters=c('ensembl_transcript_id_version'),
                                 values = rownames(exprs),
                                 mart=ensembl))

  exprsRNA <- merge(N,exprs, by.x = "ensembl_transcript_id_version", by.y=0)
  # remove it if it causes problem
  exprsRNA <- exprsRNA[exprsRNA$chromosome_name %in% c(1:100,"X","Y"),]
  aln <- GenomicAlignments::readGAlignments(file = chipMeta$files[1])
  #### get chip counts
  message("get chip counts")
  if(region=="promoter"){
    # transcript starting site
      if("chr1"%in%as.character(GenomicRanges::seqnames(aln))){
      tss <- GenomicRanges::GRanges(
        seqnames = paste0("chr", exprsRNA$chromosome_name),
        ranges = IRanges::IRanges(start = exprsRNA$transcription_start_site,
                                  width = 1))
    }else{
      tss <- GenomicRanges::GRanges(
        seqnames = paste0(exprsRNA$chromosome_name),
        ranges = IRanges::IRanges(start = exprsRNA$transcription_start_site,
                                  width = 1))
  }
    promoter <- suppressMessages(GenomicRanges::resize(
      tss, fix = "center", width = promoter.length))
    names(promoter) <- exprsRNA$ensembl_transcript_id_version

    files <- NULL
    for (i in unique(paste(chipMeta$mark,"_HM_",
                           chipMeta$condition, sep = ""))){
      a <- chipMeta[paste(chipMeta$mark,"_HM_", chipMeta$condition,
                          sep = "")==i,]$files
      names(a) <- paste(i,"_REP",1:length(a), sep = "")
      files <- c(a, files)
    }
    fileNames <- names(files)


    exprsChIP.list <- NULL
    for (f in fileNames){
      a <- bam2counts(bamFile = files[f],
                      region = promoter,
                      fragLength = fragLength)
      exprsChIP.list[[f]] <- a
    }
    exprsChIP <- as.data.frame(exprsChIP.list)
    ### match RNA-Seq and ChIP-Seq on the gene level

    match.rc <- merge(exprsRNA,
                      exprsChIP,
                      by.x="ensembl_transcript_id_version",
                      by.y=0)
    chip <- match.rc[,colnames(match.rc)%in%c("external_gene_name",
                                              colnames(exprsChIP))]
    rna <- match.rc[,colnames(match.rc)%in%c("external_gene_name",
                                             colnames(exprs))]
    res.rna <- stats::aggregate(. ~ external_gene_name,
                                data=rna,
                                FUN=sum)
    if(method == "weighted.mean"){
      res.chip <- stats::aggregate(
        . ~ external_gene_name,
        data= chip,
        FUN = function(x)stats::weighted.mean(x+0.01, x+0.01))
      matched.wtmean <- merge(res.rna, res.chip,
                              by = "external_gene_name")
      res.rna <- stats::aggregate(. ~ external_gene_name,
                                  data = res.rna,
                                  FUN = as.integer)
      res.chip <- stats::aggregate(. ~ external_gene_name,
                                   data = res.chip,
                                   FUN = as.integer)
      matched.wtmean <- stats::aggregate(. ~ external_gene_name,
                                         data = matched.wtmean,
                                         FUN = as.integer)
      return(list(res.rna = res.rna,
                  res.chip = res.chip,
                  matched.data = matched.wtmean))
    }else
    {
      res.chip <- stats::aggregate(. ~ external_gene_name,
                                   data = chip,
                                   FUN = max)
      matched.highest <- merge(res.rna, res.chip, by = "external_gene_name")
      res.rna <- stats::aggregate(. ~ external_gene_name,
                                  data = res.rna,
                                  FUN = as.integer)
      res.chip <- stats::aggregate(. ~ external_gene_name,
                                   data = res.chip,
                                   FUN = as.integer)
      matched.highest <- stats::aggregate(. ~ external_gene_name,
                                          data = matched.highest,
                                          FUN = as.integer)
      return(list(res.rna = res.rna,
                  res.chip = res.chip,
                  matched.data = matched.highest))
    }
  }else
  {
    rna <- exprsRNA[,colnames(exprsRNA)%in%c("external_gene_name",
                                             colnames(exprs))]
    res.rna <- stats::aggregate(. ~ external_gene_name,
                                data = rna,
                                FUN = sum)



    genebody.region <- exprsRNA[colnames(exprsRNA)%in%c(
      "external_gene_name",
      "chromosome_name",
      "start_position",
      "end_position")]
    # remove duplicated genes
    genebody.region <-
      genebody.region[!duplicated(genebody.region$external_gene_name),]

    if("chr1"%in%as.character(GenomicRanges::seqnames(aln))){
      genebody <- GenomicRanges::GRanges(
        seqnames = paste0("chr", genebody.region$chromosome_name),
        ranges = IRanges::IRanges(start = genebody.region$start_position,
                                  end = genebody.region$end_position))

    }else{
      genebody <- GenomicRanges::GRanges(
        seqnames = genebody.region$chromosome_name,
        ranges = IRanges::IRanges(start = genebody.region$start_position,
                                  end = genebody.region$end_position))
      }


    names(genebody) <-  genebody.region$external_gene_name

    files <- NULL
    for (i in unique(paste(chipMeta$mark,
                           "_HM_",
                           chipMeta$condition,
                           sep = ""))){
      a <- chipMeta[paste(
        chipMeta$mark,
        "_HM_",
        chipMeta$condition, sep = "")==i,]$files
      names(a) <- paste(i,"_REP",1:length(a), sep = "")
      files <- c(a, files)
    }
    fileNames <- names(files)


    exprsChIP.list <- NULL
    for (f in fileNames){
      a <- bam2counts(bamFile = files[f],
                      region = genebody,
                      fragLength = fragLength)
      exprsChIP.list[[f]] <- a
    }
    res.chip <- as.data.frame(exprsChIP.list)
    res.chip$external_gene_name <- rownames(res.chip)
    matched.genebody <- merge(res.rna, res.chip, by="external_gene_name")
    res.rna <- stats::aggregate(. ~ external_gene_name,
                                data = res.rna,
                                FUN = as.integer)
    res.chip <- stats::aggregate(. ~ external_gene_name,
                                 data = res.chip,
                                 FUN = as.integer)
    matched.genebody <- stats::aggregate(. ~ external_gene_name,
                                         data = matched.genebody,
                                         FUN = as.integer)
    return(list(res.rna = res.rna,
                res.chip = res.chip,
                matched.data = matched.genebody))
  }
}
