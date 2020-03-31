#' @title Compute the normalized number of reads (rpm) fall into specific
#' genomic region
#'
#' @aliases bam2rpm
#'
#' @description
#' \code{bam2rpm} computes the normalized number of reads (rpm) fall into
#' specific genomic region such as promoter, enhancer, genebody
#'
#' @param bamFile
#' Aligned bam file as input.
#'
#' @param region
#' The GRanges object defined by user to calculate the number of reads fall into
#' this specific region. For ChIP-Seq of histone modifications they are usually
#' promoter, enhancer and genebody regions.
#'
#' @param fragLength
#' Extend reads toward the 3'-end to the average DNA fragment size obtained
#' after DNA size selection
#'
#' @return a vector of numbers
#' @import Rsamtools
#' @importFrom GenomicRanges start end strand seqnames countOverlaps mcols
#' @importFrom GenomicAlignments readGAlignments
#' @import GenomeInfoDb
#' @examples
#' data("promoter")
#' file.bam <- system.file("extdata", "SRR925640.bam", package = "intePareto")
#' bam2rpm(bamFile = file.bam, region = promoter, fragLength = 180)
#' @export
#'
bam2rpm <- function(bamFile,
                    region,
                    fragLength=180){

  ### check input
  if (class(region)!="GRanges") stop("region must be a GRanges object")
  if (class(fragLength)!="numeric") stop("fragLength must be numeric")

  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE),
    what = "mapq")
  aln <- GenomicAlignments::readGAlignments(bamFile, param = param)
  aln <- aln[GenomicRanges::mcols(aln)$mapq > 0]
  aln <- GenomicRanges::granges(aln)
  # get the chromosome length
  length.chr <- GenomeInfoDb::seqlengths(seqinfo(aln))
  # to make sure the "region + fragLength" is not out-of-bound
  out.bound <- ifelse(
    GenomicRanges::strand(aln)=="+",
    GenomicRanges::start(aln)+ fragLength-1,
    GenomicRanges::end(aln))  >
    length.chr[as.character(GenomicRanges::seqnames(aln))]
  aln <- aln[GenomicRanges::end(aln) > fragLength&!out.bound]
  aln <- GenomicRanges::resize(aln, fragLength, fix = "start")
  # countOverlaps is strand aware, so remove strand
  GenomicRanges::strand(aln) <- "*"
  counts <- GenomicRanges::countOverlaps(region, aln)
  seqdepth <- length(aln)
  rpm <- round(counts/seqdepth * 10^6, digits = 3)
  names(rpm) <- names(region)
  return(rpm)
}
