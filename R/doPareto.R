#' @title Prioritization of genes based on Z scores
#'
#' @aliases doPareto
#'
#' @description
#' \code{doPareto} takes the Z scores of several different histone modifications
#' as input, the prioritization of genes based on Z scores can be formulated as
#' multiobjective optimization problem and solved with Pareto optimization.
#'
#' @param df_final
#' a data frame which is the output of doIntegration.
#'
#' @param objective
#' a data frame which include column of "mark" column indicates the z scores of
#' markers of histone modifications (e.g. "z.H3K4me3"), and a column named "obj"
#' indicates the direction of the operation on the z scores, one of "max" and
#' "min".
#'
#' @param nr.fronts
#' the number of the pareto fronts you want to get.
#'
#' @return a data.frame ranked by the level of pareto fronts.
#'
#' @import rPref
#' @examples
#' data("df_final")
#' objective <- data.frame(mark = c("z.H3K27ac","z.H3K4me3"),
#'                         obj=c("max","max"),stringsAsFactors=FALSE)
#' nr.fronts <- 3
#' doPareto(df_final = df_final,
#'          objective = objective,
#'          nr.fronts = nr.fronts)
#' @export
#'

doPareto <- function(df_final,
                     objective,
                     nr.fronts){
  ### check input
  checkInputParam <- function(df_final, objective, nr.fronts) {

    # simple input param checks
    if(missing(df_final) || is.null(df_final) ||
       missing(objective) || is.null(objective) ||
       missing(nr.fronts) || is.null(nr.fronts)) {
      stop("arguments df_final, objective, nr.fronts must be specified")
    }

    # check rnaMeat
    if(class(df_final)!="data.frame"){
      stop("df_final must be a data.frame")
    }
    if(class(objective)!="data.frame"){
      stop("objective must be a data.frame")
    }
    if(!c("mark")%in%colnames(objective)){
      stop('objective has to contain "mark" column to indicate z scores you
           want to operate')
    }
    if(!c("obj")%in%colnames(objective)){
      stop('objective has to contain "obj" column to indicate the objective
           of max or min')
    }
    if(!all(objective$obj %in% c("min","max"))){
      stop('the value of "obj" column of objective has to be "max" or "min"')
    }
    if(length(unique(objective$obj)) > 2) {
      stop('the value of "obj" column of objective has to be "max" or "min"')
    }
    if(!all(objective$mark %in% colnames(df_final))){
      stop('the value of "mark" column of objective has to be one of the column
           of the df_final')
    }
  }

  checkInputParam(df_final = df_final,
                  objective = objective,
                  nr.fronts = nr.fronts)

  p <- rPref::empty()
  for (i in 1:nrow(objective)) {
    if(objective[i,]$obj=="max"){
      m <- as.character(objective[i,]$mark)
      a <- rPref::high_(m)
    }else
    {
      m <- as.character(objective[i,]$mark)
      a <- rPref::low_(m)
    }
    p <- p*a
  }

  res <- rPref::psel(df = df_final, pref = p, top_level = nr.fronts)
  names(res)[names(res) == '.level'] <- 'front'
  return(res)
}
