#'The rqt class
#'
#'This class stores parameters and results of the rtq algorithms
#'
#'@section Slots:
#'    \describe{
#'    \item{\code{phenotype}:}{Phenotype (a vector of length 
#'        \code{N}, where \code{N} - number of individuals).}
#'      \item{\code{genotype}:}{Genotype - an object of class 
#'      \code{SummarizedExperiment}. Should contain one assay 
#'      (matrix, \code{N} 
#'      by \code{M} where \code{N} - number of individuals, \code{M}
#'       - number of genetic variants).}
#'      \item{\code{covariates}:}{data frame \code{N} 
#'      by \code{K} where \code{N} - number of individuals, \code{K}
#'       - number of covariates)}
#'      \item{\code{results}:}{A list of two: 
#'      test statistics (\code{Q1}, \code{Q2}, \code{Q3}), 
#'      p-values (\code{p1.Q1}, \code{p2.Q2}, \code{p3.Q3})}
#'}
#'@rdname rqt-class
setClass("rqt", slots=c(phenotype="vector", 
                        genotype="SummarizedExperiment", 
                        covariates="data.frame", 
                        results="list"))


