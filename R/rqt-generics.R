
#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param obj Object of class \code{rqt}
#' @param ... Additional parameters to pass to the function
#' @return Updated rqt object with result slot
#' @docType methods
#' @rdname rqt-geneTest
setGeneric("geneTest", function(obj, ...) standardGeneric("geneTest"))

#' This function performs a gene-level meta-analysis based on 
#' combined effect sizes.
#' 
#' @param objects List of objects of class rqt
#' @param ... Additional parameters to pass to the function
#' @return A list of two: (i) final.pvalue - 
#' a final p-value across all studies; 
#' (ii) pvalueList - p-values for each study; 
#' @docType methods
#' @rdname rqt-geneTestMeta
setGeneric("geneTestMeta", function(objects, ...) standardGeneric("geneTestMeta"))

#' This function performs an access to phenotype
#' 
#' @rdname rqt-phenotype
setGeneric("phenotype", function(obj) standardGeneric("phenotype"))

#' This function performs an access to genotype.
#' 
#' @rdname rqt-genotype
setGeneric("genotype", function(obj) standardGeneric("genotype"))

#' This function performs an access to covariates
#' 
#' @rdname rqt-covariates
setGeneric("covariates", function(obj) standardGeneric("covariates"))

#' This function performs an access to covariates
#' 
#' @rdname rqt-results
setGeneric("results", function(obj) standardGeneric("results"))
setGeneric("results<-", function(obj, value) standardGeneric("results<-"))