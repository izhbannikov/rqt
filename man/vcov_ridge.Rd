% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{vcov_ridge}
\alias{vcov_ridge}
\title{vcov_ridge: returns variance-covariance 
matrix and standard deviation
for ridge/LASSO regression object}
\usage{
vcov_ridge(x, y, rmod, verbose = FALSE)
}
\arguments{
\item{x}{Genotype matrix}

\item{y}{Phenotype}

\item{rmod}{Ridge/LASSO regression object}

\item{verbose}{Indicates verbosing output,
Default: FALSE.}
}
\value{
list(vcov, se). 
vcov: variance-covariance matrix;
se: standard deviation
}
\description{
vcov_ridge: returns variance-covariance 
matrix and standard deviation
for ridge/LASSO regression object
}
