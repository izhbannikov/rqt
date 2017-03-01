#' The rqt class constructor
#' 
#' This function generates rqt class objects
#' @param  phenotype Phenotype (a vector of length 
#'        \code{N}, where \code{N} - number of individuals).
#' @param genotype Genotype - an object of class 
#'      \code{SummarizedExperiment}. Should contain one assay 
#'      (matrix, \code{N} 
#'      by \code{M} where \code{N} - number of individuals, \code{M}
#'       - number of genetic variants).
#' @param covariates Covariates, a data frame \code{N} 
#'      by \code{K} where \code{N} - number of individuals, \code{K}
#'       - number of covariates
#' @param results A list of two: test statistics: 
#' (\code{Q1}, \code{Q2}, \code{Q3}), 
#' p-values: (\code{p1.Q1}, \code{p2.Q2}, \code{p3.Q3})
#' @return Object of class \code{rqt}
#' @examples
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' print(obj)
#' @rdname rqt-methods
#' @export
rqt <- function(phenotype=NULL, genotype=NULL, covariates=NULL, 
                results=NULL) {
  
  if(is.null(phenotype)) {
    phenotype <- c()
  }
  
  if(is.null(genotype)) {
    genotype <- matrix()
    colnames(genotype) <- "geno"
    genotype <- SummarizedExperiment(genotype)
  } 
  
  if(is.null(covariates)) {
    covariates <- data.frame()
  }
  
  if(is.null(results)) {
    results <- list()
  }
  
  new("rqt", phenotype=phenotype, 
      genotype=genotype, 
      covariates=covariates, 
      results=results)
}


#' geneTest
#' This function performs a gene-level test based on combined effect sizes.
#' 
#' @param perm Integer indicating the number of permutations 
#' to compute p-values. Default: 0.
#' @param STT Numeric indicating soft truncation threshold (STT) 
#' to convert to gamma parameter (must be <= 0.4). 
#' Needed for an optimal parameter a in Gamma-distribution. Default: 0.2. 
#' See, for example, Fridley, et al 2013: "Soft truncation thresholding 
#' for gene set analysis of RNA-seq data: Application to a vaccine study".
#' @param weight Logical value. Indicates using weights (see Lee et al 2016). 
#' Default: FALSE.
#' @param cumvar.threshold Numeric value indicating 
#' the explained variance threshold for PCA-like methods. Default: 75
#' @param method Method used to reduce multicollinerity and account for LD. 
#' Default: PLS-DA.
#' @param out.type Character, indicating a type of phenotype. 
#' Possible values: \code{D} (dichotomous or binary), 
#' \code{C} (continous or quantitative).
#' @param scaleData A logic parameter (TRUE/FALSE) indicating scaling of 
#' the genotype dataset.
#' @param asym.pval Indicates Monte Carlo approximation for p-values. Default: FALSE.
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @rdname rqt-geneTest
#' @export
#' @examples
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat",
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' res <- geneTest(obj, method="pca", out.type = "D")
#' print(res)
setMethod("geneTest", signature = "rqt", 
    function(obj, perm=0, STT=0.2, weight=FALSE, 
            cumvar.threshold=75, out.type="D", 
            method="pca", scaleData=FALSE, asym.pval=FALSE,
            verbose=FALSE) {
            # Prepare test: load distribution table and prepare #
            # some other information #
        if(cumvar.threshold > 100) {
            warning("Warning: cumvar.threshold > 100 and will be set to 100.")
            cumvar.threshold <- 100
        }
      
        # Load data #
        phenotype <- phenotype(obj)
        genotype <- assays(genotype(obj))[[1]]
        covariates <- covariates(obj)
        
        # Dimensions
        phenoSize <- length(phenotype)
        genoSize <- dim(genotype)
        
        if(weight & scaleData) {
          if(verbose) {
              print("Warining! You can not use scaling in presence of weights!")
              print("Parameter weight will be set to FALSE.")
          }
          weight <- FALSE
        }
        
        # Start the tests #
        if(perm==0){
            rslt0 <- geneTestOne(phenotype=phenotype, 
                    genotype=genotype, 
                    covariates=covariates, 
                    STT=STT, 
                    weight=weight,
                    cumvar.threshold=cumvar.threshold, 
                    out.type=out.type, method=method, 
                    scaleData = scaleData, verbose=verbose)
            
            if(asym.pval) {
                rsltMC <- do.call(rbind, 
                                  lapply(1:dim(genotype)[1], 
                                         function(k){
                    yP <- phenotype[sample(1:phenoSize, phenoSize, 
                                           replace=FALSE)]
                    t.res <- geneTestOne(phenotype=yP,genotype=genotype, 
                                 covariates=covariates,STT=STT,
                                 weight=weight,
                                 cumvar.threshold=cumvar.threshold, 
                                 out.type=out.type, method=method, 
                                 scaleData = scaleData)
                    if(is.na(t.res)) {
                        t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                               data.frame(pVal1=1,pVal2=1,pVal3=1) )
                        names(t.res) <- c("Qstatistic", "pValue")
                    }
                    tt.res <- t.res$Qstatistic
                }))
            
                if(!is.na(rslt0)) {
                    nn <- genoSize[1]+1
                    rslt <- list("Qstatistic"= data.frame(
                            Q1=rslt0$Qstatistic$Q1,
                            Q2=rslt0$Qstatistic$Q2, 
                            Q3=rslt0$Qstatistic$Q3),
                           "pValue" = data.frame(
                             pVal1 = (length(rsltMC[,1][rsltMC[,1] >= 
                                            rslt0$Qstatistic$Q1])+1)/nn,
                             pVal2 = (length(rsltMC[,2][rsltMC[,2] >= 
                                            rslt0$Qstatistic$Q2])+1)/nn,
                             pVal3 = (length(rsltMC[,3][rsltMC[,3] >= 
                                            rslt0$Qstatistic$Q3])+1)/nn),
                           beta = rslt0$beta)
                } else {
                    rslt <- list( Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA),
                                pValue=data.frame(pVal1=1,pVal2=1,pVal3=1) )
                }
            } else {
                rslt <- rslt0
            }
        } else {
            rsltPP <- do.call(rbind, lapply(1:perm, function(k){
                yP <- phenotype[sample(1:phenoSize, phenoSize,replace=FALSE)]
                t.res <- geneTestOne(phenotype=yP,genotype=genotype, 
                    covariates=covariates,STT=STT,
                    weight=weight,
                    cumvar.threshold=cumvar.threshold, 
                    out.type=out.type, method=method, 
                    scaleData = scaleData, verbose=verbose)
                if(is.na(t.res)) {
                    t.res <- list( data.frame(Q1=NA, Q2=NA, Q3=NA), 
                        data.frame(pVal1=1,pVal2=1,pVal3=1) )
                    names(t.res) <- c("Qstatistic", "pValue")
                }
                pv <- t.res$pValue
            }))
              
            rslt0 <- as.numeric(geneTestOne(phenotype=phenotype, 
                genotype=genotype, 
                covariates=covariates, 
                STT=STT, weight=weight, 
                cumvar.threshold=cumvar.threshold, 
                out.type=out.type, method=method, 
                scaleData = scaleData))
              
            if(!is.na(rslt0)) {
                rslt <- list(Qstatistic= data.frame(Q1=NA, Q2=NA, Q3=NA),
                    pValue = data.frame(
                    pVal1 = (length(rsltPP[,1][rsltPP[,1] < 
                                              rslt0$pValue[1]])+1)/(perm+1),
                    pVal2 = (length(rsltPP[,2][rsltPP[,2] < 
                                              rslt0$pValue[2]])+1)/(perm+1),
                    pVal3 = (length(rsltPP[,3][rsltPP[,3] < 
                                              rslt0$pValue[3]])+1)/(perm+1)),
                    beta = rslt0$beta)
            } else {
                rslt <- list( Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA),
                            pValue=data.frame(pVal1=1,pVal2=1,pVal3=1) )
            }
        }
        
        results(obj) <- rslt
        return(obj)
})


#' This function performs a gene-level meta-analysis based on 
#' combined effect sizes.
#' 
#' @param perm Integer indicating the number of permutations 
#' to compute p-values. Default: 0.
#' @param STT Numeric indicating soft truncation threshold (STT) 
#' to convert to gamma parameter (must be <= 0.4). 
#' Needed for an optimal parameter a in Gamma-distribution. Default: 0.2. 
#' See, for example, Fridley, et al 2013: "Soft truncation thresholding 
#' for gene set analysis of RNA-seq data: Application to a vaccine study".
#' @param weight Logical value. Indicates using weights (see Lee et al 2016). 
#' Default: FALSE.
#' @param cumvar.threshold Numeric value indicating 
#' the explained variance threshold for PCA-like methods. Default: 75
#' @param method Method used to reduce multicollinerity and account for LD. 
#' Default: PLS-DA.
#' @param out.type Character, indicating a type of phenotype. 
#' Possible values: \code{D} (dichotomous or binary), 
#' \code{C} (continous or quantitative).
#' @param scaleData A logic parameter (TRUE/FALSE) indicating scaling of 
#' the genotype dataset.
#' @param asym.pval Indicates Monte Carlo approximation for p-values. Default: FALSE.
#' @param comb.test Statistical test for combining p-values.
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of two: (i) final.pvalue - 
#' a final p-value across all studies; 
#' (ii) pvalueList - p-values for each study; 
#' @export
#' @docType methods
#' @rdname rqt-geneTestMeta
#' @examples
#'data1 <- data.matrix(read.table(system.file("extdata/phengen2.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data1[,1]
#'geno <- data1[, 2:dim(data1)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj1 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'data2 <- data.matrix(read.table(system.file("extdata/phengen3.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data2[,1]
#'geno <- data2[, 2:dim(data2)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj2 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'data3 <- data.matrix(read.table(system.file("extdata/phengen.dat",
#'                                            package="rqt"), skip=1))
#'pheno <- data3[,1]
#'geno <- data3[, 2:dim(data3)[2]]
#'colnames(geno) <- paste(seq(1, dim(geno)[2]))
#'geno.obj <- SummarizedExperiment(geno)
#'obj3 <- rqt(phenotype=pheno, genotype=geno.obj)
#'
#'res.meta <- geneTestMeta(list(obj1, obj2, obj3))
#'print(res.meta)
setMethod("geneTestMeta", signature="list", 
          function(objects, perm=0, STT=0.2, weight=FALSE, 
                   cumvar.threshold=75, out.type="D", 
                   method="pca", scaleData=FALSE, asym.pval=FALSE,
                   comb.test="wilkinson",
                   verbose=FALSE) {
            
            if(cumvar.threshold > 100) {
              warning("Warning: cumvar.threshold > 100 
                      and will be set to 100.")
              cumvar.threshold <- 100
            }
            if(class(objects) != "list") {
              stop("objects must be a list of rqt class objects!")
            }
            ### Meta-analysis ###
            numStudies <- length(objects)
            pv <- rep(NA_real_, numStudies) 
            for(i in 1:numStudies) {
              res <- geneTest(objects[[i]],
                              STT=STT, 
                              weight=weight, 
                              cumvar.threshold=cumvar.threshold, 
                              out.type=out.type, 
                              method=method, 
                              perm = perm, 
                              scaleData=scaleData,
                              asym.pval=asym.pval,
                              verbose=verbose)
              
              if(length(results(res)) != 0) {
                pv[i] <- results(res)$pValue$pVal3
              }
              
            }
            
            ### Combining p-values via some comb.test ###
            comb.res <- list()
            switch(comb.test, 
                   wilkinson={
                     # Wilkinson
                     comb.res <- wilkinsonp(pv)
                   },
                   fisher={
                     # Fisher
                     chi.comb <- sum(-2*log(pv[!is.na(pv)]))
                     df <- 2*length(pv)
                     comb.res[["p"]] <- 1-pchisq(q=chi.comb, df=df)
                   },
                   minimump={
                     # minimump
                     comb.res <- minimump(pv)
                   },
                   sump={
                     # sump
                     comb.res <- sump(pv)
                   },
                   sumlog={
                     # sumlog
                     comb.res <- sumlog(pv)
                   },
                   meanp={
                     comb.res <- meanp(pv)
                   },
                   logitp={
                     comb.res <- logitp(pv)
                   },
                   votep={
                     comb.res <- votep(pv)
                   },
                   {
                     # Wilkinson
                     comb.res <- wilkinsonp(pv)
                   }
            )
            
            #### End of combining p-values ####
            ### End of meta-analysis ###
            
            return(list(final.pvalue=comb.res[["p"]], 
                        pvalueList=pv))
          })


#' Get a given STT
#' @param L TODO
#' @param STT Numeric indicating soft truncation threshold (STT) to convert 
#' to gamma parameter (must be <= 0.4).
#' @return a TODO
get.a <- function(L,STT=0.2) {
  aa <- diff <- seq(0,1,length=200)
  asize <- length(aa)
  diff <- unlist(lapply(1:asize, 
                        function(i) { 
                          abs(get.stt(L,aa[i], STT) - STT) 
                        } ))
  return(aa[which.min(diff)])
}

get.stt <- function(L,a,STT=0.2){
  ans <- 1-pgamma(L*qgamma(1-STT,a,1),L*a,1)
  return(ans)
}


get.reg.family <- function(out.type="D") {
  if(out.type == "D") {
    reg.family <- "binomial"
  } else if(out.type == "C") {
    reg.family <- "gaussian"
  } else {
    stop(paste("Unknown out.type:", out.type))
  }
  return(reg.family)
}

# geneTestOne: the basic function for 
# gene-level association test
# Some parts of this code were adopted from the 
# supplementary code provided by Lee et al 2016 
# "Gene-set association tests for next-generation sequencing data"
#
geneTestOne <- function(phenotype, genotype, covariates, STT=0.2, weight=FALSE,
                        cumvar.threshold=75, method="pca", out.type="D", 
                        scaleData=FALSE, verbose=FALSE) {
  
  reg.family <- get.reg.family(out.type)
  
  rslt <- list(Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA), 
               pValue=data.frame(pVal1=NA, pVal2=NA, pVal3=NA), 
               beta=NA)
  
  res <- list()
  tryCatch({
    
    if(dim(genotype)[2] > 1) {
      
      # Removing constant columns #
      preddata <- data.frame(genotype[,apply(genotype, 2, 
                                             var, na.rm=TRUE) != 0])
      
      indexes <- NA
      ### Dimensionality reduction and account for LD ###
      if(method != "none") {
        tmp <- cor(preddata)
        tmp[upper.tri(tmp)] <- 0
        diag(tmp) <- 0
        preddata <- preddata[,!apply(tmp,2,function(x) any(x > 0.99))]
        res <- preprocess(data=preddata, 
                          pheno=phenotype, method=method,
                          reg.family=reg.family, 
                          scaleData=scaleData, 
                          cumvar.threshold=cumvar.threshold, 
                          out.type=out.type,
                          verbose=verbose)
      } else {
        res[["S"]] <- preddata
      }
      
      #### Regression after data preprocessing ####
      if(!(method %in% c("lasso", "ridge"))) {
        S <- res[["S"]]
        if(method == "pca") {
          indexes <- res[["indexes"]]
          #    print(dim(S))
        }
        # Build null model #
        null.model <- build.null.model(y=phenotype, x=covariates,
                                       reg.family=reg.family)
        
        res <- simple.multvar.reg(null.model=null.model, Z=S)
        
        S <- res[["S"]]
        fit <- res[["fit"]]
        if(dim(coef(summary(fit)))[2] >= 2) {
          coef.multivar <- coef(summary(fit))[-1,1:2]
        }
        
        if(length(coef.multivar) != 2) { 
          beta.multivar <- coef.multivar[,1]
          se.multivar <- coef.multivar[,2] 
        }
        if(length(coef.multivar) == 2) { 
          beta.multivar <- coef.multivar[1]
          se.multivar <- coef.multivar[2] 
        }
        
        vMat <- vcov(fit)[-1,-1]
        alpha <- 1/(se.multivar^2)
        alpha <- alpha/sum(1/(se.multivar^2))
        
        if(length(coef(fit)) > 2) {
          mean.vif <- mean(vif(fit), na.rm=TRUE)
        } else {
          mean.vif <- summary(fit)$coefficients[2,2]
        }
        
      } else {
        
        fit <- res[["fit"]]
        coef.multivar <- coef(fit)[-1]
        S <- preddata
        
        if(sum(coef.multivar) == 0) {
          print("All coefficients in lasso/ridge regression 
                are equal to 0. 
                Trying ordinary regressing instead.")
          
          
          # Build null model #
          null.model <- build.null.model(y=phenotype, x=covariates,
                                         reg.family=reg.family)
          
          res <- simple.multvar.reg(null.model=null.model, Z=S)
          S <- res$S
          fit <- res[["fit"]]
          if(dim(coef(summary(fit)))[2] >= 2) {
            coef.multivar <- coef(summary(fit))[-1,1:2]
          }
          
          
          if(length(coef.multivar)!=2){
            beta.multivar<-coef.multivar[,1]
            se.multivar <- coef.multivar[,2]
          }
          if(length(coef.multivar)==2){
            beta.multivar<-coef.multivar[1]
            se.multivar <- coef.multivar[2]
          }
          vMat <- vcov(fit)[-1,-1]
        } else {
          
          beta.multivar <- coef.multivar[coef.multivar != 0L]
          vMat <- vcov_ridge(x=as.matrix(preddata), 
                             y=phenotype, 
                             rmod=fit)$vcov[coef.multivar != 0L, coef.multivar != 0L]
          
          if(class(vMat)[1] == "dgeMatrix") {
            se.multivar <- sqrt(diag(vMat))
          } else {
            se.multivar <- sqrt(vMat)
          }
          vMat <- as.matrix(vMat)
        }
        
        alpha <- as.matrix((1/(se.multivar^2)), ncol=1)
        alpha <- alpha/sum(1/(se.multivar^2))
        
        # Temporary disabled:
        #mean.vif <- mean(vif(fit), na.rm=TRUE)
        #print(mean.vif)
        mean.vif <- NA
      }
      
      beta.pool.base <- 0
      beta.pool <- 0
      var.pool <- NA
      var.pool.base <- NA
      
      if(length(coef.multivar) != 0) {
        #=================== QTest1 ===================#
        if(weight == FALSE) {
          var.pool.base <- t(alpha) %*% vMat %*% alpha
          beta.pool.base <- t(alpha) %*% beta.multivar
          zScore.base <- beta.pool.base/sqrt(var.pool.base)
          QStat1 <- zScore.base^2
          pVal1 <- pchisq(QStat1, df=1, lower.tail=FALSE)
        } else {
          if(!(method %in% c("pca", "pls"))) {
            maf.S <- apply(S,2, function(v)mean(na.omit(v))/2)
          } else {
            maf.S <- apply(S,2, 
                           function(v) mean(na.omit(abs(v)))/2 )
          }
          
          w.S0 <- qbeta(maf.S,1,25,lower.tail=FALSE)
          
          WS <- diag(w.S0)
          if(length(beta.multivar)==1){
            WS <- w.S0
          }
          
          var.pool <- t(alpha) %*% WS %*% vMat %*% WS %*% alpha
          beta.pool <- t(alpha) %*% WS %*% beta.multivar
          z.score <- beta.pool/sqrt(var.pool)
          QStat1 <- z.score^2
          pVal1 <- pchisq(QStat1, df=1, lower.tail=FALSE)
        }
        
        #=================== QTest-2 ===================#
        QStat2.eigen <- eigen(vMat)
        U2 <- QStat2.eigen$vectors
        l2 <- QStat2.eigen$values
        na.l2 <- which(l2/mean(l2) < 0.01)
        if(length(na.l2) > 0){
          l2 <- l2[-na.l2]
          U2 <- U2[,-na.l2]
        }
        
        p2 <- pchisq((t(U2) %*% beta.multivar)^2/l2, df=1, lower.tail=FALSE)
        a <- get.a(length(p2), STT)
        q2 <- 2*(qgamma(p2, a, 1, lower.tail=FALSE))
        QStat2 <- sum(q2)
        pVal2 <- pchisq(QStat2, df=2*a*length(q2), lower.tail=FALSE)
        
        #=================== QTest-3 ===================#
        
        if(weight==FALSE){
            cov.beta <- c(t(alpha) %*% vMat)
            b.star <- beta.multivar - beta.pool.base*cov.beta/var.pool.base[1]
            vMat.star <- vMat - (cov.beta%*%t(cov.beta))/var.pool.base[1]
        } else {
            w.vMat<-WS %*% vMat %*% WS
            cov.beta <- c(t(alpha)%*%w.vMat)
            b.star <- WS %*% beta.multivar - beta.pool*cov.beta/var.pool[1]
            vMat.star <- w.vMat - (cov.beta%*%t(cov.beta))/var.pool[1]
        }
        
        if(length(beta.multivar) !=1 ) {
            QStat3.eigen <- eigen(vMat.star)
            U3 <- QStat3.eigen$vectors
            l3 <- QStat3.eigen$values
            na.l3 <- which(l3/mean(l3) < 0.001)
          
            if(length(na.l3) > 0) {
                l3<-l3[-na.l3]
                U3<-U3[,-na.l3]
            }
          
          
          q2.proj <- (t(U3) %*% b.star)^2/l3
          p2.1 <- pchisq(q2.proj, df=1, lower.tail=FALSE)
          a <- get.a(length(p2.1),STT)
          QStat2.proj <- sum(2*qgamma(p2.1,a,1,lower.tail=FALSE))
          pVal2.proj <- pchisq(QStat2.proj,df=2*a*length(l3),
                              lower.tail=FALSE)
          if(pVal2.proj == 0) {
            pVal2.proj <- 1e-8
          }
          QStat2.1 <- qchisq(pVal2.proj,df=1,lower.tail=FALSE)
          
          pi0 <- seq(0,1,by=0.1)
          pVal3.can <- QStat3 <- rep(1,11)
          pVal3.can[1] <- pVal2.proj
          QStat3[1] <- QStat2.1
          pVal3.can[11] <- pVal1
          QStat3[11] <- QStat1
          
          for(h in 2:10){
            QStat3[h] <- pi0[h]*QStat1 + (1-pi0[h])*QStat2.1
            pVal3.can[h] <- davies(QStat3[h],
                                  c(pi0[h],(1-pi0[h])),
                                  c(1,1))$Qq
            if(pVal3.can[h] <= 0 | pVal3.can[h] > 1) {
              pVal3.can[h] <- imhof(QStat3[h],
                                   c(pi0[h],(1-pi0[h])),
                                   c(1,1))$Qq
            }
            if(pVal3.can[h]<=0|pVal3.can[h]>1){
              pVal3.can[h]<-liu(QStat3[h],c(pi0[h],
                                           (1-pi0[h])),c(1,1))[1]
            }
          }
          
          QStat3final <- QStat3[which.min(pVal3.can)]
          pVal3 <- (sum(null.dist.Q3[null.dist.Q3[,1] > 
                                      -log10(min(pVal3.can)),2])+1) / 
            (sum(null.dist.Q3[,2])+1)
        }
        
        if(length(beta.multivar)==1){
          pVal3 <- pVal1
          QStat3final <- QStat1
        }
        
        rslt <- list(Qstatistic=data.frame(Q1=QStat1, Q2=QStat2, Q3=QStat3final), 
                     pValue=data.frame(pVal1,pVal2,pVal3),
                     beta=ifelse(weight, beta.pool, beta.pool.base),
                     var.pooled=ifelse(weight, var.pool, var.pool.base),
                     mean.vif=mean.vif)
      }
      
      if(length(coef.multivar)==0) { 
        rslt <- list(Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA), 
                     pValue=data.frame(pVal1=1,pVal2=1,pVal3=1),
                     beta=NA, var.pooled=NA, mean.vif=NA)
        #rslt <- NA
      }
    } else {
      # Simple regression:
      # Build null model #
      
      null.model <- build.null.model(y=phenotype, x=covariates,
                                     reg.family=reg.family)
      
      res <- simple.multvar.reg(null.model=null.model, Z=genotype)
      reg.coef <- coef(summary(res[["fit"]]))
      
      if(length(coef(res[["fit"]])) > 2) {
        mean.vif <- mean(vif(res[["fit"]]), na.rm=TRUE)
      } else {
        mean.vif <- summary(res[["fit"]])$coefficients[2,2]
      }
      
      if(dim(reg.coef)[1] == 2) {
        rslt <- list(Qstatistic=data.frame(Q1=reg.coef[2,3], 
                                           Q2=reg.coef[2,3], Q3=reg.coef[2,3]), 
                     pValue=data.frame(pVal1=reg.coef[2,4],
                                        pVal2=reg.coef[2,4],pVal3=reg.coef[2,4]),
                     beta=reg.coef[2,1],
                     var.pooled=reg.coef[2,2],
                     mean.vif=mean.vif)
      } else {
        rslt <- list(Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA), 
                     pValue=data.frame(pVal1=NA,pVal2=NA,pVal3=NA), 
                     beta=NA, var.pooled=NA, mean.vif=NA)
        #rslt <- NA
      }
      
      
    }
  },error=function(e) {
    print(e)
    rslt <- list(Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA), 
                 pValue=data.frame(pVal1=NA,pVal2=NA,pVal3=NA), 
                 beta=NA, var.pooled=NA, mean.vif=NA)
    #rslt <- NA
  }, finally=rslt)
  
  if(is.na(rslt$pValue$pVal3)) {
    rslt <- list(Qstatistic=data.frame(Q1=NA, Q2=NA, Q3=NA), 
                 pValue= data.frame(pVal1=NA,pVal2=NA,pVal3=NA), 
                 beta=NA, var.pooled=NA, mean.vif=NA)
    #rslt <- NA
  }
  
  return(rslt)
}


#' Common methods for class rqt
#' 
#' @name rqt-general
#' @rdname rqt-general
#'
#' @aliases show.rqt
#' @aliases summary.rqt
#' @aliases print.rqt
#' 
#' @title General functions of \code{rqt} 
#' such as accessors and printing.
#' 
#' @description Common methods for class rqt. 
#' This document lists a series of basic methods for the class rqt
#' 
NULL

#' This function performs an access to phenotype
#' 
#' @description A phenotype accessor
#' @param obj An object of \code{rqt} class.
#' @return phenotype returns the phenotype
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' phenotype(obj)
#' @rdname rqt-phenotype
#' @docType methods
#' @export
setMethod("phenotype", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "phenotype"))
})


#' This function performs an access to genotype.
#' 
#' @description A genotype accessor
#' @param obj An object of \code{rqt} class.
#' @return genotype returns the genotype
#' @rdname rqt-genotype
#' @docType methods
#' @export
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' genotype(obj)
setMethod("genotype", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "genotype"))
})


#' This function performs an access to covariates
#' 
#' @description An accessor to covariates
#' @param obj An object of \code{rqt} class.
#' @return covariates returns the covariates
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' covariates(obj)
#' @rdname rqt-covariates
#' @docType methods
#' @export
setMethod("covariates", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "covariates"))
})

#' This function performs an access to covariates
#' 
#' @description An accessor to results
#' @param obj An object of \code{rqt} class.
#' @return results returns the results
#' @examples 
#' data <- data.matrix(read.table(system.file("extdata/test.bin1.dat", 
#' package="rqt"), header=TRUE))
#' pheno <- data[,1]
#' geno <- data[, 2:dim(data)[2]]
#' colnames(geno) <- paste(seq(1, dim(geno)[2]))
#' geno.obj <- SummarizedExperiment(geno)
#' obj <- rqt(phenotype=pheno, genotype=geno.obj)
#' res <- geneTest(obj, method="pca", out.type = "D")
#' results(res)
#' @rdname rqt-results
#' @docType methods
#' @export
setMethod("results", signature(obj = "rqt"), function(obj) {
  return(slot(obj, "results"))
})

setReplaceMethod("results", signature(obj = "rqt"), 
                 function(obj, value) {
                     slot(obj, "results") <- value
                     obj
                 }
)

setMethod("show", signature(object = "rqt"), function(object) {
  cat("Phenotype:\n")
  print(head(phenotype(object)))
  cat("...\n\n")
  cat("Genotype:\n")
  print(head(assays(genotype(object))[[1]]))
  cat("...\n\n")
  cat("Covariates:\n")
  print(head(covariates(object)))
  cat("\n\n")
  cat("Results:\n\n")
  print(results(object))
})

setMethod("summary", signature="rqt", function(object) {
  cat("Phenotype:\n")
  print(summary(phenotype(object)))
  cat("...\n\n")
  cat("Genotype:\n")
  print(summary(assays(genotype(object))[[1]]))
  cat("...\n\n")
  cat("Covariates:\n")
  cat(summary(covariates(object)))
  cat("\n\n")
  cat("Results:\n\n")
  print(summary(results(object)))
})
