#' Importing required packages and functions
#' 
#' @importFrom methods new slot slot<-
#' @importFrom stats binomial cor glm
#' @importFrom stats pchisq pgamma prcomp 
#' @importFrom stats qbeta qchisq qgamma
#' @importFrom stats resid vcov var
#' @importFrom stats na.exclude na.omit
#' @importFrom stats prcomp
#' @importFrom stats coef predict
#' @importFrom utils head
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom ropls opls
#' @importFrom Matrix Matrix
#' @importFrom CompQuadForm davies imhof liu
#' @importFrom metap wilkinsonp minimump sump
#' @importFrom metap sumlog meanp logitp
#' @importFrom metap votep wilkinsonp
#' @importFrom car vif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay assays
#' @importFrom RUnit checkEqualsNumeric
#' @importFrom pls plsr cppls explvar selectNcomp
#'
NULL


ridge_se <- function(xs,y,yhat,my_mod, verbose=FALSE){
  # Note, you can't estimate an intercept here
  n <- dim(xs)[1]
  k <- dim(xs)[2]
  sigma_sq <- sum((y - yhat)^2)/ (n-k)
  lam <- my_mod$lambda.min
  if(is.null(my_mod$lambda.min)) {
    lam <- 0
  }
  i_lams <- Matrix(diag(x=1,nrow=k,ncol=k),sparse=TRUE)
  xpx <- t(xs) %*% xs
  xpxinvplam <- solve(xpx + lam*i_lams)
  var_cov <- sigma_sq * (xpxinvplam %*% xpx %*% xpxinvplam)
  se_bs <- sqrt(diag(var_cov))
  
  if(verbose) 
    print('NOTE: These standard errors are very biased.')
  
  return(list(vcov=var_cov, se=se_bs))
}

#' vcov_ridge: returns variance-covariance 
#' matrix and standard deviation
#' for ridge/LASSO regression object
#' @param x Genotype matrix
#' @param y Phenotype
#' @param rmod Ridge/LASSO regression object
#' @param verbose Indicates verbosing output,
#' Default: FALSE.
#' @return list(vcov, se). 
#' vcov: variance-covariance matrix;
#' se: standard deviation
vcov_ridge <- function(x, y,  rmod, verbose=FALSE) {
    # Predictions
    r_yhat   <- predict(rmod,newx=x,s='lambda.min')
    ro_yhat  <- predict(rmod,newx=x)
    # Variance-covariance matrix and Standard Erros
    rmod_ses <- ridge_se(x,y,r_yhat,rmod,verbose=verbose)
    return(rmod_ses)
}

#' Preprocess input data with Principal Component Analysis method (PCA)
#' @param data An input matrix with values of 
#' independent variables (predictors).
#' @param pheno A phenotype - column-vector, needed for LASSO/ridge and 
#' \code{NULL} by default.
#' @param method A dimensionality reduction method.
#' Default: \code{pca}.
#' @param reg.family A regression family. 
#' Default: \code{"binomial"}.
#' @param scaleData A logical variable, indicates wheither or 
#' not scaling should be performed. Default: \code{FALSE}.
#' @param cumvar.threshold A threshold value for explained variance.
#' Default: \code{75}
#' @param out.type An output (phenotype) type. 
#' Default: \code{"D"}
#' @param penalty Value of penalty parameter for LASSO/ridge regression. 
#' Default: \code{0.001}
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of one: "S" - a data frame of predictor values.
preprocess <- function(data, pheno=NULL,
                       method="pca",
                       reg.family="binomial", 
                       scaleData=FALSE, 
                       cumvar.threshold=75,
                       out.type="D",
                       penalty=0.001,
                       verbose=FALSE) {
    
    
    switch(method, 
        pca={
            return(preprocessPCA(data, scaleData, cumvar.threshold, verbose))
        },
        pls={
            return(preprocessPLS(data, pheno, scaleData, cumvar.threshold, out.type))
        },
        lasso={
            return(preprocessLASSO(data, pheno, reg.family, penalty))
        },
        ridge={
            return(preprocessRidge(data, pheno, reg.family, penalty))
        },
        {
            stop("Unknown method provided.")
        }
  )
}


#' Applies linear of logistic regregression to the data.
#' @param null.model A fitted null model
#' @param Z A genotype matrix
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of two: "S" - a dataframe with predictors and "fit"
#' - an object returned by "glm" function.
simple.multvar.reg <- function(null.model, Z, verbose=FALSE) {
    # Fit regression according to provided null model #
    fit <- glm(null.model ~ ., data=data.frame(Z))
    na.S <- which(is.na(coef(fit)[-1]) == TRUE)

    if( length(na.S) > 0 & (dim(data.frame(Z))[2] > 1) ){
        S <- as.matrix(Z[,-na.S])
        fit <- glm(null.model ~ . ,data=data.frame(S))
    } else {
        S <- Z
    }
    return(list(S=S, fit=fit, na.S=na.S))
}

#' Applies linear of logistic regregression to the data.
#' @param y A vector with values of dependent variable (outcome).
#' @param x A data.frame of covariates.
#' @param reg.family A regression family. 
#' Can be either "binomial" or "gaussian."
#' @param verbose Indicates verbosing output. Default: FALSE.
#' @return A list of two: "S" - a dataframe with predictors and "fit"
#' - an object returned by "glm" function.
build.null.model <- function(y, x, reg.family="binomial", verbose=FALSE) {
    if(reg.family == "binomial") {
        if(length(x) != 0) {
            fit <- glm(y ~ ., data=data.frame(x), 
                           na.action=na.exclude,
                         family = binomial) 
            #family = binomial(link=logit)),TRUE)
        } else {
            fit <- glm(y ~ 1, na.action=na.exclude,
                         family = binomial) 
            #family = binomial(link=logit)),TRUE)
        }
    } else if(reg.family == "gaussian") {
        if(length(x) != 0) {
            fit <- glm(y ~ ., data=data.frame(x), 
                   na.action=na.exclude,
                   family = reg.family)
        } else {
            fit <- glm(y ~ 1, 
                     na.action=na.exclude,
                     family = reg.family)
        }
    } else {
        stop(paste("Unknown reg.family:", reg.family))
    }
  
    return(resid(fit))
}


preprocessPCA <- function(data, scaleData, cumvar.threshold, verbose) {
    ct <- cumvar.threshold
    res.pca <- prcomp(data, scale=scaleData)
    # Eigenvalues
    eig <- (res.pca$sdev)^2
    # Variances in percentage
    variance <- eig*100/sum(eig)
    # Cumulative variances
    cumvar <- cumsum(variance)
    eig.table <- data.frame(eig = eig, 
                          variance = variance, 
                          cumvariance = cumvar)
  
    ######### Filtering by threshold ##############
    if(length(eig.table$cumvar[eig.table$cumvar <= ct]) == 0) {
        if(verbose) {  
            print("Warning: cumvar.threshold will be set to")
            print(paste("the first component of cum. var.:", 
                  eig.table$cumvar[1]))
        }
        cumvar.threshold <- eig.table$cumvar[1]
    }
  
    #S <- res.pca$x[,which(eig.decathlon2.active$cumvar <= 
    #    cumvar.threshold)] %*% 
    #t(res.pca$rotation[,which(eig.decathlon2.active$cumvar <= 
    #    cumvar.threshold)])
    ct <- cumvar.threshold
    S <- res.pca$x[,which(eig.table$cumvar <= ct)]
  
    ### And add the center (and re-scale) back to data ###
    #if(scale){
    #   S <- scale(S, center = FALSE , scale=1/res.pca$scale)
    #}
    #if(center){
    #  S <- scale(S, center = -1 * res.pca$center, scale=FALSE)
    #}
  
    indexes <- which(eig.table$cumvar <= ct)
    
    return(list(S = S, indexes=indexes, model="PCA"))
}


preprocessPLS <- function(data, pheno, scaleData, cumvar.threshold, out.type) {
  
    
  
    if(scaleData){
        data.scaled <- scale(data, center = TRUE)
    } else {
        data.scaled <- data
    }
  
    n.snp <- dim(data.scaled)[2]
    
    #if(!is.null(cumvar.threshold)) {
    #    ct <- cumvar.threshold/100
    #    npred <- round(n.snp*ct)
    #    numcomp <- ifelse(n.snp < 10, n.snp, npred)
    #}
  
    if(out.type == "D") { # PLS-DA
        stop("Sorry, PLS for dichotomous trait is not yet supported.")
        ## Uncomment the code below if you want to use ropls and PLS-DA ##
        #model <- opls(x = data.scaled, y=as.factor(pheno), 
        #                  predI=numcomp, 
        #                  plotL = FALSE, 
        #                  log10L=FALSE, 
        #                  algoC = "nipals", 
        #                  silent = TRUE)
        #
        #if(inherits(model, "try-error")) {
        #    for(i in 2:6) {
        #        ct <- ct/i
        #        npred <- round(dim(data.scaled)[2]*ct)
        #
        #        model <- opls(x = data.scaled, y=as.factor(pheno), 
        #                      predI=npred, 
        #                      plotL = FALSE, 
        #                      log10L=FALSE, 
        #                      algoC = "nipals", 
        #                      silent = TRUE)
        #
        #        if(!inherits(model, "try-error")) {
        #            break
        #        }
        #    }
        #}
        
        #S <- model@scoreMN #%*% t(model@loadingMN)
        #Y <- model@uMN %*% t(model@cMN)
        
    } else if(out.type == "C") { # PLS
        #model <- opls(x = data.scaled, y=pheno, 
        #              predI=numcomp, 
        #              plotL = FALSE, 
        #              log10L=FALSE, 
        #              algoC = "nipals", 
        #              silent = TRUE)
    
        #if(inherits(model, "try-error")) {
        #    cumvar.threshold <- cumvar.threshold/2
        #    npred <- round(dim(data.scaled)[2]*ct)
        #    model <- opls(x = data.scaled, y=as.factor(pheno), 
        #                  predI=npred, plotL = FALSE, 
        #                  log10L=FALSE, algoC = "nipals", 
        #                  silent = TRUE)
        #}
        
        dd <- cbind(pheno, data.scaled)
        temp.model <- plsr(pheno ~ ., data=dd, validation = "CV")
        if(is.null(cumvar.threshold)) {
            numcomp <- selectNcomp(temp.model)
            if(numcomp == 0) numcomp <- 1
        } else {
            ct <- cumvar.threshold/100
            numcomp <- round(n.snp*ct)
            
            #bbb <- sort(explvar(temp.model), decreasing = TRUE)
            #tt <- bbb[1]
            #if(tt >= cumvar.threshold) {
            #    numcomp = 1
            #} else {
            #    for(i in 2:length(bbb)) {
            #        tt <- tt + bbb[i]
            #        if(tt >= cumvar.threshold) {
            #            numcomp <- i
            #            break
            #        }
            #    }
            #    numcomp <- i
            #}
            
        }
        
        if(numcomp == 1) numcomp <- 2
        
        model <- plsr(pheno ~ ., numcomp, data=dd, validation = "none")
        
        # X
        nrow.scores <- dim(model[["scores"]])[1]
        ncol.scores <- dim(model[["scores"]])[2]
        nrow.loadings <- dim(model[["loadings"]])[1]
        ncol.loadings <- dim(model[["loadings"]])[2]
        # Y
        nrow.yscores <- dim(model[["Yscores"]])[1]
        ncol.yscores <- dim(model[["Yscores"]])[2]
        nrow.yloadings <- dim(model[["Yloadings"]])[1]
        ncol.yloadings <- dim(model[["Yloadings"]])[2]
        # , byrow = TRUE
        S <- matrix(model[["scores"]], nrow=nrow.scores, ncol=ncol.scores) #%*% t(matrix(model[["loadings"]], nrow=nrow.loadings, ncol=ncol.loadings))
        Y <- matrix(model[["Yscores"]], nrow=nrow.yscores, ncol=ncol.yscores) %*% t(matrix(model[["Yloadings"]], nrow=nrow.yloadings, ncol=ncol.yloadings))
        
    }
  
    
    
    return(list(S = S, Y = Y, model=model))
}

preprocessLASSO <- function(data, pheno, reg.family, penalty=0.001) {
    #### LASSO ####
    tryCatch({
        if(reg.family == "binomial") {
            pheno <- as.factor(pheno)
        }
        fit <- glmnet(x=as.matrix(data),
                         alpha=1, # LASSO
                         y=pheno, 
                         family=reg.family)
        #fit <- cv.glmnet(x=as.matrix(data),
        #              alpha=1, # LASSO
        #              y=pheno, 
        #              family=reg.family,
        #              nfolds=10)
    
    }, error=function(e) {
        print(e)
    })
    S <- data[,which(coef(fit, s=penalty)[-1] != 0)]
    return(list(S=S, fit=fit, model="LASSO"))
}

preprocessRidge <- function(data, pheno, reg.family, penalty=0.001) {
    tryCatch({
        if(reg.family == "binomial") {
            pheno <- as.factor(pheno)
        }
        fit <- glmnet(x=as.matrix(data),
                         alpha=0, # Ridge
                         y=pheno, 
                         family=reg.family)
    }, error=function(e) {
        print(e)
    })
    S <- data[,which(coef(fit, s=penalty)[-1] != 0)]
    return(list(S=S, fit=fit, model="ridge"))
}