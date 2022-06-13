#' Constructor method for fitting a T-PLS model with given data X and Y.
#'
#' @param X Numerical matrix of predictors. Typically single-trial betas where each column is a voxel and row is observation
#' @param Y Variable to predict. Binary 0 and 1 in case of classification, continuous variable in case of regression
#' @param NComp (Optional) Number of PLS components to compute. Default is 25.
#' @param W (Optional) Observation weights. By default, all observations have equal weight.
#' @param nmc (Optional) 'no mean centering'. Default is 0. If 1, T-PLS will skip mean-centering. This option is only provided in case you already mean-centered the data and want to save some memory usage.
#'
#' @return A TPLS object that contains the following attributes. Most of the time, you won't need to access the attributes.
#' \itemize{
#'     \item \code{NComp}: The number of components you specified in the input
#'     \item \code{W}: Normalized version of the observation weights (i.e., they sum to 1)
#'     \item \code{MtrainX}: Column mean of X. Weighted mean if W is given.
#'     \item \code{MtrainY}: Mean of Y. Weighted mean if W is given.
#'     \item \code{scoreCorr}: Correlation between Y and each PLS component. Weighted correlation if W is given.
#'     \item \code{pctVar}: Proportion of variance of Y that each component explains.
#'     \item \code{betamap}: v-by-NComp matrix of TPLS coefficients for each of the v variables, provided at each model with NComp components.
#'     \item \code{threshmap} : v-by-NComp matrix of TPLS threshold values (0~1) for each of the v variables, provided at each model with NComp components.
#' }
#' See vignettes for tutorial
#' @export

TPLS <- function(X,Y,NComp=25,W=NULL,nmc=0){
  # input checking
  X = as.matrix(X); TPLSinputchecker(X,"X","mat",NULL,NULL,1); n = nrow(X); v = ncol(X)
  Y = as.matrix(Y); TPLSinputchecker(Y,"Y","colvec",NULL,NULL,1)
  TPLSinputchecker(NComp,"NComp","scalar",NULL,1,0,1)
  if(is.null(W)){W <- matrix(1, nrow=length(Y), ncol=1)} # if weight is not provided, every observation has equal weight
  else{W= as.matrix(W);TPLSinputchecker(W,"W","colvec",NULL,0)};
  TPLSinputchecker(nmc,"nmc","scalar")
  if(nrow(Y)!=n || nrow(W)!=n){stop("X, Y, and W should have equal number of rows")}
  W = W/sum(W) # normalize weight sum to 1

  # Mean-Center variables as needed by SIMPLS algorithm
  MtrainX = t(W)%*%X; MtrainY = c(t(W)%*%Y) # calculating weighted means of X and Y
  if(nmc==0){ # do mean centering
    X = X-rep(MtrainX,rep.int(n,v)); Y = Y-MtrainY # subtract means
  }else{
    if(any(abs(MtrainX) > 1e-04)){warning("Skipped mean centering, but X does not seem to be mean-centered. Results may be invalid")}
  }

  # allocate memories
  pctVar = matrix(nrow=NComp,ncol=1); scoreCorr = matrix(nrow=NComp,ncol=1) # percent of variance of Y each component explains, weighted correlation between Y and current component
  betamap = matrix(nrow=v,ncol=NComp); threshmap = matrix(0.5,nrow=v,ncol=NComp); zmap = matrix(nrow=v,ncol=NComp) # output variables
  B = matrix(nrow=NComp,ncol=1); P2 = matrix(nrow=n,ncol=NComp); C = matrix(nrow=v,ncol=NComp); sumC2 = matrix(0,nrow=v,ncol=1); r = Y; V = matrix(nrow=v,ncol=NComp) #interim variables
  WT = t(W); WTY2 = c(WT%*%(Y^2)); W2 = W^2 # often-used variables in calculation

  # perform Arthur-modified SIMPLS algorithm
  Cov = t(t(W*Y)%*%X); normCov = norm(Cov,type="2") # weighted covariance
  for (i in 1:NComp){
    message(paste("Calculating Comp #",i))
    P = X%*%Cov; norm_P = c(sqrt(WT%*%(P^2)))# this is the component and its weighted stdev
    P = P / norm_P; B[i] = (normCov^2)/norm_P; C[,i] = Cov / norm_P # normalize component, beta, and back-projection coefficient
    pctVar[i] = (B[i]^2)/WTY2; scoreCorr[i] = sqrt(pctVar[i]);

    # Update the orthonormal basis with modified Gram Schmidt
    vi = t(t(W*P)%*%X) #weighted covariance between X and current component
    vi = vi - V[,1:i-1]%*%(t(V[,1:i-1])%*%vi) # orthogonalize with regards to previous components
    vi = vi / norm(vi,type="2"); V[,i] = vi; # normalize and add to orthonormal basis matrix
    Cov = Cov - vi%*%(t(vi)%*%Cov); Cov = Cov - V[,1:i]%*%(t(V[,1:i])%*%Cov); normCov = norm(Cov,type="2") # Deflate Covariance using the orthonormal basis matrix

    # back-projection
    sumC2 = sumC2 + C[,i]^2; P2[,i] = P^2; r = r - P*B[i] # some variables that will facilitate computation later
    if(i>1){ # not the first component
      betamap[,i] = C[,1:i]%*%B[1:i] # back-projection of coefficients
      se = sqrt( t(P2[,1:i])%*%(W2*(r^2))  ) # Huber-White Sandwich estimator (assume no small T bias)
      zmap[,i] = (C[,1:i]%*%(B[1:i]/se)) / sqrt(sumC2) # back-projected z-statistics
      threshmap[,i] = (v- rank( abs(zmap[,i]) ) )/v # convert into thresholds between 0 and 1
    }else{betamap[,i] = C[,1:i]*B[1:i]} # back-projection of coefficients

    # check if there's enough covariance to milk
    if(normCov < 10*.Machine$double.eps){
      message("All Covariance between X and Y has been explained. Stopping..."); break
    }else{
      if(pctVar[i] < 10*.Machine$double.eps){ # Proportion of Y variance explained is small
        message("New PLS component does not explain more covariance. Stopping..."); break
      }
    }
  }


  TPLSmdl <- list("NComp" = NComp, "W" = W, "MtrainX" = MtrainX, "MtrainY" = MtrainY, "scoreCorr"= scoreCorr, "pctVar"=pctVar, "betamap" = betamap, "threshmap"= threshmap, "zmap"=zmap)
  class(TPLSmdl) <- "TPLS"
  TPLSmdl
}

#' Method for extracting the T-PLS predictor at a given compval and threshval
#'
#' @param TPLSmdl A TPLS object created from using function \code{TPLS}
#' @param compval Vector of number of components to use in final predictor. Providing a vector will provide multiple betamaps (e.g., c(3,4,5) will provide three betamaps each with 3, 4, and 5 PLS components)
#' @param threshval Threshold number between 0 and 1 (inclusive) for thresholding the betamap. This must be a scalar.
#' @return
#' \itemize{
#'     \item \code{bias}: The intercept of the extracted model. Vector of intercepts if compval is a vector.
#'     \item \code{betamap}: Column vector of betamap. Matrix of betamaps if compval is a vector.
#' }
#' @export

makePredictor <- function(TPLSmdl,compval,threshval){
  TPLSinputchecker(compval,"compval","vec",TPLSmdl$NComp,1,0,1)
  TPLSinputchecker(threshval,"threshval","scalar",1,0)
  if(threshval==0){
    betamap = TPLSmdl$betamap[,compval] * 0
  }else{
    betamap = TPLSmdl$betamap[,compval] * (TPLSmdl$threshmap[,compval]<=threshval)
  }
  bias = TPLSmdl$MtrainY - TPLSmdl$MtrainX %*% betamap # post-fitting of bias
  return(list("bias" = bias, "betamap" = betamap))
}


#' Method for making predictions on a testing dataset testX
#'
#' @param TPLSmdl A TPLS object created from using function \code{TPLS}
#' @param compval Vector of number of components to use in final predictor. Providing a vector will provide multiple predictions (e.g., c(3,4,5) will provide three prediction columns each with 3, 4, and 5 PLS components)
#' @param threshval Threshold number between 0 and 1 (inclusive) for thresholding the betamap. This must be a scalar.
#' @param testX Data that you want to predict the Y of
#' @return
#' \itemize{
#'     \item \code{score}: Column vector of prediction scores. Matrix of scores if compval is a vector.
#' }
#' @export

TPLSpredict <- function(TPLSmdl,compval,threshval,testX){
  TPLSinputchecker(testX,"testX")
  preds = makePredictor(TPLSmdl,compval,threshval)
  score = do.call(rbind,replicate(nrow(testX),preds$bias,simplify=FALSE)) + testX%*%preds$betamap
  return(score)
}
