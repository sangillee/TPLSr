#' Fit a TPLS model to data
#'
#' @param X n-by-v data matrix of real numbers. Rows correspond to observations (trials) and columns to variables (e.g., fMRI voxels).
#' @param Y n-by-1 Vector of real numbers. Can be binary (0/1) for classification model, or can be continuous.
#' @param NComp Maximum number of partial least squares component you want to use. Default is 50, and this is on the safe side for fMRI.
#' @param W n-by-1 vector of positive observation weights.
#' @param nmc A switch to skip mean-centering. Default is off (0). Only turn it on (1) when the data is already mean-centered and you want to save memory space by not creating another copy of the data for mean-centering.
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
#' @examples
#' # Fit TPLS model to example data
#' @export

TPLS <- function(X,Y,NComp=50,W=0,nmc=0){
  if(any(W < 0)){stop("Observation weights should be non-negative")}
  if(length(W) == 1){W <- matrix(1, nrow=length(Y), ncol=1)} # if weight is not provided, every observation has equal weight
  W = W/sum(W) # normalize weight sum to 1
  n = nrow(X); v = ncol(X)

  # would like to potentially include support for single precision computing, but not sure how to yet...

  # Mean-Center variables as needed by SIMPLS algorithm
  MtrainX = t(W)%*%X; MtrainY = c(t(W)%*%Y) # Weighted means of X and Y
  if(nmc==0){X = X-rep(MtrainX,rep.int(n,v)); Y = Y-MtrainY # if no switch is given to skip mean centering
  }else{
    message('mean centering disabled')
    if(mean(abs(MtrainX)) > 1e-04){warning("X does not seem to be mean-centered. Results may not be valid")}
  }

  # allocate memories for output variables, interim variables, and calculate often used variables
  scoreCorr = matrix(nrow=NComp,ncol=1); betamap = matrix(nrow=v,ncol=NComp); threshmap = cbind(matrix(0.5,nrow=v,ncol=1),matrix(nrow=v,ncol=NComp-1))
  B = matrix(nrow=NComp,ncol=1); P2 = matrix(nrow=n,ncol=NComp); C = matrix(nrow=v,ncol=NComp); sumC2 = matrix(0,nrow=v,ncol=1); r = Y; V = matrix(nrow=v,ncol=NComp)
  WT = t(W); WYT = t(W*Y); WTY2 = c(WT%*%(Y^2)); W2 = W^2 # often-used variables

  # perform Arthur-modified SIMPLS algorithm
  Cov = t(WYT%*%X) # weighted covariance
  for (i in 1:NComp){
    message(paste("Calculating Comp #",i))
    P = X%*%Cov # this is the component, before normalization
    norm_P = c(sqrt(WT%*%(P^2))) # weighted standard deviation of component as normalization constant
    P = P / norm_P; B[i] = (norm(Cov,type="2")^2)/norm_P; C[,i] = Cov / norm_P # normalize component, beta, and back-projection coefficient

    # Update the orthonormal basis with modified Gram Schmidt
    vi = t(t(W*P)%*%X) #weighted covariance between X and current component
    vi = vi - V[,1:i-1]%*%(t(V[,1:i-1])%*%vi) # orthogonalize with regards to previous components
    vi = vi / norm(vi,type="2"); V[,i] = vi; # normalize and add to orthonormal basis matrix
    Cov = Cov - vi%*%(t(vi)%*%Cov); Cov = Cov - V[,1:i]%*%(t(V[,1:i])%*%Cov) # Deflate Covariance using the orthonormal basis matrix

    # back-projection
    sumC2 = sumC2 + C[,i]^2; P2[,i] = P^2; r = r - P*B[i] # some variables that will facilitate computation later
    if(i>1){ # not the first component
      betamap[,i] = C[,1:i]%*%B[1:i] # back-projection of coefficients
      se = sqrt( t(P2[,1:i])%*%(W2*(r^2))  ) # Huber-White Sandwich estimator (assume no small T bias)
      absz = abs( (C[,1:i]%*%(B[1:i]/se)) / sqrt(sumC2) ) # absolute value of back-projected z-statistics
      threshmap[,i] = (v- rank(absz) )/v # convert into thresholds between 0 and 1
    }else{betamap[,i] = C[,1:i]*B[1:i]} # back-projection of coefficients
  }

  pctVar = (B^2) / WTY2 # Compute the percent of variance of Y each component explains
  scoreCorr = sqrt(pctVar) # weighted correlation between Y and current component


  TPLSmdl <- list("NComp" = NComp, "W" = W, "MtrainX" = MtrainX, "MtrainY" = MtrainY, "scoreCorr"= scoreCorr, "pctVar"=pctVar, "betamap" = betamap, "threshmap"= threshmap)
  class(TPLSmdl) <- "TPLS"
  TPLSmdl
}

# dummy generic
makePredictor <- function(object,compval,threshval){
  UseMethod("makePredictor")
}

#' Extracts a predictor (betamap and intercept) from a TPLS model at a given number of components and given threshold value
#'
#' @param TPLSmdl A TPLS object created from using function \code{TPLS}
#' @param compval The number of components you want in your model. Providing a vector will provide multiple betamaps (e.g., c(3,4,5) will provide three betamaps each with 3, 4, and 5 PLS components)
#' @param threshval Threshold number between 0 and 1 (inclusive) for thresholding the betamap. This must be a scalar.
#' @return
#' \itemize{
#'     \item \code{bias}: The intercept of the extracted model. Vector of intercepts if compval is a vector.
#'     \item \code{betamap}: Column vector of betamap. Matrix of betamaps if compval is a vector.
#' }
#' @examples
#' # See examples for TPLS
#' @export

makePredictor.TPLS <- function(TPLSmdl,compval,threshval){
  if(length(threshval)>1){stop("only one threshold value should be used")}
  betamap = TPLSmdl$betamap[,compval]
  if(threshval<1){
    betamap = betamap*(TPLSmdl$threshmap[,compval]<=threshval)
  }
  bias = TPLSmdl$MtrainY - TPLSmdl$MtrainX %*% betamap # post-fitting of bias
  return(list("bias" = bias, "betamap" = betamap))
}


#' Make predictions about given data \code{testX} by using an extracted TPLSmdl with \code{compval} components and \code{threshval} threshold.
#'
#' @param TPLSmdl A TPLS object created from using function \code{TPLS}
#' @param compval The number of components you want in your model. Providing a vector will provide multiple predictions (e.g., c(3,4,5) will provide three prediction columns each with 3, 4, and 5 PLS components)
#' @param threshval Threshold number between 0 and 1 (inclusive) for thresholding the betamap. This must be a scalar.
#' @return
#' \itemize{
#'     \item \code{score}: Column vector of prediction scores. Matrix of scores if compval is a vector.
#' }
#' @examples
#' # See examples for TPLS
#' @export

predict.TPLS <- function(TPLSmdl,compval,threshval,testX){
  if(length(threshval)>1){stop("only one threshold value should be used")}
  if(threshval == 0){
    score = matrix(TPLSmdl$MtrainY,nrow=nrow(testX),ncol=length(compval))
  }else{
    preds = makePredictor(TPLSmdl,compval,threshval)
    score = matrix(preds$bias,nrow=nrow(testX),ncol=length(compval)) + testX%*%preds$betamap
  }
  return(score)
}
