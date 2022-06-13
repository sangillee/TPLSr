#' Constructor method for fitting a cross-validation T-PLS model
#'
#' @param X Numerical matrix of predictors. Typically single-trial betas where each column is a voxel and row is observation
#' @param Y Variable to predict. Binary 0 and 1 in case of classification, continuous variable in case of regression
#' @param CVfold Cross-validation testing fold information. Can either be a vector or a matrix, the latter being more general.
#'               Vector: n-by-1 vector. Each element is a number ranging from 1 ~ numfold to identify which testing fold each observation belongs to
#'               Matrix: n-by-numfold matrix. Each column indicates the testing data with 1 and training data as 0.
#'               Example: For leave-one-out CV, Vector would be 1:n, Matrix form would be eye(n)
#'               Matrix form is more general as it can have same trial be in multiple test folds
#' @param NComp (Optional) Number of PLS components to compute. Default is 25.
#' @param W (Optional) Observation weights. Optional input. By default, all observations have equal weight.
#'          Can either be a n-by-1 vector or a n-by-nfold matrix where each column is observation weights in that CV fold
#' @param nmc (Optional) 'no mean centering'. See TPLS for more detail.
#'             Turning this on will skip mean centering on all cross validation folds, so they should all be mean-centered already
#' @return A TPLS_cv object that contains the following attributes. Most of the time, you won't need to access the attributes.
#' \itemize{
#'     \item \code{NComp}: The number of components you specified in the input
#'     \item \code{numfold}: Total number of cross-validation folds
#'     \item \code{CVfold}: A matrix of indicators for testing data for each cross validation fold in each column
#'     \item \code{cvMdls} : A vector of TPLS models, one for each fold.
#' }
#' See vignettes for tutorial
#' @export

TPLS_cv <- function(X,Y,CVfold,NComp=25,W=NULL,nmc=0){
  # input checking
  X = as.matrix(X); TPLSinputchecker(X,"X","mat",NULL,NULL,1)
  Y = as.matrix(Y); TPLSinputchecker(Y,"Y","colvec",NULL,NULL,1)
  TPLSinputchecker(CVfold,"CVfold")
  TPLSinputchecker(NComp,"NComp","scalar",NULL,1,0,1)
  if(is.null(W)){W <- matrix(1, nrow=length(Y), ncol=1)} # if weight is not provided, every observation has equal weight
  else{W= as.matrix(W);TPLSinputchecker(W,"W",NULL,NULL,0)};
  TPLSinputchecker(nmc,"nmc","scalar")
  if(nrow(Y)!=nrow(X) || nrow(W)!=nrow(X)){stop("X, Y, and W should have equal number of rows")}
  CVfold = prepCVfold(CVfold); numfold = ncol(CVfold) # convert CVfold into matrix form, if not already
  if(ncol(W)==1){W = W[,rep(1,numfold)]} # convert weights into matrix form, if not already

  cvMdls <- matrix(list(), nrow=numfold, ncol=1)
  for (i in 1:numfold){
    message(paste("Fold #",i))
    train = CVfold[,i] == 0
    cvMdls[[i]] = TPLS(X[train,],Y[train],NComp,W[train,i],nmc)
  }

  TPLScvmdl <- list("NComp" = NComp, "numfold" = numfold, "CVfold" = CVfold, "cvMdls" = cvMdls)
  class(TPLScvmdl) <- "TPLS_cv"
  TPLScvmdl
}



#' prepare CV fold data into a matrix form, which is more generalizable
#' @noRd

prepCVfold <- function(inCVfold){
  if(is.vector(inCVfold)){
    return(prepCVfold(as.matrix(inCVfold)))
  }

  if(ncol(inCVfold)==1){ # single column matrix
    uniqfold = unique(inCVfold); nfold = length(uniqfold);
    CVfold = matrix(0,length(inCVfold),nfold)
    for (i in 1:nfold) {
      CVfold[,i] = 1*(inCVfold == uniqfold[i])
    }
  }else{ # matrix format
    nfold = ncol(inCVfold); CVfold = inCVfold;
    if(any(CVfold !=0 & CVfold != 1)){
      stop("Matrix form CVfold should only have binary elements 0 and 1 with each column noting the test observations in 1 and training observations in 0")
    }
  }
  return(CVfold)
}
