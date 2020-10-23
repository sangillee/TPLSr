#' Fit a TPLS model to data with cross validation
#'
#' @param X n-by-v data matrix of real numbers. Rows correspond to observations (trials) and columns to variables (e.g., fMRI voxels).
#' @param Y n-by-1 Vector of real numbers. Can be binary (0/1) for classification model, or can be continuous.
#' @param foldid A vector of values between 1 and number of folds identifying what fold each observation is in.
#' @param NComp Maximum number of partial least squares component you want to use. Default is 50, and this is on the safe side for fMRI.
#' @param W n-by-1 vector of positive observation weights.
#' @return A TPLS_cv object that contains the following attributes. Most of the time, you won't need to access the attributes.
#' \itemize{
#'     \item \code{NComp}: The number of components you specified in the input
#'     \item \code{numfold}: Total number of cross-validation folds
#'     \item \code{testfold}: A vector indice that should be the same as foldid, if it was provided accurately.
#'     \item \code{cvMdls} : A vector of TPLS models, one for each fold.
#' }
#' @examples
#' # Fit TPLS model to example data
#' @export

TPLS_cv <- function(X,Y,foldid,NComp=50,W=0){
  # input handling
  if(any(W < 0)){stop("Observation weights should be non-negative")}
  if(length(W) == 1){W <- matrix(1, nrow=length(Y), ncol=1)} # if weight is not provided, every observation has equal weight

  # identifying folds and preparing structure
  uniqfold = unique(foldid)
  numfold = length(uniqfold)
  testfold = foldid
  cvMdls <- matrix(list(), nrow=numfold, ncol=1)
  for (i in 1:numfold){
    message(paste("Fold #",i))
    train = foldid != uniqfold[i]
    testfold[!train] = i
    cvMdls[[i]] = TPLS(X[train,],Y[train],NComp,W[train])
  }

  TPLScvmdl <- list("NComp" = NComp, "numfold" = numfold, "testfold" = testfold, "cvMdls" = cvMdls)
  TPLScvmdl
}
