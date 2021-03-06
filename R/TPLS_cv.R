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
#' # Fit example TPLS data with a TPLS model using cross-validation
#' # Load example data (included with package).
#' X = TPLSdat$X # single trial brain image of subjects pressing left/right buttons
#' Y = TPLSdat$Y # binary variable that is 1 if right button is pushed, 0 if left button is pushed
#' subj = TPLSdat$subj # 1, 2, or 3, depending on who the subject is
#' run = TPLSdat$run # 1, 2, ..., 8, depending on the scan run of each subject
#'
#' # Fit the model, using 3-fold cross-validation at the subject level
#' # (i.e., train on two subjects, test on 1, repeat three times)
#' TPLScvmdl <- TPLS_cv(X,Y,subj)
#'
#' # Evaluate the tuning parameters via cross-validation.
#' # We'll test 1~50 components and thresholding from 0 to 1 in 0.05 increments.
#' # Also include subfold information.
#' # This allows for calculation of correlation at the run-level instead of at the subject level.
#' cvstats <- evalTuningParam(TPLScvmdl,"pearson",X,Y,1:50,seq(0,1,0.05),subfold=run)
#'
#' # plot the tuning parameter surface.
#' # It'll show the point of best performance (and also point of 1SE performance).
#' # The plot is interactive, so spin it around
#' plotTuningSurface(cvstats)
#'
#' # These are the tuning parameters of best performance
#' cvstats$compval_best # 8 components
#' cvstats$threshval_best # 0.1 thresholding (leave only 10% of all voxels)
#'
#' # Now build a new TPLS model, using all the data, using the best tuning parameters
#' TPLSmdl <- TPLS(X,Y,NComp=cvstats$compval_best)
#'
#' # Extract the prediction betamap that gave the best CV performance
#' betamap <- makePredictor(TPLSmdl,cvstats$compval_best,cvstats$threshval_best)
#'
#' # This is the intercept
#' betamap$bias
#'
#' # These are the coefficients for the original variables
#' betamap$betamap
#'
#' # Project the betamap into brain-space so that we can look at it.
#' mask = TPLSdat$mask # mask 3D image of the brain from which X was extracted from
#' brainimg = mask*1 # make a copy
#' brainimg[mask] = betamap$betamap # put the betamap into the brain image
#' fig1 <- plot_ly(z = brainimg[,15,], type = "heatmap") # looking at a slice of the brain image
#' fig2 <- plot_ly(z = 1*mask[,15,], type = "heatmap") # a slice of the brain mask for reference
#' fig <- subplot(fig1, fig2)
#' fig
#'
#' # Figures show a coronal section of the brain (but flipped right 90 degrees).
#' # on the left, you should see the bilateral motor cortex coefficients with opposing signs.
#' # This is just a simple visual demonstration. You should use other packages to output
#' # coefficients into a nifti file and view them in a separate viewer.
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
