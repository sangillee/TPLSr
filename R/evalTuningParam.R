#' Evaluating cross-validation performance of a TPLS_cv model at compvec and threshvec
#'
#' @param TPLScvmdl TPLS_cv model created from \code{TPLS_cv}
#' @param type CV performance metric type. One of LLbinary, negMSE, Pearson, Spearman, AUC, ACC.
#' @param X The same X as used in \code{TPLScvmdl}.
#' @param Y The SAME Y as used in \code{TPLScvmdl}.
#' @param compvec Vector of number of components to test in cross-validation.
#' @param threshvec Vector of threshold level (0 ~ 1) to test in cross-validation.
#' @param subfold (Optional) vector of subdivision within testing fold to calculate performance. For example scan run division within subject.
#'
#' @return A evalTuningParam object that contains the following attributes.
#' \itemize{
#'     \item \code{type}: Cross validation performance measure type, as specified in the input
#'     \item \code{threshval}: Same as the input threshvec
#'     \item \code{compval}: Same as the input compvec
#'     \item \code{perfmat}: Performance measure 3D matrix: length(compvec)-by-length(threshvec)-by-numfold
#'     \item \code{perf_best}: Best CV performance out of all combinations of compvec and threshvec
#'     \item \code{compval_best}: Number of components that gave the best performance (i.e., perf_best)
#'     \item \code{threshval_best}: Threshold level that gave the best performance (i.e., perf_best)
#'     \item \code{perf_1se} : Performance of the most parsimonious model (least number of coefficients) that is within 1 standard error of perf_best.
#'     \item \code{compval_1se} : Number of components that gave perf_1se
#'     \item \code{threshval_1se} : Threshold level that gave perf_1se
#'     \item \code{best_at_threshold} : a 3-column matrix; first column is max performance at threshold, second column is threshold values, third column is number of components for the best model at threshold
#' }
#' @import plotly
#' @export

evalTuningParam <- function(TPLScvmdl,type=c("Pearson","negMSE","ACC","AUC","LLbinary","Spearman"),X,Y,compvec,threshvec,subfold=NULL){
  # input checking
  if(is.null(subfold)){subfold = rep(1,length(Y))}
  X = as.matrix(X); TPLSinputchecker(X,"X","mat",NULL,NULL,1)
  Y = as.matrix(Y); TPLSinputchecker(Y,"Y","colvec",NULL,NULL,1)
  TPLSinputchecker(compvec,"compvec","vec",TPLScvmdl$NComp,1,0,1); compvec = sort(compvec);
  TPLSinputchecker(threshvec,"threshvec","vec",1,0); threshvec = sort(threshvec);
  TPLSinputchecker(subfold,"subfold",'vec')

  # Perform CV prediction and measure performance
  perfmat = array(NA,dim=c(length(compvec),length(threshvec),TPLScvmdl$numfold))
  for (i in 1:TPLScvmdl$numfold){
    message(paste("Fold #",i))
    testCVfold = TPLScvmdl$CVfold[,i] == 1
    Ytest = Y[testCVfold]
    testsubfold = subfold[testCVfold]
    uniqtestsubfold = unique(testsubfold)
    for (j in 1:length(threshvec)){
      predmat = TPLSpredict(TPLScvmdl$cvMdls[[i]],compvec,threshvec[j],X[testCVfold,])
      smallperfmat = matrix(nrow=length(compvec),ncol=length(uniqtestsubfold))
      for (k in 1:length(uniqtestsubfold)){
        subfoldsel = testsubfold == uniqtestsubfold[k]
        smallperfmat[,k] = util_perfmetric(predmat[subfoldsel,],Ytest[subfoldsel],type)
      }
      perfmat[,j,i] = rowMeans(smallperfmat)
    }
  }

  # find the point of maximum CV performance
  avgperfmat = apply(perfmat,c(1,2),mean) # mean performance
  perf_best = max(avgperfmat,na.rm = TRUE) # best mean performance
  bestind = which(avgperfmat==perf_best,TRUE) # coordinates of best point
  compval_best = compvec[bestind[1,1]]; threshval_best = threshvec[bestind[1,2]] # component and threshold of best point
  standardError = stats::sd(perfmat[bestind[1,1],bestind[1,2],])/sqrt(dim(perfmat)[3]) # standard error of best point
  ind1se = which(avgperfmat[,1:bestind[1,2],drop = FALSE]>=(perf_best-standardError),TRUE) # coordinates of 1SE point
  perf_1se = avgperfmat[ind1se[1,1],ind1se[1,2]] # performance of 1SE point
  compval_1se = compvec[ind1se[1,1]]; threshval_1se = threshvec[ind1se[1,2]]
  best_at_threshold = matrix(NA,nrow=length(threshvec),ncol=3)
  for (i in 1:length(threshvec)) {
    tempmax = max(avgperfmat[,i]); tempind = which(avgperfmat[,i]==tempmax,TRUE)
    best_at_threshold[i,] = c(tempmax,threshvec[i],compvec[tempind[1]])
  }

  # output list
  cvstats <- list("type" = type, "threshval" = threshvec, "compval" = compvec, "perfmat" = perfmat,
              "perf_best"= perf_best, "compval_best"=compval_best, "threshval_best" = threshval_best,
              "perf_1se"= perf_1se, "compval_1se"=compval_1se, "threshval_1se"=threshval_1se,
              "best_at_threshold"=best_at_threshold)
  class(cvstats) <- "evalTuningParam"
  cvstats
}


#' Plots the tuning surface of TPLS
#'
#' @param object : evalTuningParam object
#' @export

plotTuningSurface <- function(object){
  meansurf = apply(object$perfmat,c(1,2),mean)
  fig <- mygrid(object$threshval,object$compval,meansurf)
  axy <- list(title = "Number of PLS components"); axx <- list(title = "Proportion of Voxels Left"); axz <- list(title = object$type)
  fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  fig <- fig %>% add_markers(object$threshval_best,object$compval_best,object$perf_best, name="Max Perf")
  fig <- fig %>% add_markers(object$threshval_1se,object$compval_1se,object$perf_1se, name ="1SE Perf")
  fig
}

util_perfmetric <- function(predmat,testY,type){
  switch(type,
         LLbinary={
           if(binarycheck(testY)!=1){stop("LL binary can be only calculated for binary measures")}
           predmat[testY!=1,] = 1-predmat[testY!=1,] # flip probability
           predmat[predmat>1] = 1; predmat[predmat<=0] = .Machine$double.xmin # take care of probability predictions outside of range
           Perf = colMeans(log(predmat)) # taking the mean so that we have per-trial average LL, which would make each fold count equally (not weighted by number of trials)
         },
         negMSE={
           Perf = -colMeans((predmat-testY)^2)
         },
         ACC={
           if(binarycheck(testY)!=1){stop("Accuracy can be only calculated for binary measures")}
           Perf = colMeans( (1*(predmat>0.5))==testY )
         },
         AUC={
           if(binarycheck(testY)!=1){stop("AUC can be only calculated for binary measures")}
           n = length(testY); num_pos = sum(testY==1); num_neg = n - num_pos
           if (num_pos>0 && num_pos < n){
             ranks = apply(predmat,2,rank)
             Perf = (  colSums(ranks[testY==1,]) - num_pos*(num_pos+1)/2 ) / (num_pos*num_neg)
           }else{
             Perf = matrix(0.5,1,ncol(predmat))
           }
         },
         Pearson={
           suppressWarnings({result = stats::cor(predmat, testY, method = "pearson")})
           result[is.na(result)] = 0 # this probably happened because your prediction (or your Y) has no variation.
           return(result)
         },
         Spearman={
           suppressWarnings({result = stats::cor(predmat, testY, method = "spearman")})
           result[is.na(result)] = 0 # this probably happened because your prediction (or your Y) has no variation.
           return(result)
         }
  )
}

binarycheck <- function(Y){
  if(any(Y!=0 & Y!=1)){return(0)}else{return(1)}
}

mygrid <- function(x,y,z){
  # stupid code just to add a grid to 3D surface! Hope there's a better way
  fig <- plot_ly()
  fig <- fig %>% add_surface(x = x, y = y, z = z,lighting = list(roughness = 0.9))
  gridx <- c(); gridy <- c(); gridz <- c()
  for (i in 1:nrow(z)){
    if(i %% 2 == 1){
      gridx <- c(gridx,x); gridz <- c(gridz,z[i,])
    }else{
      gridx <- c(gridx,rev(x)); gridz <- c(gridz,rev(z[i,]))
    }
    gridy <- c(gridy,rep.int(y[i],length(x)))
  }
  fig <- fig %>% add_trace(x = gridx, y = gridy, z = gridz,type = 'scatter3d', mode = 'lines',line = list(width = 3, color = 'rgb(120, 120, 120)'), showlegend = FALSE)
  gridx <- c(); gridy <- c(); gridz <- c()
  for (i in 1:ncol(z)){
    if(i %% 2 == 1){
      gridy <- c(gridy,y); gridz <- c(gridz,z[,i])
    }else{
      gridy <- c(gridy,rev(y)); gridz <- c(gridz,rev(z[,i]))
    }
    gridx <- c(gridx,rep.int(x[i],length(y)))
  }
  fig <- fig %>% add_trace(x = gridx, y = gridy, z = gridz,type = 'scatter3d', mode = 'lines',line = list(width = 3, color = 'rgb(120, 120, 120)'), showlegend = FALSE)
  fig <- fig %>% add_trace(x = ~x, y = ~y[max.col(t(z))], z = ~apply(z,2,max),type = 'scatter3d', mode = 'lines+markers',
                           line = list(width = 6, color = 'rgb(50, 200, 255)'), marker = list(size = 3.5, color = 'rgb(50,50,50)'), name = 'best at threshold')
  return(fig)
}
