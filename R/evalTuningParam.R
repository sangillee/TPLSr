#' Evaluate TPLS tuning parameters using cross validation
#'
#' @param TPLScvmdl TPLS_cv model created from \code{TPLS_cv}
#' @param type Cross validation performance measure type. One of 'pearson', 'spearman', or 'AUC'
#' @param X The SAME X that was used to create the \code{TPLScvmdl}. If it's not the same, the function may not work or the results will be completely off
#' @param Y The SAME Y that was used to create the \code{TPLScvmdl}.
#' @param compvec Vector containing the number of components you want to assess CV performance for (e.g., c(3,4,5) will provide CV performance of 3, 4, and 5 component TPLS model at various thresholds)
#' @param threshvec Vector containing the thresholding level betweeon 0 and 1 you want to assess CV performance for (e.g., seq(0,1,0.1) will provide CV performance of TPLS models at thresholds of 0, 0.1, 0.2, ... ,1)
#' @param subfold Optional vector containing smaller data division within folds. For example, if the cross-validation was done at the subject level, with each testing fold being a subject,
#' subfold can be the run number of the scan of each person. This allows for calculation of average CV metric at the run level instead of at the subject level.
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
#' }
#' @examples
#' # see examples under TPLS_cv as you'd need a TPLS_cv object to run this function
#' @import plotly
#' @export

evalTuningParam <- function(TPLScvmdl,type=c("pearson","spearman","AUC"),X,Y,compvec,threshvec,subfold=NULL){
  if(is.null(subfold)){subfold=TPLScvmdl$testfold} # if no subdivision provided, same as CV fold

  # Perform CV prediction and performance measurement
  threshvec = sort(threshvec); compvec = sort(compvec); # sorted from low to high
  perfmat = array(NA,dim=c(length(compvec),length(threshvec),TPLScvmdl$numfold))
  for (i in 1:TPLScvmdl$numfold){
    message(paste("Fold #",i))
    testCVfold = TPLScvmdl$testfold==i
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
  out <- findBestPerf(perfmat)
  compval_best = compvec[out$row_best]; threshval_best = threshvec[out$col_best]
  compval_1se = compvec[out$row_1se]; threshval_1se = threshvec[out$col_1se]

  # output list
  out <- list("type" = type, "threshval" = threshvec, "compval" = compvec, "perfmat" = perfmat,
              "perf_best"= out$perf_best, "compval_best"=compval_best, "threshval_best" = threshval_best,
              "perf_1se"= out$perf_1se, "compval_1se"=compval_1se, "threshval_1se"=threshval_1se)
  class(out) <- "evalTuningParam"
  out
}


#' Plots the tuning surface of TPLS
#'
#' @param object : evalTuningParam object
#' @examples
#' # See examples for TPLS_cv
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

findBestPerf <- function(perfmat){
  avgperfmat = apply(perfmat,c(1,2),mean); perf_best = max(avgperfmat,na.rm = TRUE)
  bestind = which(avgperfmat==perf_best,TRUE)
  row_best = bestind[1,1]; col_best = bestind[1,2]
  standardError = stats::sd(perfmat[row_best,col_best,])/dim(perfmat)[3] # finding the standard error of the best point
  candidates = avgperfmat[,1:col_best]>(perf_best-standardError) # finding points whose metric is higher than perf_max minus 1 SE
  ind1se = which(candidates,TRUE)
  row_1se = ind1se[1,1]; col_1se=ind1se[1,2]; perf_1se = avgperfmat[row_1se,col_1se]
  return(list("perf_best"=perf_best,"row_best"=row_best,"col_best"=col_best,"perf_1se"=perf_1se,"row_1se"=row_1se,"col_1se"=col_1se))
}

util_perfmetric <- function(predmat,testY,type){
  if(type == 'auc'){
    n = length(testY); num_pos = sum(testY==1); num_neg = n - num_pos
    if (num_pos>0 && num_pos < n){
      ranks = apply(predmat,2,rank)
      return( (  colSums(ranks[testY==1,]) - num_pos*(num_pos+1)/2 ) / (num_pos*num_neg) )
    }
  }else{
    suppressWarnings({result = stats::cor(predmat, testY, method = type)})
    result[is.na(result)] = 0 # this probably happened because your prediction (or your Y) has no variation.
    return(result)
  }
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
  fig <- fig %>% add_trace(x = x[max.col(z)], y = y, z = apply(z,1,max),type = 'scatter3d', mode = 'lines+markers',
                           line = list(width = 6, color = 'rgb(50, 200, 255)'), marker = list(size = 3.5, color = 'rgb(50,50,50)'), name = 'best at component')
  fig <- fig %>% add_trace(x = ~x, y = ~y[max.col(t(z))], z = ~apply(z,2,max),type = 'scatter3d', mode = 'lines+markers',
                           line = list(width = 6, color = 'rgb(255, 200, 50)'), marker = list(size = 3.5, color = 'rgb(50,50,50)'), name = 'best at threshold')
  return(fig)
}
