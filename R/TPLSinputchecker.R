#' TPLSinputchecker
#'
#' Checks for input errors on TPLS, TPLS_cv, and evalTuningParam
#' @noRd

TPLSinputchecker <- function(input,name,type=NULL,maxval=NULL,minval=NULL,variation=0,integercheck=0){
  if(any(!is.numeric(input))){stop(paste(name,"should be numeric", sep=" "))} # numeric check
  if(any(is.nan(input))){stop(paste("NaN found in",name,sep=" "))} # nan check
  if(any(is.infinite(input))){stop(paste("Non finite value found in",name,sep=" "))} # inf check

  if(!is.null(type)){
    n = nrow(input); v = ncol(input)
    switch(type,
           mat={
              if(v<3){stop(paste(name,"should have at least 3 columns", sep=" "))}
              if(n<3){stop(paste(name,"should have at least 3 observations", sep=" "))}
           },
           vec={
              if(!is.vector(input)){stop(paste(name,"should be a vector", sep=" "))}
           },
           colvec={
              if(v!=1){stop(paste(name,"should be a column vector", sep=" "))}
           },
           scalar={
              if(!is.atomic(input) | length(input) != 1){stop(paste(name,"should be a scalar", sep=" "))}
           },
           {
              stop("Unexpected input type checking requested")
           }
    )
  }

  if(!is.null(maxval)){
    if(any(input>maxval)){stop(paste(name,"should be less than or equal to",maxval, sep=" "))}
  }

  if(!is.null(minval)){
    if(any(input<minval)){stop(paste(name,"should be greater than or equal to",maxval, sep=" "))}
  }

  if(variation==1){
    if(stats::sd(input)==0){stop(paste("there is no variation in",name, sep=" "))}
  }

  if(integercheck==1){
    test<- all.equal(input,as.integer(input))
    if(test!=TRUE){stop(paste(name,"should be integer",sep=" "))}
  }
}
