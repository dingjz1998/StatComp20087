#' @title A permutation function using R
#' @description Change the label of a certain variable to calculate its importance
#' @param dt The initial data frame
#' @param i  The label of variable to be permuted
#' @return A permuted data frame
#' @useDynLib StatComp20087
#' @examples 
#' \dontrun{
#' dt=iris
#' i=3
#' permudata(dt,i)
#' } 
#' @export
permudata=function(dt,i){
  n=nrow(dt)
  index=sample(n, size = n, replace = FALSE)
  temp=dt[index,i]
  dt[i]=temp
  return(dt)
} 

#' @importFrom stats na.omit
#' @importFrom Rcpp evalCpp
#' @import randomForest
#' @import MASS
#' @import DAAG
#' @import bootstrap
#' @import boot
#' @import Ball
#' @import RANN
#' @import energy 
#' @import microbenchmark
#' @importFrom graphics hist lines
#' @importFrom stats cov qchisq rnorm runif var


#' @title A function using R to calculate the permutation importance 
#' @description Calculate the permutation importance of variables with random forest algorithm. This function can be used to 
#' calculate the PI for classification and regression problems. Pay attention that in the input data frame, the dependent variable 'Y' should be placed in the second column.
#' @param dt The initial data frame
#' @param ntree  The number of decision trees in the forest
#' @return The permutation importance vector of each variable\code{PI}
#' @useDynLib StatComp20087
#' @examples 
#' \dontrun{
#' dt=longley
#' ntree=500
#' rfpi(dt,ntree)
#' } 
#' @export
rfpi=function(dt,ntree){
  dt=na.omit(dt)
  n=nrow(dt)
  p=ncol(dt)
  temp=numeric(p-1)
  set.seed(65432)
  rf_ntree0=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
  mean0=mean(rf_ntree0$mse)
  for (i in 2:p){
    datai=permudata(dt,i)  
    sol=randomForest(datai[,1]~., data=datai, ntree=ntree,important=TRUE,proximity=TRUE)
    temp[i-1]=mean(sol$mse)
  }
  PI=temp-mean0
  return(PI)
}


#' @title Non recursive feature elimination (NRFE)
#' @description A non recursive function using random forest to select variables and calculate the error.This function can be used to 
#' calculate the PI for classification and regression problems. Pay attention that in the input data frame, the dependent variable 'Y' should be placed in the second column. 
#' @param dt The initial data frame
#' @param ntree  The number of decision trees in the forest
#' @return A vector \code{err}
#' @useDynLib StatComp20087
#' @examples 
#' \dontrun{
#' dt=longley
#' ntree=800
#' NRFE(dt,ntree)
#' } 
#' @export
NRFE=function(dt,ntree){
  PI=rfpi(dt,ntree) 
  n=nrow(dt)
  p=ncol(dt)
  dele=order(PI)
  err=numeric(p-1)
  set.seed(32145)
  mse0=mean(randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)$mse)
  err[1]=mse0
  for (i in 2:(p-1)){
    delev=dele[1]+1
    dt=dt[,-delev]
    sol=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
    err[i]=mean(sol$mse)
    PI=PI[-(delev-1)]
    dele=order(PI)
  }
  return(err)
}

#' @title Recursive feature elimination (RFE)
#' @description A recursive function using random forest to select variables and calculate the error. This function can be used to 
#' calculate the PI for classification and regression problems. Pay attention that in the input data frame, the dependent variable 'Y' should be placed in the second column. 
#' @param dt The initial data frame
#' @param ntree  The number of decision trees in the forest
#' @return A vector \code{err}
#' @useDynLib StatComp20087
#' @examples 
#' \dontrun{
#' dt=longley
#' ntree=500
#' RFE(dt,ntree)
#' } 
#' @export
RFE=function(dt,ntree){
  n=nrow(dt)
  p=ncol(dt)
  err=numeric(p-1)
  set.seed(32145)
  mse0=mean(randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)$mse)
  err[1]=mse0
  for (i in 2:(p-1)){
    PI=rfpi(dt,ntree) 
    dele=order(PI)
    delev=dele[1]+1
    dt=dt[-delev]
    sol=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
    err[i]=mean(sol$mse)
  }
  return(err)
}


