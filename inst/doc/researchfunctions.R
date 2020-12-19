## ----eval=FALSE---------------------------------------------------------------
#  permudata=function(dt,i){
#    n=nrow(dt)
#    index=sample(n, size = n, replace = FALSE)
#    temp=dt[index,i]
#    dt[i]=temp
#    return(dt)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(randomForest)
#  rfpi=function(dt,ntree){
#  dt=na.omit(dt)
#  n=nrow(dt)
#  p=ncol(dt)
#  temp=numeric(p-1)
#  set.seed(65432)
#  rf_ntree0=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
#  mean0=mean(rf_ntree0$mse)
#  for (i in 2:p){
#  datai=permudata(dt,i)
#  sol=randomForest(datai[,1]~., data=datai, ntree=ntree,important=TRUE,proximity=TRUE)
#  temp[i-1]=mean(sol$mse)
#  }
#  PI=temp-mean0
#  return(PI)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NRFE=function(dt,ntree){
#  PI=rfpi(dt,ntree)
#  n=nrow(dt)
#  p=ncol(dt)
#  dele=order(PI)
#  err=numeric(p-1)
#  set.seed(32145)
#  mse0=mean(randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)$mse)
#  err[1]=mse0
#  for (i in 2:(p-1)){
#    delev=dele[1]+1
#    dt=dt[,-delev]
#    sol=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
#    err[i]=mean(sol$mse)
#    PI=PI[-(delev-1)]
#    dele=order(PI)
#  }
#  return(err)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  RFE=function(dt,ntree){
#    n=nrow(dt)
#    p=ncol(dt)
#    err=numeric(p-1)
#    set.seed(32145)
#    mse0=mean(randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)$mse)
#    err[1]=mse0
#    for (i in 2:(p-1)){
#      PI=rfpi(dt,ntree)
#      dele=order(PI)
#      delev=dele[1]+1
#      dt=dt[-delev]
#      sol=randomForest(dt[,1]~., data=dt, ntree=ntree,important=TRUE,proximity=TRUE)
#      err[i]=mean(sol$mse)
#    }
#    return(err)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20087)
dt=longley
ntree=1000
NRFE(dt,ntree)
RFE(dt,ntree)

