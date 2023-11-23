#library(FSelector)# there is no package called ‘randomForest’
library(tidyverse)
library(xgboost)
library(CoxBoost)
library(glmnet)
library(survivalsvm)
library(rms)
library(pROC)
library(survivalROC)
library(coefplot)
library(Hmisc)
library(mlr3filters)
library(mlr3)
library(praznik)
library(SurvMetrics)
library(pec)
library(caret)
library(survXgboost)
library(randomForestSRC)
library(ranger)


coxboost.fn<- function(train.dat, test.dat){
 
  coxboost.obj <- CoxBoost(time= train.dat[,"y"], 
                           status= train.dat[,"event"], 
                           x=as.matrix(train.dat),penalty=100, stepno = 100)
  # BS
  t.interest = sort(test.dat$y[test.dat$event == 1])
  coxboost.obj.pred <- predict(coxboost.obj, 
                                 newdata=as.matrix(test.dat),
                                 newtime=test.dat[,"y"],
                                 newstatus=test.dat[,"event"],times=t.interest, type = "risk")
  
  brier_coxboost<-IBS(Surv(test.dat$y, test.dat$event), sp_matrix =  coxboost.obj.pred, t.interest)
  # C
  coxboost.obj.pred <- t(predict(coxboost.obj, 
                                 newdata=as.matrix(test.dat),
                                 newtime=test.dat[,"y"],
                                 newstatus=test.dat[,"event"]))
  test.dat <- within(test.dat, {
    ## Create a survival vector
    Surv <- Surv(y, event)
  })
  
  C.index<-rcorrcens(formula = Surv ~ I(-1 * coxboost.obj.pred), data = test.dat)
  return(list(BS=brier_coxboost,C.index=C.index))
  
}

###################################################################################
##random forest
#Overall out of bag prediction error. For classification this is the fraction of missclassified samples, for probability estimation the Brier score, 
#for regression the mean squared error and for survival one minus Harrell's C-index.


RSF.fn<- function(train.dat, test.dat, mtry=sqrt(25), nodesize=5 , ntree=1000){
  
  set.seed(123) # question, if seed is not set, c index will keep changing
  rfsrc.obj <- rfsrc(Surv(y, event)~., train.dat)
  t.interest = sort(test.dat$y[test.dat$event == 1])
  rsf.prediction = predictSurvProb(rfsrc.obj,newdata=test.dat,times=t.interest)
  brier_rsf = IBS(Surv(test.dat$y, test.dat$event), sp_matrix = rsf.prediction, t.interest)
  
  pred.rf <- predict.rfsrc(rfsrc.obj,test.dat)
  C.index<- rcorrcens(formula = Surv(y, event) ~ I(-1 * rowSums(pred.rf$chf)), data = test.dat)
  
  return(list(BS=brier_rsf,C.index=C.index))
}

##ranger(fast implementation of random forest)
ranger.fn <- function(train.dat, test.dat,mtry=sqrt(25), min.node.size=5,num.trees=1000){
  set.seed(123) # question, if seed is anot set, c index will keep changing
  ranger.obj <- ranger(Surv(y, event)~., data = train.dat, splitrule = "maxstat",importance = "permutation")
  ranger.prediction = predict(ranger.obj,data = test.dat, na.action="na.impute")
  index<- train.dat$event==1 
  Y <- sort(unique(train.dat$y[index]))
  #Y <- sort(unique(Y))
  match.ls <- match(Y, ranger.prediction$unique.death.times)
  brier_ranger = IBS(Surv(test.dat$y, test.dat$event), sp_matrix = ranger.prediction$survival[,match.ls], Y)
  
  pred <- predict(ranger.obj, test.dat)
  C.index <- rcorrcens(formula = Surv(y, event) ~ I(-1 * rowSums(pred$chf)), data = test.dat)
  
  return(list(BS=brier_ranger,C.index=C.index))
}

##################################################################################
##penalized cox
pen.cox.fn<-function(train.dat, test.dat){
  x<-as.matrix(train.dat)
  y<-Surv(train.dat$y,train.dat$event)
  cvfit <- cv.glmnet(x,y,family= "cox",alpha = 0, nfolds=5, type.measure = "C")
  nonzero.coef<-as.data.frame(extract.coef(cvfit))
  nonzero.coef<-nonzero.coef$Coefficient
  train.data<-train.dat[,names(train.dat) %in% nonzero.coef]
  fitcox = coxph(Surv(y, event)~., data = train.data, x = TRUE)
  t.interest = sort(test.dat$y[test.dat$event == 1])
  mat_cox = predictSurvProb(fitcox, test.dat, t.interest)
  brier.cox = IBS(Surv(test.dat$y, test.dat$event), sp_matrix = mat_cox, t.interest)
  
  pred <- predict(cvfit, newx = as.matrix(test.dat))
  test.dat <- within(test.dat, {
    ## Create a survival vector
    Surv <- Surv(y, event)
  })
  C.index <- rcorrcens(formula = Surv ~ I(-1 * pred), data = test.dat)# C= 0.994, Dxy=0.987 # I revised by using the same function

  return(list(BS=brier.cox,C.index=C.index))
}

#################################################################################
##xgboost                     
xgboost.fn <- function(train.dat, test.dat){
  
  set.seed(123) # again, we need to set the seed
  train.dat<-train.dat[is.finite(rowSums(train.dat)),]
  test.dat<-test.dat[is.finite(rowSums(test.dat)),]
  label <- ifelse(train.dat$event == 1, train.dat$y, -train.dat$y)
  label_val<- ifelse(test.dat$event == 1, test.dat$y, -test.dat$y)

  x_train <- as.matrix(train.dat[, !names(train.dat) %in% c("y", "event")])
  x_test<- as.matrix(test.dat[, !names(test.dat) %in% c("y", "event")])
  x_label <- label
  x_val <- xgb.DMatrix(as.matrix(test.dat[, !names(test.dat) %in% c("y", "event")]),
                       label = label_val)
  surv_xgboost_model <- xgb.train.surv(
    params = list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.05 # larger eta leads to algorithm not converging, resulting in NaN predictions
    ), data = x_train, label = x_label,
    watchlist = list(val2 = x_val),
    nrounds = 1000, early_stopping_rounds = 30
  )
  
  t.interest = sort(test.dat$y[test.dat$event == 1])
  xgboost.prediction <- predict(object = surv_xgboost_model, newdata = x_test, type = "surv", times = t.interest)
  brier.xgboost.tree<- IBS(Surv(test.dat$y, test.dat$event), sp_matrix = xgboost.prediction,t.interest)
  
  set.seed(123) # again, we need to set the seed
  train.dat<-train.dat[is.finite(rowSums(train.dat)),]
  test.dat<-as.data.frame(test.dat[is.finite(rowSums(test.dat)),])
  label <- ifelse(train.dat$event == 1, train.dat$y, -train.dat$y)
  label_val<- ifelse(test.dat$event == 1, test.dat$y, -test.dat$y)
  #
  
  y.xbg.pred <- predict(object=surv_xgboost_model, newdata=x_test, type = "risk")
  test.dat<- within(test.dat, {
    ## Create a survival vector
    Surv <- Surv(y, event)
  })
  cindex.rs<-rcorrcens(formula = Surv ~ I(-1 * y.xbg.pred), data = test.dat)
  
  return(list(BS=brier.xgboost.tree,C.index=cindex.rs))
}

####xgboost.linear
xgboost.linear.fn<-function(train.dat, test.dat){
 
  set.seed(123) # again, we need to set the seed
  train.dat<-train.dat[is.finite(rowSums(train.dat)),]
  test.dat<-as.data.frame(test.dat[is.finite(rowSums(test.dat)),])
  
  label <- ifelse(train.dat$event == 1, train.dat$y, -train.dat$y)
  label_val<- ifelse(test.dat$event == 1, test.dat$y, -test.dat$y)
  
  x_train <- as.matrix(train.dat[, !names(train.dat) %in% c("y", "event")])
  x_test<- as.matrix(test.dat[, !names(test.dat) %in% c("y", "event")])
  x_label <- label
  x_val <- xgb.DMatrix(as.matrix(test.dat[, !names(test.dat) %in% c("y", "event")]),
                       label = label_val)
  surv_xgboost_model <- xgb.train.surv(
    params = list(
      booster= "gblinear",
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.0001 # larger eta leads to algorithm not converging, resulting in NaN predictions
    ), data = x_train, label = x_label,
    watchlist = list(val2 = x_val),
    nrounds = 1000, early_stopping_rounds = 30
  )
  
  t.interest = sort(test.dat$y[test.dat$event == 1])
  xgboost.prediction <- predict(object = surv_xgboost_model, newdata = x_test, type = "surv", times = t.interest)
  brier.xgboost.linear<- IBS(Surv(test.dat$y, test.dat$event), sp_matrix = xgboost.prediction,t.interest)
  
  #
  test.dat$y.xbg.linear.pred=predict(surv_xgboost_model, x_test)
  test.dat <- within(test.dat, {
    ## Create a survival vector
    Surv <- Surv(y, event)
  })
  test.dat.nif<- test.dat[!is.infinite(test.dat$y.xbg.linear.pred),]
  cindex.rs<-rcorrcens(formula = Surv ~ I(-1 * y.xbg.linear.pred), data = test.dat.nif)
  
  return(list(BS=brier.xgboost.linear,C.index=cindex.rs))
}

########################################################################################################
#####################################Feature selection##################################################
########################################################################################################

###random forest
rfsrc.FN<- function(dt, top.num=30){
  rfsrc=rfsrc(Surv(y, event)~., dt,ntree = 1500, nsplit = 10, mtry = sqrt(1000), nodesize=15, importance = TRUE)
  vimp=as.data.frame(vimp(rfsrc)$importance)  
  vimp=rownames_to_column(vimp, var = "Gene")
  vimp.reorder=vimp[order(-vimp(rfsrc)$importance),]
  rfsrc.top.dt=vimp.reorder[1:top.num,"Gene"]
  return(rfsrc.top.dt)
}

rfsrc.md.FN<- function(dt,top.num=25){
  var.select.rs=var.select(Surv(y, event) ~ ., dt, method = "md", ntree=1000, nsplit=10, nodesize=3, splitrule="logrank")
  return(var.select.rs$topvars[1:top.num])
}

rfsrc.vh.FN<- function(dt,top.num=25){
  var.select.rs=var.select(Surv(y, event) ~ ., dt, method = "vh", ntree=1000, nodesize=3, splitrule="logrank", nsplit=10, nrep=50, K=5, nstep=1)
  return(var.select.rs$topvars[1:top.num])
}

###univariate
cox.FN<- function(dt,p.cutoff=0.05, top.num=25){
  cph.OS.DPP <- NULL
  
  fs.names<- colnames(within(dt, rm("event", "y")))
  for (i in fs.names){   ###3th column is gene
    cph.OS <- summary(coxph(Surv(y, event) ~ dt[,i],data=dt))
    
    cph.OS.DPP <- rbind(cph.OS.DPP,c(i,cph.OS$coefficients[1],cph.OS$coefficients[5]))
  }
  colnames(cph.OS.DPP) <- c("genes","Est","pvalue")
  cph.OS.DPP <- as.data.frame(cph.OS.DPP)
  cph.OS.DPP$Est <- as.numeric(cph.OS.DPP$Est)
  cph.OS.DPP$pvalue<-as.numeric(cph.OS.DPP$pvalue)
  sig.gene=cph.OS.DPP %>% filter(pvalue <=0.05) %>% arrange(desc(Est))%>% dplyr::select(genes) 
  cox.sig.gene=sig.gene[1:top.num,]
  return(cox.sig.gene)
}

###lasso
lasso.FN<- function(dt,top.num=25,alpha = 0, nfolds=5){
  
  lasso.x<-as.matrix(dt)
  lasso.y<-Surv(dt[,"y"],dt[,"event"])
  lasso.cvfit <- cv.glmnet(lasso.x,lasso.y,family= "cox")
  lasso.coef=extract.coef(lasso.cvfit) 
  lasso.coef=lasso.coef[-c(1:2),]
  lasso.coef.reorder=lasso.coef%>% arrange(desc(Value))
  lasso.coef.reorder.top.dt=lasso.coef.reorder[1:top.num,"Coefficient"] 
  return(lasso.coef.reorder.top.dt)
}

lasso1 <- function(out,lam,numCol){
  lasso.pred  <- predict(out, type="coefficients",s=lam)[1:numCol,]
  return(lasso.pred)
}

lasso.FN.1<-function(dt){
  x1 <- model.matrix(Surv(y,event)~.,dt)[,-1]
  y1 <- Surv(dt$y,dt$event)
  lambdas <- 10^seq(10, -10, by = -.1)
  
  #j=0
  lassoCheck1 = 0
  while(lassoCheck1 == 0){
    cvfit1 <- withCallingHandlers(cv.glmnet(x1,y1,alpha=1,family = "cox", nfolds=5, lambda= lambdas,type.measure = "C"), 
                                  warning=function(w){invokeRestart("muffleWarning")})
    bestlam1=cvfit1$lambda.1se
    
    # muffle warnings, get  coefficients 
    call.coef1<- withCallingHandlers(lasso1(cvfit1, bestlam1, ncol(x1)), warning=function(w){invokeRestart("muffleWarning")})
    
    Lasso.top.dt <- names(call.coef1[call.coef1 != 0])
    lassoCheck1 <- length(Lasso.top.dt)
    
    return(Lasso.top.dt)
  }
  
  lasso.FN.1()
  
}


mim.FN<- function(dt,thread=1){
  X<- within(dt, rm("event", "y"))
  mim.results<- MIM(X=X, Y=dt$event,  25)
  return(names(mim.results$selection))
}

mrmr.FN<- function(dt,thread=1){
  X<- within(dt, rm("event", "y"))
  mrmr.results<- MRMR(X=X, Y=dt$event,  25)
  return(names(mrmr.results$selection))
}

Gen.ranking <- function(features.top,feature.name,order=TRUE)
{
  rank <- c()
  if (order)
  {
    for (name in feature.name)
    { 
      if (is.element(name,features.top))
        rank <- append(rank,which(features.top==name))
      else
        rank <- append(rank,length(feature.name))
    }
  }
  else
  {
    for (name in feature.name)
    { 
      if (is.element(name,features.top))
        rank <- append(rank,length(features.top)/2) #ave ranking
      else
        rank <- append(rank,length(feature.name))
    }
  }
  
  return(rank)
}


spglasso <- function(X,y,n_perm,group,penalty,cut.off,block.info)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X)
  Spi <- NULL
  
  group.pseudo <- group+nrow(block.info)
  pp <- rep(nrow(block.info)+p, times=nrow(block.info))
  penalty.pseudo <- pp*sqrt(block.info$size)
  for(k in 1:n_perm){
    X_ko1 <- X[c(sample(nrow(X))), ]
    colnames(X_ko1) <- paste(colnames(X),"_pc",sep="")
    ii <- 1
    m1 <- grpsurv(X=cbind(X, X_ko1), y=y, group=c(group,group.pseudo), penalty="grLasso",  alpha=1,group.multiplier=c(penalty,penalty.pseudo) ,nlambda=100)
    glasso <- as.matrix(m1$beta)
    while (max(abs(glasso[(p+1):(2*p),ii]))==0 & ii<dim(glasso)[2]){
      ii <- ii+1 
    }
    selected_lasso <-  which(abs(glasso[,ii-1])>0)
    Spi <- c(Spi,selected_lasso)
  }
  Spi <- names(Spi)
  freq <- table(Spi)/n_perm
  out <- names(freq)[as.numeric(freq)>=cut.off]
  #out
  return(out)
}

FS.ensemble <- function(dt,rank.sum,rho.t=0.75,K=50,tau=0.5)
{
  X <- dt[,-c(1,2)]
  #Y <- as.factor(dt$y)
  y <- dt[,c(1,2)]
  
  ##removing constant columns in X
  sd.X <- apply(X,2,sd)
  var0.idx <- which(sd.X==0)
  if (length(var0.idx)!=0)
  {
    X <- X[,-var0.idx]
  }
  
  ##generate correlation blocks####
  cor.x <- abs(cor(X,method="spearman"))
  hc <- hclust(as.dist(1-cor.x)) 
  names.fix <- colnames(cor.x)
  k <- 0
  min.cor <- 0
  while (min.cor<rho.t)
  {
    k <- k+1
    memb <- cutree(hc, k)
    mincor.list <- c()
    cor.block <- list()
    for(i in 1:k)
    {
      subset <- which(memb==i)
      abs.cor.subset <- abs(cor.x[subset,subset])
      mincor.list <- append(mincor.list,min(abs.cor.subset))
      cor.block[[i]] <- names(subset)
    }
    min.cor <- min(mincor.list)
  }
  #generate block-info
  rank.min <- apply(rank.sum,2,median)
  block.info <- data.frame(cbind(1:k,as.numeric(lapply(cor.block,function(x) min(rank.min[x])))))
  names(block.info) <- c("bid","rank.min")
  #penalize using the length directly
  block.info$size <- as.numeric(lapply(cor.block, length))
  
  #penalize using the length of large rank variables in the block
  block.info$psize <- as.numeric(lapply(cor.block, function(x) ifelse(max(rank.min[x])<=50, 1, sum(rank.min[x]>= mean(rank.min[x])))))
  
  block.info$rank.final <- rank(block.info$rank.min,ties.method="max")
  block.info <- block.info[order(block.info$rank.min),]
  row.names(block.info) <- NULL
  cor.block.order <- cor.block[block.info$bid]
  block.all.order <- unlist(cor.block.order)

  X <- X[block.all.order]
  group <- unlist(mapply(rep, 1:nrow(block.info), block.info$size))
  penalty <- block.info$rank.final*sqrt(block.info$size)
  out <- spglasso(X,y,K, group, penalty,tau,block.info)
  out
  return(out)
  
}

