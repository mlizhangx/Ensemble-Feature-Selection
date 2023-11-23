library(mlr)
library(tidyverse)
library(openxlsx)
library(randomForestSRC)
library(ranger)
library(mRMRe)
library(xgboost)
library(CoxBoost)
library(glmnet)
library(survivalsvm)
library(rms)
library(survivalROC)
library(pROC)
library(coefplot)
library(Hmisc)
library(caret)
library(praznik)
library(stringr)
library(reshape2)
library(grpreg)
library(dcurves)
library(tictoc)
#library(mlr3filters)

############################Colorectal##############################
load("./Colorectal.RData")
source("./source_function_PL.R")
# set.seed(123)
cph.OS.fs <- NULL
for (i in 25:19971){
  cph.OS <- summary(coxph(Surv(OS.time.month,OS) ~ dt_merge[,i],data=dt_merge))
  cph.OS.fs <- rbind(cph.OS.fs,c(colnames(dt_merge)[i],cph.OS$coefficients[1],cph.OS$coefficients[5]))
}

colnames(cph.OS.fs) <- c("genes","Est","pvalue")
cph.OS.fs <- as.data.frame(cph.OS.fs)
cph.OS.fs$pvalue<-as.numeric(cph.OS.fs$pvalue)
cph.OS.fs$Est <- as.numeric(cph.OS.fs$Est)


sig.gene <- cph.OS.fs[cph.OS.fs$pvalue < 0.05 & !is.na(cph.OS.fs$pvalue),]# Li
## checking # among 19947 genes,521 have no pvalue, 2303 significant genes
dt.merge.subset <- dt_merge[,c("OS.time.month","OS",sig.gene$genes)]
colnames(dt.merge.subset) <- make.names(colnames(dt.merge.subset),unique = T)

dt.merge.subset.norm=cbind(dt.merge.subset[,1:2], as.data.frame(scale(dt.merge.subset[,3:ncol(dt.merge.subset)])))
colnames(dt.merge.subset.norm)[1:2] <- c("y", "event")
dt.merge.subset.norm.nna<-dt.merge.subset.norm[!is.na(dt.merge.subset.norm$event),]
dt.merge.subset.norm.nna$event= ifelse(dt.merge.subset.norm.nna$event==1,0,1)
dt.merge.subset.norm.nna= dt.merge.subset.norm.nna[dt.merge.subset.norm.nna$y!=0,]

####################################################################################################################
################################################Real data############################################################
rf.top.dt <- list()
cox.top.dt <- list()
lasso.top.dt <- list()
rf.md.top.dt<- list()
rf.vh.top.dt<- list()
mim.top.dt<- list()
mrmr.top.dt<- list()
ensemble.top.dt<- list()
ensemble2.top.dt<-list()

# save classification results
ClassResults <-data.frame(matrix(nrow=0,ncol=0))

# create dataframe to hold evaluation of FS method results
FSResults <- as.data.frame(matrix(NA,nrow=0,ncol=length(name_evaluate)))
colnames(FSResults) <- name_evaluate

# tic("Brier scores heatmap")
n_sim = 5

i =1

#n_sim <- 2
while(i <= n_sim) {
  
  dt <- dt.merge.subset.norm.nna
  # remove the row 93 due to outlier, drop NA any record that has NA
  dt <- dt[-93,] %>% drop_na() %>% dplyr::select(-c("LELP1.149018"))
  flds <- createFolds(dt$y, k = 5)

  # save feature names for evaluation criteria
  feature.name <-  names(dt[,-c(1,2)])
  n_gene <- ncol(dt)-2
  rank.sum <- as.data.frame(matrix(NA,nrow=4,ncol=n_gene))
  rownames(rank.sum) <- c("RF", "Uni.Cox",  "LASSO","MIM")
  colnames(rank.sum) <- feature.name
  
  ##RF#######################
  rf.top.dt[[i]] <- rfsrc.FN(dt[-flds[[i]],])  
  rank.rf <- Gen.ranking(rf.top.dt[[i]],feature.name)
  rank.sum["RF",] <- rank.rf
  
  #Uni.Cox################
  cox.top.dt[[i]] <- cox.FN(dt[-flds[[i]],])
  rank.cox <- Gen.ranking(cox.top.dt[[i]],feature.name)
  rank.sum["Uni.Cox",] <- rank.cox
  
  ######LASSO
  lasso.top.dt[[i]] <-lasso.FN.1(dt[-flds[[i]],])
  rank.lasso <- Gen.ranking(lasso.top.dt[[i]],feature.name)
  rank.sum["LASSO",] <- rank.lasso
  
  #RF minimal depth
  rf.md.top.dt[[i]] <- rfsrc.md.FN(dt[-flds[[i]],])
  
  #RF variable hunting
  rf.vh.top.dt[[i]]<- rfsrc.vh.FN(dt[-flds[[i]],])
  #MIM
  mim.top.dt[[i]] <- mim.FN(dt[-flds[[i]],])
  rank.mim <- Gen.ranking(mim.top.dt[[i]],feature.name)
  rank.sum["MIM",] <- rank.mim
  
  #MRMR
  mrmr.top.dt[[i]] <- mrmr.FN(dt[-flds[[i]],])
  #Ensemble###############
  ensemble.top.dt[[i]]<- FS.ensemble(dt[-flds[[i]],], rank.sum[1:3,])
  ensemble2.top.dt[[i]]<- FS.ensemble(dt[-flds[[i]],], rank.sum)
  #######################  classification/Prediction  ########################
  
  #retrieve first sim features
  method.set <- c(1:9)
  for(v in method.set){
    # v=1, v=2, v=3, v=4, v=5
    # grab 1 of the FS methods features, perform all classifiers
    # go through each FS list of features'
    
    if(v ==1){feature=rf.top.dt; fsname="RF"}
    if(v ==2){feature=cox.top.dt; fsname="Uni.cox"}
    if(v ==3){feature=lasso.top.dt; fsname="LASSO"}
    if(v ==4){feature=rf.md.top.dt; fsname="RF.MD"}
    if(v ==5){feature=rf.vh.top.dt; fsname="RF.VH"}
    if(v ==6){feature=mim.top.dt; fsname="MIM"}
    if(v ==7){feature=mrmr.top.dt; fsname="MRMR"}
    if(v ==8){feature=ensemble.top.dt; fsname="Ensemble1"}
    if(v ==9){feature=ensemble2.top.dt; fsname="Ensemble2"}
    
    fs <- as.character(feature[[i]])
    fs <- na.omit(fs)
    # obtain dataset with selected features, if theres only 1 feature or more than 1 feature.
    if(length(fs) == 1){
      Modeldata <- data.frame(dt[-flds[[i]],c(fs)],y=dt[-flds[[i]],"y"], event= dt[-flds[[i]], "event"])
      colnames(Modeldata) <-c(fs,"y", "event")
      
      Testdata <- data.frame(dt[flds[[i]],c(fs)],y=dt[flds[[i]],"y"], event= dt[flds[[i]], "event"])
      colnames(Testdata) <-c(fs,"y", "event")
    } else if (length(fs) == 0) {
      
      next}
      
      else {
      Modeldata <- data.frame(dt[-flds[[i]],c(fs)],y=dt[-flds[[i]],"y"], event= dt[-flds[[i]], "event"])
      Testdata <- data.frame(dt[flds[[i]],c(fs)], y=dt[flds[[i]],"y"], event= dt[flds[[i]], "event"])
      }
    
    
    train.dat <- Modeldata
    test.dat <- Testdata
    
    # PROBLEMS: xgboost, xgboost, Lasso,coxboost, ranger, rf
    xgb_results <- c(xgboost.fn(train.dat, test.dat))
    xgblinear_results <- c(xgboost.linear.fn(train.dat, test.dat))
    lasso_results<- c(pen.cox.fn(train.dat, test.dat))
    rf_results <- c(RSF.fn(train.dat, test.dat))
    ranger_results <- c(ranger.fn(train.dat, test.dat))
    coxboost_results <- c(coxboost.fn(train.dat, test.dat))
    
    ClassResults[i, paste(fsname,"xgb","Brier",sep="_")] <- xgb_results[[1]][1] 
    ClassResults[i, paste(fsname,"xgblinear","Brier",sep="_")] <- xgblinear_results[[1]][1] 
    ClassResults[i,paste(fsname,"lasso","Brier",sep="_")] <- lasso_results[[1]][1]
    ClassResults[i, paste(fsname,"rf","Brier",sep="_")] <- rf_results[[1]][1]
    ClassResults[i,paste(fsname,"ranger","Brier",sep="_")] <- ranger_results[[1]][1] 
    ClassResults[i,paste(fsname,"coxboost","Brier",sep="_")] <- coxboost_results[[1]][1] 
    
    ClassResults[i, paste(fsname,"xgb","C.index",sep="_")] <- xgb_results[[2]][1] 
    ClassResults[i, paste(fsname,"xgblinear","C.index",sep="_")] <- xgblinear_results[[2]][1] 
    ClassResults[i,paste(fsname,"lasso","C.index",sep="_")] <- lasso_results[[2]][1]
    ClassResults[i, paste(fsname,"rf","C.index",sep="_")] <- rf_results[[2]][1]
    ClassResults[i,paste(fsname,"ranger","C.index",sep="_")] <- ranger_results[[2]][1] 
    ClassResults[i,paste(fsname,"coxboost","C.index",sep="_")] <- coxboost_results[[2]][1] 
    
  }  
  
  
  i = i +1
  
}

save.image("./Results11212023.RData")

