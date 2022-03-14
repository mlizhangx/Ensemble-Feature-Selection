#Ensemble feature selection
#loading packages
library(caret) #svmrfe
library(vita) #vita
library(Boruta) #boruta
library(FSelector) #information gain
library(gglasso) #group lasso

#######################
####define functions####
########################


####Generate rankings####
#input1: selected features
#input2: the complete list of feature names
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
        rank <- append(rank, length(features.top)/2) #max ranking
      else
        rank <- append(rank,length(feature.name))
    }
  }

  return(rank)
}

######Feature selection ensemble##################
#input1: data, first column as Y, the rest as X
#input2: rank summary

FS.ensemble <- function(dt,rank.sum,rho.t=0.75,K=50,tau=0.5)
{
  X <- dt[,-1]
  Y <- as.factor(dt$y)
  
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
  rank.min <- apply(rank.sum,2,min)
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
  #block.info[1:30,]
  #
  y <- 2*as.numeric(Y)-3
  X <- X[block.all.order]
  group <- unlist(mapply(rep, 1:nrow(block.info), block.info$size))
  penalty <- block.info$rank.final*sqrt(block.info$psize)
  out <- pglasso(X,y,K, group, penalty,tau,block.info)
  #out
  return(out)
  
}

#y is +1/-1

pglasso <- function(X,y,n_perm,group,penalty,cut.off,block.info)
{
  #n_perm <- 10
  X <- as.matrix(X)
  #y <- 2*as.numeric(dt$y)-3
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X)
  Spi <- NULL
  
  group.pseudo <- group+nrow(block.info)
  pp <- rep(nrow(block.info)+p, times=nrow(block.info))
  penalty.pseudo <- pp*sqrt(block.info$size)
  for(k in 1:n_perm){
    X_ko1 <- X[c(sample(nrow(X))), ]
    ii <- 1
    m1 <- gglasso(x=cbind(X, X_ko1),y=y,group=c(group,group.pseudo),loss="logit",pf=c(penalty,penalty.pseudo))
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
svm.RFE.FN <- function(dt){
  dt$y <- as.factor(dt$y)
  p <- ncol(dt)-1
  rfe.ctrl <- rfeControl(functions=caretFuncs, 
                         method = "cv",
                         number = 5,
                         verbose = F)
  subsets <- seq(floor(0.01*p),floor(0.1*p),by=3)

  svm.rfe <- rfe(y~.,dt,
                 sizes = subsets,
                 rfeControl = rfe.ctrl, metric = "Accuracy", method="svmRadial", verbose = FALSE)
  return(predictors(svm.rfe))
}


vita.FN <- function(dt){
  # https://cran.r-project.org/web/packages/vita/vita.pdf
  # cv permutation variable importance 
  dt$y <- as.factor(dt$y)
  vita_results<- CVPVI(X=dt[,-1], y=dt[,"y"],k=5)
  # novel test approach
  vita_nta <- NTA(vita_results$cv_varim)
  #vita_p_summary <- summary(vita_p,pless = 0.1)
  vita_p<- vita_nta$pvalue
  # order p-values
  vita.ord<- vita_p[order(vita_p),]
  # find p-values==0
  vita.sig <-vita.ord[vita.ord ==0]
  return(names(vita.sig))
}





boruta.FN <- function(dt){
  #Boruta (wrapper with random forest)
  dt$y <- as.factor(dt$y)
  Boruta.Ozone <- Boruta(y ~ ., data = dt)
  # obtain features
  Boruta.out <- names(Boruta.Ozone$finalDecision)[which(Boruta.Ozone$finalDecision!="Rejected")]
  return(Boruta.out)
}











####################
## main function ###
####################

#loading test data
#setwd("~/Dropbox/research/TCR/Comparative study/package")
dt <- read.csv("data.csv",header=T)


feature.name <- names(dt[,-1])
rank.sum <- as.data.frame(matrix(NA,nrow=4,ncol=ncol(dt)-1))
rownames(rank.sum) <- c("SVM_RFE", "Boruta",  "Info.Gain","Vita")
colnames(rank.sum) <- feature.name


##svmRFE#######################
svmRFE.top.dt <- svm.RFE.FN(dt)
rank.sum["SVM_RFE",] <- Gen.ranking(svmRFE.top.dt,feature.name)


####Vita#############
Vita.top.dt <- vita.FN(dt)
rank.sum["Vita",] <- Gen.ranking(Vita.top.dt,feature.name,order=FALSE)

######Boruta
Boruta.top.dt <-boruta.FN(dt)
rank.sum["Boruta",] <- Gen.ranking(Boruta.top.dt,feature.name,order=FALSE)


#Information gain
Info.top.dt <- info.FN(dt)
rank.sum["Info.Gain",] <- Gen.ranking(Info.top.dt,feature.name)



#ensemble
ensemble.top.dt <- FS.ensemble(dt,rank.sum)
ensemble.top.dt

