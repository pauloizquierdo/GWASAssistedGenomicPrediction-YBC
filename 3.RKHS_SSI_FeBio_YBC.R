#!/usr/bin/env Rscript
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-200

rm(list=ls())

library(SFSI)
library(BGLR)

load('YBC_GWASAssistedGP.RData')

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
JOBS <- expand.grid(rep=1:100, trait=1:2)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

rownames(YBC_GWASAssistedGP$pheno) <- YBC_GWASAssistedGP$pheno$taxa

FeBio <- c("FeBio_MI18", "FeBio_MI19")
index <- which(colnames(YBC_GWASAssistedGP$pheno) %in% FeBio)

FeBio_traits <- YBC_GWASAssistedGP$pheno[,index]

Y <- FeBio_traits
M <- YBC_GWASAssistedGP$geno
M <- scale(M)/sqrt(ncol(M))
G <- tcrossprod(M)
D <- as.matrix(dist(M,method="euclidean"))^2 #euclidian distance
D <- D/mean(D)

index <- grep('YBC', rownames(Y))
G <- G[index,index]
D <- D[index,index]
Y <- Y[index,]
M <- M[index,]

trait_name <- colnames(Y)[trait]
y <- Y[,trait]

index <- which(!is.na(y))
G <- G[index,index]
D <- D[index,index]
y <- y[index]
M <- M[index,]
Y <- Y[index,]

n <- length(y)

set.seed(rep)
tst <- sample(1:n, 0.3*n) 
trn <- seq_along(1:n)[-tst]

# Calculate variance components ratio using training data

yNA = y
yNA[tst] = NA
fm0 = fitBLUP(yNA,K=G)
theta = fm0$varE/fm0$varU
h2 = fm0$varU/(fm0$varU + fm0$varE)
b = fm0$b # intercept

# Obtain an 'optimal' lambda by repeating the CV several times
fm1 = SSI.CV(y,K=G,trn=trn,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn, tst=tst,lambda=lambda)
SSI = summary(fm2)$optCOR

## KA
h <- c(.02,1,5) # bandwidth kernels
KList<-list() 

for(i in 1:length(h)){ # Avering kernel
  KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
}

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_",trait_name,"_rep_",rep))

KA = cor(y[tst],fmKA$yHat[tst],
         use = 'pairwise.complete.obs')


############ Models with QTLs
# Obtain an 'optimal' lambda by repeating the CV several times

index <- which(colnames(YBC_GWASAssistedGP$pheno) %in% FeBio)
FeBio_traits <- YBC_GWASAssistedGP$pheno[,index]

Y = FeBio_traits
M = YBC_GWASAssistedGP$geno

indexQTL <- which(colnames(M) %in% 
                    c('S07_29169848')) ## QTl for yield

M_new <- M[,-indexQTL]
QTL <- as.matrix(M[,indexQTL])

M_new=scale(M_new)/sqrt(ncol(M_new))
G_qtl = tcrossprod(M_new)
D_qtl<-as.matrix(dist(M_new,method="euclidean"))^2 #euclidian distance
D_qtl<-D_qtl/mean(D_qtl)

index <- grep('YBC', rownames(Y))
G_qtl <- G_qtl[index,index]
D_qtl <- D_qtl[index,index]
M_new <- M_new[index,]
Y <- Y[index,]
QTL <- QTL[index,]
y <- Y[,trait]

index <- which(!is.na(y)) # remove missing values
G_qtl <- G_qtl[index,index]
D_qtl <- D_qtl[index,index]
QTL <- QTL[index]
y <- y[index]

#### SSI + QTL

fm4 = SSI.CV(y, b=QTL, K=G_qtl,trn=trn,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm4)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm5 = SSI(y,b=QTL,K=G_qtl,trn=trn, tst=tst,lambda=lambda)
SSI_QTL = summary(fm5)$optCOR


### KA - QTL
## KA
h <- c(.02,1,5) # bandwidth kernels
KList <-list()
KList <- list(list(X=QTL,model="FIXED"))

for(i in 1:length(h)){
  KList[[i+1]]<- list(K=exp(-h[i]*D_qtl),model='RKHS')
}

fmKA_qtl<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_qtl_",trait_name,"_rep_",rep))

KA_qtl = cor(y[tst],fmKA_qtl$yHat[tst],
         use = 'pairwise.complete.obs')


results = c(SSI,SSI_QTL,KA,KA_qtl)

names(results) = c('SSI','SSI_MSE', 'SSI_df','SSI_lambda',
                   'SSI_qtl','SSI_MSE_qtl', 'SSI_df_qtl','SSI_lambda_qtl',
                   "KA", "KA_qtl")

save(results,
     file=paste0("results_",trait_name,"_rep_",rep,".RData"))





  
