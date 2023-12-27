#!/usr/bin/env Rscript
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-1200

rm(list=ls())

library(SFSI)
library(BGLR)

load('YBC_GWASAssistedGP.RData')

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
JOBS <- expand.grid(rep=1:100, trait=1:12)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

phenotype <- YBC_GWASAssistedGP$pheno
rownames(phenotype) <- phenotype$taxa

traits <- c("Yield_MI18","Yield_MI19","Yield_NE18", "Yield_NE19",
            "Fe_MI18", "Fe_MI19","Fe_NE18", "Fe_NE19", 
            "Zn_MI18","Zn_MI19","Zn_NE18","Zn_NE19")

phenotype <- phenotype[,traits]  


Y = phenotype
M = YBC_GWASAssistedGP$geno
M = scale(M)/sqrt(ncol(M))
G = tcrossprod(M)
D<-as.matrix(dist(M,method="euclidean"))^2 #euclidian distance
D<-D/mean(D)

index <- grep("YBC",rownames(Y))
G <- G[index,index]
D <- D[index,index]
Y <- Y[index,]

trait_name <- colnames(Y)[trait]
y <- Y[,trait]

index <- which(!is.na(y))# remove missing values
G <- G[index,index]
D <- D[index,index]
y <- y[index]

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

### G-Blup
GBLUP = cor(fm0$u[tst],y[tst])


# Obtain an 'optimal' lambda by repeating the CV several times
fm1 = SSI.CV(y,K=G,trn=trn,
              nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn, tst=tst,lambda=lambda)
SSI = summary(fm2)$optCOR

## KA
h<- c(.02,1,5) # bandwidth kernels
KList<-list() 

for(i in 1:length(h)){ # Avering kernel
  KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
}

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_",trait_name,"_rep_",rep))

KA <- cor(y[tst],fmKA$yHat[tst],
         use = 'pairwise.complete.obs')


## ==============

results = c(GBLUP,SSI, KA)

names(results) = c('GBLUP','SSI',
                   'SSI_MSE', 'SSI_df','SSI_lambda', "KA")


save(results,
     file=paste0("results_",trait_name,"_rep_",rep,".RData"))