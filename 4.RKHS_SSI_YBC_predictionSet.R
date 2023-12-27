#!/usr/bin/env Rscript
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --constraint=intel16|intel18
#SBATCH --array=1-300

rm(list=ls())

library(SFSI)
library(BGLR)

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
JOBS <- expand.grid(rep=1:100, trait=1:3)

trait <- as.vector(JOBS[job,"trait"])
rep <- as.vector(JOBS[job,"rep"])

load('YBC_GWASAssistedGP.RData')
phenotype <- YBC_GWASAssistedGP$pheno
rownames(phenotype) <- phenotype$taxa

traits <- c('Yield_MI19','Fe_MI19','Zn_MI19') # traits to be predicted
phenotype <- phenotype[,traits]

M = YBC_GWASAssistedGP$geno
M=scale(M)/sqrt(ncol(M))
G = tcrossprod(M)
D<-as.matrix(dist(M,method="euclidean"))^2 #euclidian distance
D<-D/mean(D)

index <- which(!is.na(phenotype[,trait]))
Y = phenotype[index,]
G = G[index,index]
D = D[index,index]

trait_name <- colnames(Y)[trait]
y <- Y[,trait]
n <- length(y)

set.seed(rep)

# Define training and testing sets
index <- grep('YBC', rownames(Y))
trn0 <- index
tst0 <- seq_along(1:n)[-trn0]

Y1608 = grep('Y1608', rownames(Y))
Y1608_10 = sample(Y1608, ceiling(length(Y1608)*0.1))
Y1608_20 = sample(Y1608, ceiling(length(Y1608)*0.2))
Y1608_30 = sample(Y1608, ceiling(length(Y1608)*0.3))

Y1609 = grep('Y1609', rownames(Y))
Y1609_10 = sample(Y1609, ceiling(length(Y1609)*0.1))
Y1609_20 = sample(Y1609, ceiling(length(Y1609)*0.2))
Y1609_30 = sample(Y1609, ceiling(length(Y1609)*0.3))

Y1612 = grep('Y1612', rownames(Y))
Y1612_10 = sample(Y1612, ceiling(length(Y1612)*0.1))
Y1612_20 = sample(Y1612, ceiling(length(Y1612)*0.2))
Y1612_30 = sample(Y1612, ceiling(length(Y1612)*0.3))

Y1701 = grep('Y1701', rownames(Y))
Y1701_10 = sample(Y1701, ceiling(length(Y1701)*0.1))
Y1701_20 = sample(Y1701, ceiling(length(Y1701)*0.2))
Y1701_30 = sample(Y1701, ceiling(length(Y1701)*0.3))

Y1702 = grep('Y1702', rownames(Y))
Y1702_10 = sample(Y1702, ceiling(length(Y1702)*0.1))
Y1702_20 = sample(Y1702, ceiling(length(Y1702)*0.2))
Y1702_30 = sample(Y1702, ceiling(length(Y1702)*0.3))

Y1703 = grep('Y1703', rownames(Y))
Y1703_10 = sample(Y1703, ceiling(length(Y1703)*0.1))
Y1703_20 = sample(Y1703, ceiling(length(Y1703)*0.2))
Y1703_30 = sample(Y1703, ceiling(length(Y1703)*0.3))

trn10 = c(trn0, Y1608_10,Y1609_10, Y1612_10,Y1701_10,
          Y1702_10,Y1703_10)
  
trn20 = c(trn0, Y1608_20,Y1609_20, Y1612_20,Y1701_20,
          Y1702_20,Y1703_20)
  
trn30 = c(trn0, Y1608_30,Y1609_30, Y1612_30,Y1701_30,
          Y1702_30,Y1703_30)

tst10 <- seq_along(1:n)[-trn10]
tst20 <- seq_along(1:n)[-trn20]
tst30 <- seq_along(1:n)[-trn30]


# Calculate variance components ratio using training data
yNA = y
yNA[tst0] = NA
fm0 = fitBLUP(yNA,K=G)
GBLUP_0 = cor(fm0$u[tst0],y[tst0])

yNA = y
yNA[tst10] = NA
fm0 = fitBLUP(yNA,K=G)
GBLUP_10 = cor(fm0$u[tst0],y[tst0])

yNA = y
yNA[tst20] = NA
fm0 = fitBLUP(yNA,K=G)
GBLUP_20 = cor(fm0$u[tst0],y[tst0])

yNA = y
yNA[tst30] = NA
fm0 = fitBLUP(yNA,K=G)
GBLUP_30 = cor(fm0$u[tst0],y[tst0])

##### SSI

#0%
fm1 = SSI.CV(y,K=G,trn=trn0,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn0, tst=tst0,lambda=lambda)
SSI_0 = summary(fm2)$optCOR

#10%
fm1 = SSI.CV(y,K=G,trn=trn10,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn10, tst=tst10,lambda=lambda)
SSI_10 = summary(fm2)$optCOR

#20%
fm1 = SSI.CV(y,K=G,trn=trn20,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn20, tst=tst20,lambda=lambda)
SSI_20 = summary(fm2)$optCOR

#30%
fm1 = SSI.CV(y,K=G,trn=trn30,
             nCV=10,name="5 5CV", nfolds = 10)

lambda = summary(fm1)$optCOR["lambda"]

# Fit the index with the obtained lambda
fm2 = SSI(y,K=G,trn=trn30, tst=tst30,lambda=lambda)
SSI_30 = summary(fm2)$optCOR

## KA
yNA = y
yNA[tst0] = NA

h<- c(.02,1,5) # bandwidth kernels
KList<-list() 

for(i in 1:length(h)){ # Avering kernel
  KList[[i]]<-list(K=exp(-h[i]*D),model='RKHS')
}

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_preset0_",trait_name,"_rep_",rep))

KA_0 = cor(y[tst0],fmKA$yHat[tst0],
         use = 'pairwise.complete.obs')

yNA = y
yNA[tst10] = NA

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_preset10_",trait_name,"_rep_",rep))

KA_10 = cor(y[tst10],fmKA$yHat[tst10],
           use = 'pairwise.complete.obs')

#== 20 % 
yNA = y
yNA[tst20] = NA

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_preset20_",trait_name,"_rep_",rep))

KA_20 = cor(y[tst20],fmKA$yHat[tst20],
            use = 'pairwise.complete.obs')

# ======= 30%
yNA = y
yNA[tst30] = NA

fmKA<-BGLR(y=yNA,ETA=KList,
           nIter=12000,burnIn=2000,
           saveAt=paste0("KA_preset30_",trait_name,"_rep_",rep))

KA_30 = cor(y[tst30],fmKA$yHat[tst30],
            use = 'pairwise.complete.obs')


results = c(KA_0,KA_10,KA_20,KA_30,
            SSI_0,SSI_10,SSI_20,SSI_30)

names(results) = c('KA_0','KA_10','KA_20','KA_30',
                   'SSI_0', 'SSI0_MSE', 'SSI0_df','SSI0_lambda',
                   'SSI_10','SSI10_MSE', 'SSI10_df','SSI10_lambda',
                   'SSI_20','SSI20_MSE', 'SSI20_df','SSI20_lambda',
                   'SSI_30','SSI30_MSE', 'SSI30_df','SSI30_lambda'
                   )
save(results,
     file=paste0("results_predset_",trait_name,"_rep_",rep,"_2023.RData"))
