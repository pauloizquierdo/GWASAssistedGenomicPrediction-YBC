
# Genome-wide association and genomic prediction for Fe and Zn concentration and Fe bioavailability in a yellow bean collection of dry beans
# Paulo Izquierdo1, Rie Sadohara1, Jason Wiesinger2, Raymond Glahn2, Carlos Urrea3, Karen Cichy1,4

# 1 Department of Plant, Soil and Microbial Sciences, Michigan State University, East Lansing, MI, USA.
# 2 USDA – ARS, Robert W. Holley Center for Agriculture and Health, Ithaca, NY, USA.
# 3 Department of Agronomy and Horticulture, Panhandle Research and Extension Center, University of Nebraska-Lincoln, Scottsbluff, NE, USA.
# 4 USDA-ARS, Sugarbeet and Bean Research Unit, East Lansing, MI, USA.


library(tidyverse)
library(readxl)
library(phia)
library(lme4)
library(agridat)
library(rgl)
library(tidyverse)
library(plyr)
library(plot3Drgl)
library(statgenGxE)
library(corrplot)
library(RColorBrewer)
library(beanplot)
library(GAPIT3)
library(pheatmap)

rm(list = ls())

######### Load data
load("YBC_GWASAssistedGP.RData")

###============== Correlations

cor_pearson <- cor(YBC_GWASAssistedGP$pheno[,-c(1,20)], y = NULL, 
              use = "pairwise.complete.obs",
              method = c("pearson"))

testRes <- cor.mtest(YBC_GWASAssistedGP$pheno[,-c(1,20)], conf.level = 0.95)

colfunc <- colorRampPalette(c("red", 
                              'white', 
                              "cornflowerblue"))


corrplot(cor_pearson, method="color", 
         col = colfunc(10),
          addCoef.col = "black", 
          tl.cex = 1,
          number.cex = 0.9, 
          tl.col="black", pch.cex = 1,
         # p.mat=testRes$p,
          type="lower",  
          diag = F, 
          #Add coefficient of correlation
          #tl.srt=0, #Text label color and rotation
          #sig.level = c(0.001,0.01,0.05), 
          insig='label_sig',
          outline = T, 
          tl.offset = 0.5
          )

###============== Bean plots

# beanplot minerals

par(mfrow =c(3,2)) 
par(mar=c(2,4,1,1))

beanplot(YBC_GWASAssistedGP$pheno$Fe_MI18,
         YBC_GWASAssistedGP$pheno$Fe_MI19,
         ylab = "Fe (µg/g)", 
         log="", 
         ylim = c(20,120),
         names = 'MI',
         side="both",col = list("gold", "grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80")) 

beanplot(YBC_GWASAssistedGP$pheno$Fe_NE18,
         YBC_GWASAssistedGP$pheno$Fe_NE19,
         ylab = "Fe (µg/g)", 
         ylim = c(20,120),
         log="", 
         names = 'NE',
         side="both",col = list("gold", "grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80")) 
         
beanplot(YBC_GWASAssistedGP$pheno$Zn_MI18, 
         YBC_GWASAssistedGP$pheno$Zn_MI19,
         ylab = "Zn (µg/g)", 
         log="", 
         names = 'MI',
         side="both",col = list("gold", "grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80")) 
  
         
beanplot(YBC_GWASAssistedGP$pheno$Zn_NE18, 
         YBC_GWASAssistedGP$pheno$Zn_NE19,
         ylab = "Zn (µg/g)", 
         log="", 
         names = 'NE',
         side="both",col = list("gold", "grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80")) 


beanplot(YBC_GWASAssistedGP$pheno$FeBio_MI18, 
         YBC_GWASAssistedGP$pheno$FeBio_MI19,
         ylab = "FeBio (% ferritin/control)", 
         log="", 
         names = 'MI',
         side="both",col = list("gold", "grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80")) 


# Agronomic traits
par(mfrow =c(1,2)) 

beanplot(YBC_GWASAssistedGP$pheno$Yield_MI18,YBC_GWASAssistedGP$pheno$Yield_MI19,
         YBC_GWASAssistedGP$pheno$Yield_NE18,YBC_GWASAssistedGP$pheno$Yield_NE19,
         ylab = "Yield (kg/ha)",
         names = c("Mi", "NE"),
         log="", 
         side="both",
         col = list("gold", c("grey80")),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80"))

beanplot(YBC_GWASAssistedGP$pheno$SW_MI18,YBC_GWASAssistedGP$pheno$SW_MI19,
         YBC_GWASAssistedGP$pheno$SW_NE18,YBC_GWASAssistedGP$pheno$SW_NE19,
         ylab = "SW (g)", 
         names = c("Mi", "NE"),
         log="", 
         side="both",
         col = list("gold","grey80"),
         what=c(0,1,1,0))
         legend("topright", bty="n",c("2018", "2019"),
         fill = c("gold", "grey80"))
###


###============== Variance component analysis 

long_pheno <- gather(YBC_GWASAssistedGP$pheno, traits,
                     value, Yield_MI18:Zn_NE19, 
                     factor_key = T)                   

long_pheno <- long_pheno %>% separate(traits,c('trait', 'loc_year'), sep=('_'))

# split location and year using the last 2 characters
long_pheno$location <- substr(long_pheno$loc_year, 1, 2)
long_pheno$year <- substr(long_pheno$loc_year, 3, 4)

long_pheno$taxa = as.factor(long_pheno$taxa)
long_pheno$trait  = as.factor(long_pheno$trait)

Zn = long_pheno %>% filter(trait == 'Zn') 
Fe = long_pheno %>% filter(trait == 'Fe')
YD = long_pheno %>% filter(trait == 'Yield')
FeBio = long_pheno %>% filter(trait == 'FeBio')
SW = long_pheno %>% filter(trait == 'SW') 

traits <- list(Fe,Zn,YD,SW,FeBio)
names(traits) <- c("Fe","Zn","YD","SW", "FeBio")

resultsmodel <- list()

# Run model for Fe,Zn,YD,SW - Location:Year = TRUE (2 years, 2 locations)
for (i in 1:4) {
  TD <- statgenSTA::createTD(data = as.data.frame(traits[[i]]), 
                                genotype = 'taxa', 
                                year =  'year', 
                                loc = "location",
                                trial = 'loc_year'
                                )
  resultsmodel[[i]] <- gxeVarComp(TD = TD, 
                           trait = "value",
                           locationYear = TRUE)
}

# Run model for FeBio - Location:Year = False (2 years, 1 locations)

TD <- statgenSTA::createTD(data = as.data.frame(traits[[5]]), 
                                genotype = 'taxa', 
                                year =  'year', 
                                loc = "location",
                                trial = 'loc_year'
                                )

resultsmodel[[5]] <- gxeVarComp(TD = TD, 
                           trait = "value",
                           locationYear = FALSE)

names(resultsmodel) <- c("Fe","Zn","YD","SW", "Fe_bio")


#### Residuals
predictions_gxe <- list()

for (i in 1:length(resultsmodel)) {
  predictions_gxe[[i]] <- predict(resultsmodel[[i]], predictLevel = "trial")
}
names(predictions_gxe) <- c("Fe","Zn","YD","SW", "Fe_bio")
colnames(predictions_gxe[[1]])

# estimate the residuals
residuals <- list()

for(i in 1:length(predictions_gxe)){
  traits[[i]]$id_trial <- paste(traits[[i]]$taxa, 
                                   traits[[i]]$loc_year, sep = "_")

  predictions_gxe[[i]]$id_trial <- paste(predictions_gxe[[i]]$genotype, 
                                         predictions_gxe[[i]]$trial, 
                                         sep = "_")                                   
  predGenoTrial <- predictions_gxe[[i]][,3:4]
  predGenoTrial <- merge(traits[[i]], predGenoTrial, by = "id_trial")
  head(predGenoTrial)
  residuals[[i]] <- predGenoTrial$value - predGenoTrial$predictedValue
}

# Plot residuals - Variance component analysis 

# remove NA
for(i in 1:length(residuals)){
  residuals[[i]] <- residuals[[i]][!is.na(residuals[[i]])]
}

names(residuals) <- c("Fe","Zn","Yield","SW", "FeBio")

par(mfrow = c(5,2))
for(i in 1:length(residuals)){
  hist(residuals[[i]], breaks=50, main = paste0("Histogramns of residuals - ", names(residuals)[i]), xlab=NULL)
  qqnorm(residuals[[i]], main = paste0("Normal Q-Q Plot - ", names(residuals)[i]))
  qqline(residuals[[i]], col="red")
}

###============== Shapiro test

shapiro_table <- data.frame(Statistic = numeric(ncol(YBC_GWASAssistedGP$pheno[-c(1,20)])), 
                            P_Value = numeric(ncol(YBC_GWASAssistedGP$pheno[,-c(1,20)])), 
                            Normality_Status = character(ncol(YBC_GWASAssistedGP$pheno[,-c(1,20)])), 
                            stringsAsFactors = FALSE)

rownames(shapiro_table) <- colnames(YBC_GWASAssistedGP$pheno[-c(1,20)])

for (i in 1:ncol(YBC_GWASAssistedGP$pheno[-c(1,20)])) {
  test_result <- shapiro.test(YBC_GWASAssistedGP$pheno[-c(1,20)][,i])
  shapiro_table[i,1] <- test_result$statistic
  shapiro_table[i,2] <- round(test_result$p.value, 4)
  if (shapiro_table[i,2] < 0.05) {
    shapiro_table[i,3] <- "not normal"
  } else {
    shapiro_table[i,3] <- "normal"
  }
}
# Data that is no-normal distributed
# Yield_NE18, Yield_NE19, SW_MI19, FeBio_MI18, FeBio_MI19, Zn_MI19, Zn_NE18, Zn_NE19


# qqplot for all traits
par(mfrow = c(5,4))

for(i in 2:11){
  hist(YBC_GWASAssistedGP$pheno[,i], breaks=10, 
       main = colnames(YBC_GWASAssistedGP$pheno)[i], 
       xlab=NULL)
  qqnorm(YBC_GWASAssistedGP$pheno[,i], main = colnames(YBC_GWASAssistedGP$pheno)[i])
  qqline(YBC_GWASAssistedGP$pheno[,i], col="red")
}

for(i in 12:19){
  hist(YBC_GWASAssistedGP$pheno[,i], breaks=10, 
       main = colnames(YBC_GWASAssistedGP$pheno)[i], 
       xlab=NULL)
  qqnorm(YBC_GWASAssistedGP$pheno[,i], main = colnames(YBC_GWASAssistedGP$pheno)[i])
  qqline(YBC_GWASAssistedGP$pheno[,i], col="red")
}

# Altgougth the data is not normal distributed (shapiro test < 0.05), all traits but FeBio are close to normal distribution.

###============== Transform non-normal distributed data

# Rank-based inverse normal transform (INT)

phenotype <- YBC_GWASAssistedGP$pheno[,-c(1,20)]
rownames(phenotype) <- YBC_GWASAssistedGP$pheno$taxa
#create matrix for transformed phenotype
INT_phenotype <- matrix(NA, nrow = nrow(phenotype), ncol = ncol(phenotype))
rownames(INT_phenotype) <- rownames(phenotype)
colnames(INT_phenotype) <- colnames(phenotype)[1:ncol(phenotype)]
dim(INT_phenotype)
alpha <- 1/2 #  for larger sample sizes (>200), an alpha of 1/2 is generally used. 

for(i in 1:ncol(phenotype)){
    n <- sum(!is.na(phenotype[,i]))
    ranks <- rank(phenotype[,i])
    transfor <- qnorm((ranks - alpha) / (n - 2 * alpha + 1))
    INT_phenotype[,i] <- transfor
}

# shapiro test transformed phenotype

shapito_table_INT <- data.frame(Statistic = numeric(ncol(INT_phenotype)), 
                            P_Value = numeric(ncol(INT_phenotype)), 
                            Normality_Status = character(ncol(INT_phenotype)), 
                            stringsAsFactors = FALSE)

rownames(shapito_table_INT) <- colnames(phenotype)

for (i in 1:ncol(INT_phenotype)) {
  test_result <- shapiro.test(INT_phenotype[,i])
  shapito_table_INT[i,1] <- test_result$statistic
  shapito_table_INT[i,2] <- round(test_result$p.value, 4)
  if (shapito_table_INT[i,2] < 0.05) {
    shapito_table_INT[i,3] <- "not normal"
  } else {
    shapito_table_INT[i,3] <- "normal"
  }
}
# All traits are normal distributed after INT

###============== GWAS
INT_phenotype <- as.data.frame(INT_phenotype)
phenotype$taxa <- rownames(phenotype)
INT_phenotype$taxa <- rownames(INT_phenotype)

phenotype <- phenotype[,c(19,1:18)] # move taxa to the first column
INT_phenotype <- INT_phenotype[,c(19,1:18)]

# Run GWAS only in YBC population
index_ybc <- grep("YBC", phenotype$taxa)
phenotype <- phenotype[index_ybc,]
INT_phenotype <- INT_phenotype[index_ybc,]

# GWAS - Untransformed data
setwd("GWAS_untransformed")

myGAPIT <- GAPIT(
  Y = phenotype[,1:3],
  G = YBC_GWASAssistedGP$hapmap, 
  PCA.total=3,
  model= c('FarmCPU')
)

# read all files with GWAS.Results.csv 
files_gwas <- list.files( pattern = "GWAS.Results.csv", 
                         full.names = TRUE)

# read all files with GWAS.Results.csv
gwas_results <- list()
for (i in 1:length(files_gwas)) {
  gwas_results[[i]] <- read.csv(files_gwas[i], header = TRUE, sep = ",", stringsAsFactors = FALSE)
  gwas_results[[i]]$trait <- str_extract(files_gwas[[i]], "(?<=FarmCPU\\.).*?(?=\\.GWAS)")
}

# rbind all gwas results
gwas_results <- do.call(rbind, gwas_results)

# remove all SNPs with pvalue > bonferoni threshold (0.05/number of SNPs)
gwas_results_fil <- gwas_results %>% 
  filter(P.value < 0.05/(nrow(YBC_GWASAssistedGP$hapmap)-1))

# GWAS - Transformed data
setwd("../GWAS_INT")

myGAPIT <- GAPIT(
  Y = INT_phenotype[,c(1,3:4)],
  G = YBC_GWASAssistedGP$hapmap, 
  PCA.total=3,
  model= c('FarmCPU'),
)

# read all files with GWAS.Results.csv 
files_gwas_int <- list.files( pattern = "GWAS.Results.csv", 
                         full.names = TRUE)

# read all files with GWAS.Results.csv
gwas_results_int <- list()

for (i in 1:length(files_gwas_int)) {
  gwas_results_int[[i]] <- read.csv(files_gwas_int[i], header = TRUE, sep = ",", stringsAsFactors = FALSE)
  gwas_results_int[[i]]$trait <- str_extract(files_gwas_int[[i]], "(?<=FarmCPU\\.).*?(?=\\.GWAS)")
}

# rbind all gwas results
gwas_results_int <- do.call(rbind, gwas_results_int)

# remove all SNPs with pvalue > bonferoni threshold (0.05/number of SNPs)
gwas_results_int_fil <- gwas_results_int %>% 
  filter(P.value < 0.05/(nrow(YBC_GWASAssistedGP$hapmap)-1))


###============== Population structure

# All samples
 PCA_allsamples <- prcomp(YBC_GWASAssistedGP$geno, scale = TRUE)
 sum_pca_all <- summary(PCA_allsamples)
 sum_pca_all$importance[,1:5]


 PCA_allsamples <- PCA_allsamples$x[,1:5]
 rownames(YBC_GWASAssistedGP$pheno) <- YBC_GWASAssistedGP$pheno$taxa

PCA_allsamples <- merge(YBC_GWASAssistedGP$pheno[,c(1,20)], PCA_allsamples, by = 'row.names')

PCA_allsamples$Subpop_Sadohara_etal2022 <- factor(PCA_allsamples$Subpop_Sadohara_etal2022,
                    levels = c("Andean","MA","Admix"))

pca1 <- round(sum_pca_all$importance[2,1]*100,1)
pca2 <- round(sum_pca_all$importance[2,2]*100,1)

ggplot(PCA_allsamples,aes(PC1,PC2, fill=Subpop_Sadohara_etal2022)) +
  geom_point(size=3, alpha = 0.8, shape= 21, ) +
  labs(x = paste0("PC1 (",pca1,"%)"), size= 15, y= paste0("PC2 (",pca2,"%)")) +
  theme(axis.title=element_text(size=14))+ 
  theme_light()  + 
  scale_fill_discrete(name = "Gene pool",
                      labels = c('Andean', 
                                 'Middle American', 'Admixture'),
                      na.translate = FALSE)

# Kinship all samples
M <- as.matrix(YBC_GWASAssistedGP$geno)
G <- tcrossprod(scale(M))/ncol(M) # genomic relationship

Andean_acc <- grep("YBC",rownames(G))
AY_acc <- grep("Y1",rownames(G))

trn_r <- factor(rep(c("Training"), length(Andean_acc)))
tst_r <- factor(rep(c("Prediction"), length(AY_acc)))

annotation_col_row <- data.frame(Set = c(trn_r, tst_r))

rownames(annotation_col_row) <-  rownames(G)

callback = function(hc, mat){
sv = svd(t(mat))$v[,1]
dend = reorder(as.dendrogram(hc), wts = sv)
as.hclust(dend)
}

colfunc <- colorRampPalette(c("khaki1", "red"))

ann_colors = list(Set = c(Training = 'cornflowerblue', 
                          Prediction = 'greenyellow'))

pheatmap(G, clustering_callback = callback,
         color = colfunc(100),
         width = 8, height = 8,
         annotation_col = annotation_col_row,
         annotation_row = annotation_col_row,
         annotation_colors = ann_colors,
         #euclidean = T,
         show_rownames = F, show_colnames = F,
         #clustering_distance_rows = 'euclidean',
         # cluster_cols = 'euclidean',
         #cluster_rows = FALSE, cluster_cols = F,
         gaps_row = 275,
         gaps_col = 275) 

# Andean accessions - PCA

index_Andean <- which(YBC_GWASAssistedGP$pheno$Subpop_Sadohara_etal2022  == "Andean")
andean_geno <- YBC_GWASAssistedGP$pheno[index_Andean,]

index_Andean <- which(rownames(YBC_GWASAssistedGP$geno) %in% andean_geno$taxa)
Andean <- YBC_GWASAssistedGP$geno[index_Andean,]

Andean <- Andean[,apply(Andean,2,function(x) length(unique(x)))>1] # remove monomorphic SNPs

pca_andean <- prcomp(Andean, scale = TRUE)
sum_pca <- summary(pca_andean)
pca1_an <- round(sum_pca$importance[2,1]*100,1)
pca2_an <- round(sum_pca$importance[2,2]*100,1)

# extract pcs
PC10_andean <- as.data.frame(pca_andean$x[,1:5])
PC10_andean$Parents <- rownames(PC10_andean)
parents <- c("YBC113","YBC122","YBC126","YBC129","YBC137","YBC190")
index_parents <- which( !PC10_andean$Parents %in%  parents )
PC10_andean$Parents[index_parents] <- "Andean"

PC10_andean$Parents <- factor(PC10_andean$Parents, 
                              levels = c("Andean","YBC113","YBC122","YBC126","YBC129","YBC137","YBC190"))
# order by parents
PC10_andean <- PC10_andean[order(PC10_andean$Parents),]
PC10_andean[198,] <- PC10_andean[193,] # move this genotype to the end to avoid overlapping

ggplot(PC10_andean, aes(PC1, PC2, fill=Parents)) +
  geom_point(size=3, alpha = 0.9, shape= 21) +
  labs(x = paste0("PC1 (",pca1_an,"%)"), y= paste0("PC1 (",pca2_an,"%)")) +
  theme(axis.title = element_text(size = 14)) + 
  theme_light() +
  scale_fill_manual(name="Parents", 
                    values=c("Andean"="grey99", 
                             "YBC113"="#dd5050", 
                             "YBC122"="#3434def9", 
                             "YBC126"="black", 
                             "YBC129"="#faec72", 
                             "YBC137"="#00ff1a", 
                             "YBC190"="#91247b"),
                    labels=c("YBC113", "YBC122", "YBC126", "YBC129", "YBC137", "YBC190"),
                    breaks=c("YBC113", "YBC122", "YBC126", "YBC129", "YBC137", "YBC190"))

###============== Phenotype - families

phenotype <- YBC_GWASAssistedGP$pheno

index_ybc <- grep("YBC",rownames(phenotype))
pheno <- phenotype
pheno$Family <- "NA"
pheno$Family[index_ybc] <- "YBC"
pheno$Family[grep("Y1608",rownames(pheno))] <- "Y1608"
pheno$Family[grep("Y1609",rownames(pheno))] <- "Y1609"
pheno$Family[grep("Y1612",rownames(pheno))] <- "Y1612"
pheno$Family[grep("Y1701",rownames(pheno))] <- "Y1701"
pheno$Family[grep("Y1702",rownames(pheno))] <- "Y1702"
pheno$Family[grep("Y1703",rownames(pheno))] <- "Y1703"

colnames(pheno)

index <- grep("MI19", colnames(pheno))
pheno <- pheno[,c(index,21)]
pheno <- pheno[,-3] # YA dont have FeBio
pheno <- pheno %>% filter(Family != "NA")

colnames(pheno) <- c("Yield (kg/ha)", "SW (g)", "Fe (µg/g)", "Zn (µg/g)", "Family")

pheno_tidy <- gather(pheno, key = "Trait", value = "value", 1:4)

my_colors <-  c("#999999", "#E69F00",
 "#56B4E9", "#009E73", "#0072B2",
 "#D55E00","#F0E442")

pheno_tidy$Trait <- factor(pheno_tidy$Trait, 
                           levels = c( "Fe (µg/g)", "Zn (µg/g)","Yield (kg/ha)", "SW (g)"))


ggplot(pheno_tidy, 
   aes(x=Family, y=value, fill=Family)) + 
    geom_boxplot() +
    facet_grid(Trait ~ ., scales = "free_y") +
    xlab(label = "") + 
    ylab(label="") +
    theme(axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, 
             hjust = 1)) +
    scale_fill_manual(values = my_colors)


###============== Seed Color vs minerals

# Data
# Seed color data was obtained from Sadohara et al., 2021 (https://doi.org/10.1002/tpg2.20173). 
# Cooking time, population structure, and country of origin information were sourced from Sadohara et al., 2022."

YBC_color <- read_csv("YBC_phenotype_Color_CT.csv")
commoncol <- sort(table(YBC_color$Seed_type), decreasing = TRUE)


# chose the 9 most common seed colors
YBC_color <- YBC_color %>% filter(Seed_type %in% names(commoncol[1:10]))

YBC_color$Seed_type =  factor(YBC_color$Seed_type,
                      levels = c("Beige and brown stripes",
                                 "Amarillo (dk)", "Canary",
                                 "Green-yellow", "Manteca",
                                 "Mayocoba", "Brown",
                                 "Beige","Amarillo (lt)",
                                 "White"))

# Seed type vs FeBio MI2018
ggplot(YBC_color, aes(x = reorder(Seed_type,FeBio_MI18,
                              na.rm=TRUE), 
                  y = FeBio_MI18, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("FeBio - MI 2018 (% of control)") +
  xlab(NULL) +
  ylim(c(10,135)) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))

# Seed type vs FeBio MI2019
ggplot(YBC_color, aes(x = reorder(Seed_type,FeBio_MI19,
                              na.rm=TRUE), 
                  y = FeBio_MI19, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("FeBio - MI 2019 (% of control)") +
  xlab(NULL) +
  ylim(c(10,135)) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))

# Seed type vs Fe MI2018

ggplot(YBC_color, aes(x = reorder(Seed_type,Fe_MI18,
                              na.rm=TRUE), 
                  y = Fe_MI18, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("Fe - MI 2018 (µg/g)") +
  xlab(NULL) +
  ylim(c(40,120)) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))

# Seed type vs Fe MI2019

ggplot(YBC_color, aes(x = reorder(Seed_type,Fe_MI19,
                              na.rm=TRUE), 
                  y = Fe_MI19, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("Fe - MI 2019 (µg/g)") +
  xlab(NULL) +
  ylim(c(40,120)) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))

# Seed type vs Zn MI2018

ggplot(YBC_color, aes(x = reorder(Seed_type,Zn_MI18,
                              na.rm=TRUE), 
                  y = Zn_MI18, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("Zn - MI 2018 (µg/g)") +
  xlab(NULL) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))

# Seed type vs Zn MI2019
ggplot(YBC_color, aes(x = reorder(Seed_type,Zn_MI19,
                              na.rm=TRUE), 
                  y = Zn_MI19, 
                  fill = Seed_type)) +
  geom_boxplot(show.legend = FALSE) +
  ylab("Zn - MI 2019 (µg/g)") +
  xlab(NULL) +
  theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 20, 
                                      hjust = 1)) +
  scale_fill_manual(values = c('beige',
                              'darkorange2', 
                              'gold',
                              'lavenderblush3',
                              'khaki1', 'yellow',
                              'darkgoldenrod4',
                              'burlywood1',
                              'lightgoldenrod',
                              'white'))