#Set directory
setwd("C:/Users/dcho0009/Desktop/multiomics_trial")

#Load libraries
library(DEqMS)
library(dplyr)

#Proteomics
#Load data
prot_s <- read.csv('proteomics_SCx.csv', stringsAsFactors = F)

#Cleaning data
#removing the first column, making second column rownames
prot_s[,2] <- gsub('_RAT','', prot_s[,2])
prot_s.c <- prot_s[,-(1:2)] 
rownames(prot_s.c) <- prot_s[,2]

#Need to change the labels of the samples to reflect strain
strain <- c("GAERS.1", "GAERS.2", "GAERS.3", "GAERS.4", "GAERS.5", "GAERS.6", "NEC.1", "NEC.2", "NEC.3", "NEC.4", "NEC.5", "NEC.6")

colnames(prot_s.c) <- strain

#Checking for any NA's or 0s
range(prot_s.c) #3.925928e+02 8.137162e+07

#Filtering out low variance proteins (keeping top 5000 most variable proteins)
var_p_s <- apply(prot_s.c, 1, var) # X is of size number of samples x number of probes
hist(var_p_s) # gives you an idea of the variance
#Appending variance to protein dataframe
var_prot_s <- cbind(prot_s.c, var_p_s)
colnames(var_prot_s)[13] <- "var"
p_s_top <- var_prot_s[var_prot_s[,13] >= quantile(var_prot_s[,13], 0.07),] #selecting top 93% most variable proteins
#Removing the variance column as we don't need it anymore
p_s_final <- p_s_top[,-13]

boxplot(p_s_final, ylim = c(4.541467e+02, 3e+05))

#transforming
p_log <- log2(p_s_final)

boxplot(p_log, main = "Proteomics - log")
#title("Proteomics - log")

#scaling data
p_scale <- as.data.frame(scale(p_log))

colMeans(p_scale)
apply(p_scale, 2, sd)

boxplot(p_scale)
title("Proteomics - log + scale")

#Distribution of protein
p_prot <- as.data.frame(t(p_scale))
boxplot(p_prot[,1:20])
title("Protein 1 - 20 raw distribution")

#replacing outliers
#https://stackoverflow.com/questions/13339685/how-to-replace-outliers-with-the-5th-and-95th-percentile-values-in-r
#https://rpubs.com/frasermyers/627592
cap <- function(x){
  quantiles <- quantile(x, c(0.003, 0.997)) #quantiles of interest from variable
  x[x < quantiles[1]] <- quantiles[1] #replace points below 0.3% with the 0.3 quantile value
  x[x > quantiles[2]] <- quantiles[2] #replace points above 99.7% with the 99.7 quantile value
  x # return the whole dataset
}

#p_G1 <- cap(p_scale$GAERS.1)
p_final <- as.data.frame(apply(p_scale, 1, cap))
boxplot(p_final[,1:20])
title("Protein 1 - 20 Winsorised Distribution")

#Checking if function has worked
apply(p_scale, 2, quantile, probs = c(0.003, 0.997))
#apply(p_final, 2, quantile, probs = c(0.003, 0.997))
apply(p_final, 2, range)

#Using winzorise function
#library(DescTools)
#win <- function(x){
#  Winsorize(x,probs = c(0.003, 0.997))
#  x
#}
#p_win <- as.data.frame(apply(p_scale, 1, win))
#checking for similarities
#p_same <- semi_join(p_prot, p_win)


#Checking normality
#shapiro.test(p_final$GAERS.1)

par(mfrow = c(2, 3))
hist(p_final$GAERS.1, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.1")
hist(p_final$GAERS.2, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.2")
hist(p_final$GAERS.3, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.3")
hist(p_final$GAERS.4, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.4")
hist(p_final$GAERS.5, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.5")
hist(p_final$GAERS.6, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values GAERS.6")


par(mfrow = c(2, 3))
hist(p_final$NEC.1, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.1")
hist(p_final$NEC.2, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.2")
hist(p_final$NEC.3, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.3")
hist(p_final$NEC.4, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.4")
hist(p_final$NEC.5, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.5")
hist(p_final$NEC.6, xlab = "Log+scaled+outlier replaced intensity value", 
     ylab = "Number of proteins", main = "Distribution of intensity values NEC.6")

#transposing data
p_trn <- t(p_scale)

#Can make PCA plots as well
library(ggplot2)
library(ggforce)
library(tidyverse)
library(ggrepel)

group <- c("GAERS","GAERS","GAERS","GAERS","GAERS","GAERS","NEC","NEC","NEC","NEC","NEC","NEC")
p_group <- p_final %>%
  as.data.frame() %>%
  add_column(group, .before = "A0A023IKF6") 

p_pca  <- prcomp(p_group[,-1]) #PCA components
p_comp <- as.data.frame(p_pca$x) #isolate components
new_p <- cbind(p_group,p_comp[,c(1,2)]) #Extract components and bind to original data
par(mfrow = c(1, 1))
ggplot(new_p, aes(x=PC1, y=PC2, col = group, fill = group)) +
  stat_ellipse(geom = "polygon", col = "black", alpha =0.5) +
  geom_point(shape=21, col="black") +
  geom_text_repel(aes(label = rownames(p_group)), colour = 'black') +
  ggtitle("Proteomics")

#Save final file
write.csv(p_final, "prot_SCx_log_scale_OutReplaceRow.csv")


#Metabolomics
#Load in data
met_s <- read.csv('metabolomics_SCx.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')

#Cleaning
#For SCx, we are not including GAERS_629 and GAERS_630
met_s <- met_s[, -(5:6)]
met_s <- met_s[-1,]
met_s.c <- met_s[,-1]
rownames(met_s.c) <- met_s[,1]

#Renaming columns
strain <- c("GAERS.1", "GAERS.2", "GAERS.3", "GAERS.4", "GAERS.5", "GAERS.6", "NEC.1", "NEC.2", "NEC.3", "NEC.4", "NEC.5", "NEC.6")
colnames(met_s.c) <- strain

#Converting columns to numeric
met_s.c[,1:12] <- as.numeric(as.matrix(met_s.c[,1:12]))
sapply(met_s.c, class)

#Transformation
m_log10 <- log10(met_s.c + 1)
range(m_log10)

boxplot(m_log10)
title("Metabolomic - log transformation")

m_log_t <- as.data.frame(t(m_log10))
boxplot(m_log_t[,1:20])
title("Metabolite 1 - 2 Distribution")

#replacing outliers
cap <- function(x){
  quantiles <- quantile(x, c(0.003, 0.997))
  x[x < quantiles[1]] <- quantiles[1]
  x[x > quantiles[2]] <- quantiles[2]
  x
}

m_final <- as.data.frame(apply(m_log10, 1, cap))
boxplot(m_final)
title("Metabolomic - Outlier Replacement")

apply(m_log10, 2, quantile, probs = c(0.003, 0.997))

#apply(p_final, 2, quantile, probs = c(0.003, 0.997))
apply(m_final, 2, range)

#Testing for normality
#shapiro.test()

par(mfrow = c(2, 3))
hist(m_final$GAERS.1, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.1")
hist(m_final$GAERS.2, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.2")
hist(m_final$GAERS.3, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.3")
hist(m_final$GAERS.4, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.4")
hist(m_final$GAERS.5, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.5")
hist(m_final$GAERS.6, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values GAERS.6")


par(mfrow = c(2, 3))
hist(m_final$NEC.1, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.1")
hist(m_final$NEC.2, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.2")
hist(m_final$NEC.3, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.3")
hist(m_final$NEC.4, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.4")
hist(m_final$NEC.5, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.5")
hist(m_final$NEC.6, xlab = "Log+outlier replaced intensity value", 
     ylab = "Number of metabolites", main = "Distribution of intensity values NEC.6")

#Transposing dataset
m_trn <- as.data.frame(t(m_log10))

#PCA plot
group <- c("GAERS","GAERS","GAERS","GAERS","GAERS","GAERS","NEC","NEC","NEC","NEC","NEC","NEC")
m_group <- m_final %>%
  as.data.frame() %>%
  add_column(group, .before = "C3H5NO2") 

m_pca  <- prcomp(m_group[,-1]) #PCA components
m_comp <- as.data.frame(m_pca$x) #isolate components
new_m <- cbind(m_group,m_comp[,c(1,2)]) #Extract components and bind to original data
par(mfrow = c(1, 1))
ggplot(new_m, aes(x=PC1, y=PC2, col = group, fill = group)) +
  stat_ellipse(geom = "polygon", col = "black", alpha =0.5) +
  geom_point(shape=21, col="black") +
  geom_text_repel(aes(label = rownames(m_group)), colour = 'black') +
  ggtitle("Metabolomics")

#Saving
write.csv(m_final, 'met_SCx_log_noOutReplace.csv')
