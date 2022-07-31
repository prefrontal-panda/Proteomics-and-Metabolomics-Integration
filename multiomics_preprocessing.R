#This is for the pre-processing of proteomics and metabolomics data

##Set directory
setwd("C:/Users/dcho0009/Desktop/multiomics_trial")

#Because of the presence of missing values in the metabolomics data, we will have to perform imputation.
#The NIPALS algorithm implemented by the mixOmics package will be used.
#library(mixOmics)
#We also need to normalise the data. For proteomics, the DEqMS package will be used.
#For metabolomics, we will use quantile normalisation via the preprocessCore package
#BiocManager::install("") 
library(DEqMS)
library(preprocessCore)
library(dplyr)

##Proteomics
#Loading in the data
prot_s <- read.csv('proteomics_SCx.csv', stringsAsFactors = F)

#Cleaning up the data
#removing the first column, making second column rownames
prot_s[,2] <- gsub('_RAT','', prot_s[,2])
prot_s.c <- prot_s[,-(1:2)] 
rownames(prot_s.c) <- prot_s[,2]

#Need to change the labels of the samples to reflect strain
strain <- c("GAERS.1", "GAERS.2", "GAERS.3", "GAERS.4", "GAERS.5", "GAERS.6", "NEC.1", "NEC.2", "NEC.3", "NEC.4", "NEC.5", "NEC.6")

colnames(prot_s.c) <- strain

#Checking for any NA's or 0s
range(prot_s.c)
#3.925928e+02 8.137162e+07

#Filtering out low variance proteins
#It's a big dataset - keep the top 5000 most variable proteins.
var_p_s <- apply(prot_s.c, 1, var) # X is of size number of samples x number of probes
hist(var_p_s, xlim = c(0,7)) # gives you an idea of the variance
#Appending variance to protein dataframe
var_prot_s <- cbind(prot_s.c, var_p_s)
colnames(var_prot_s)[13] <- "var"
# https://www.biostars.org/p/291453/
p_s_top <- var_prot_s[var_prot_s[,13] >= quantile(var_prot_s[,13], 0.07),] #selecting top 93% of most variable proteins
#Removing the variance column as we don't need it anymore
p_s_final <- p_s_top[,-13]

boxplot(p_s_final, ylim = c(4.541467e+02, 3e+05))

#transform the data
p_sl <- log2(p_s_final)

p_sl_t <- t(p_sl)
write.csv(p_sl_t, "prot_SCx_log.csv")

boxplot(p_sl)
title("Prot_log transformed")

#then normalise the data for analysis
p_sn <- equalMedianNormalization(p_sl)

boxplot(p_sn)
title("Proteomics log2 + Median normalised")

#Then, transpose the data file
prot_s.t <- t(p_sn)

#Finally, save the file. 
#We save as a csv first then turn it into a tsv via tab delimited .txt file then change the file type to .tsv
write.csv(prot_s.t, 'prot_SCx_log_norm_filt.csv')
          

##Metabolomics
#Loading
met_s <- read.csv('metabolomics_SCx.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
#met_t <- read.csv('metabolomics_Tha.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')

#Removing the samples that are not included in the analysis and making rownames
#For SCx, we are not including GAERS_629 and GAERS_630
met_s <- met_s[, -(5:6)]
met_s <- met_s[-1,]
met_s.c <- met_s[,-1]
rownames(met_s.c) <- met_s[,1]
#For thalamus, GAERS_629. This has already been removed.

#met_t.c <- met_t[-1,]
#rownames(met_t.c) <- met_t.c[,1]
#met_t.c <- met_t.c[,-1]

#Renaming columns
strain <- c("GAERS.1", "GAERS.2", "GAERS.3", "GAERS.4", "GAERS.5", "GAERS.6", "NEC.1", "NEC.2", "NEC.3", "NEC.4", "NEC.5", "NEC.6")

colnames(met_s.c) <- strain

#Converting to numeric
#char_columns_scx <- sapply(met_s.c, is.character) # Identify character columns
#met_s.c[ , char_columns_scx] <- as.data.frame(apply(met_s.c[ , char_columns_scx], 2, as.numeric)) #Recode characters as numeric
#sapply(met_s.c, class) # Print classes of all columns
#This is an easier method
met_s.c[,1:12] <- as.numeric(as.matrix(met_s.c[,1:12]))
sapply(met_s.c, class)

#Imputation
range(met_s.c)

boxplot(met_s.c, ylim = c(0, 5e+6)) #right skewed data
title("Metabolomics raw")
#met_s.c[met_s.c == 0] <- NA #Converting 0s to NAs
#range(met_t.c)
#met_t.c[met_t.c == 0] <- NA
#need to remove the samples with NA (do you want to remove the ones with 50-50 of NAs?
#m.sr <- met_s.c[-which(rowMeans(is.na(met_s.c)) > 0.5), ]

#Normalisation
m_norm <- as.matrix(met_s.c)
m_norm <- normalize.quantiles(m_norm)
colnames(m_norm) <- colnames(met_s.c)
boxplot(m_norm, ylim = c(0, 5e+6))

data <- file.path(getwd(), "met.csv")
design <- file.path(getwd(), "design matrix.csv")

normalyzer(jobName = "met_norm", designPath = design, dataPath = data, outputDir = getwd())

m_norm <- as.matrix(met_s.c)
library(car)
m_norm_YJ <- powerTransform(m_norm, family = 'yjPower')


#Transformation
#m_sl <- log2(m.sr)
#m_tl <- log2(m.tr)
m_log2 <- log2(met_s.c +1)

#Scaling 
#Unit variance scaling (UV-scaling) = (value - mean) / stdev
#m_s <- (m_q - mean(m_q))/sd(m_q)
m_s.scale <- scale(met_s.c, center = F, scale = TRUE) #scale utilises the same formula

#Imputation using NIPALS
#creating a barplot of the explained variance of the principal components
#nipals.tune = nipals(m_s.scale, ncomp = 10)$eig
#barplot(nipals.tune, xlab = 'Principal component', ylab = 'Explained variance')
#title('Principal Components')
#imputing
#m.imp <- impute.nipals(m_s.scale, ncomp = 10)
#sum(is.na(m.imp)) #checking for any NA
#save for manual checking
#write.csv(m.sr, "metabolomics_SCx_imputed.csv")
#write.csv(m.tr, "metabolomics_Tha_imputed.csv")

#Transposing
m_transp <- t(m_s.scale)
#rownames(m_transp) <- strain
colnames(m_transp) <- rownames(met_s.c)

#Saving
write.csv(m_transp, 'met_SCx_norm_scale.csv') 
        #'metabolomics_SCx_clean.csv')
#write.csv(m_scale.tha.t, 'metabolomics_Tha_clean.csv')
