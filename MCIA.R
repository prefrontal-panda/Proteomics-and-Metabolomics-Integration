#MCIA Integration

setwd("C:/Users/dcho0009/Desktop/multiomics_trial")

#library
library(omicade4)
#library(ggrepel)

#Read in files
prot <- read.csv('prot_SCx_log_scale_OutReplaceRow.csv', stringsAsFactors = F)
rownames(prot) <- prot[,1]
prot <- prot[,-1]
prot <- as.data.frame(t(prot))

met <- read.csv('met_SCx_log_OutReplaceRow.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
rownames(met) <- met[,1]
met <- met[,-1]
met <- as.data.frame(t(met))

#Hierarchical Clustering
tree <- function(x){
  d <- dist(t(x))
  cluster <- hclust(d)
  dend <- as.dendrogram(cluster)
  plot(dend)
}

tree(prot)
title("Proteome Clustering")

tree(met)
title("Metabolome Clustering")

#MCIA Data Exploration
#Making a df.list object
samp.ls <- list("Proteins" = prot, "Metabolites" = met)

#Performing MCIA analysis
sam_mic <- mcia(samp.ls, cia.nf = 10)
class(sam_mic)

strain <- colnames(samp.ls$Proteins) #phenotype factor

plot(sam_mic, axes = 1:2, phenovec = strain, sample.lab = T, df.color = 1:2, gene.nlab = 5)

#Can we plot the individual windows by themselves?



#Finding molecules associated with a sample
NEC_2_mol <- selectVar(sam_mic, axis1 = c(0,1), axis2 = c(-4,-3))
NEC_2_mol

#Plotting molecules of interest
test <- plotVar(sam_mic, nlab = 10) #labelling top 10 of each mol on each axis
prot_NEC2 <- plotVar(sam_mic, var = 'A0A023IKF6', var.lab = T)
  
