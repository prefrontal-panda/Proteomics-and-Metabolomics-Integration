#COMBI trial

#Unconstrained integration
#library
library(combi)

#Read in datasets
prot <- read.csv('prot_SCx_log_scale_OutReplaceRow.csv', stringsAsFactors = F)
rownames(prot) <- prot[,1]
prot <- as.matrix(prot[,-1])

met <- read.csv('met_SCx_log_OutReplaceRow.csv', stringsAsFactors = F, fileEncoding = 'UTF-8-BOM')
rownames(met) <- met[,1]
met <- as.matrix(met[,-1])

#Metadata file
meta <- read.csv("metadata.csv")
meta[,2] <- as.factor(meta[,2])

#Making list
unconst <- combi(list("proteomics" = prot, "metabolomics" = met), 
                      distribution = c("gaussian", "gaussian"), 
                      compositional = c(F,F),
                      logTransformGaussian = F)

unconst #printing information

#Plot features and samples
plot(unconst, samDF = meta, samCol = meta[,2]) #, featurePlot = 'points')

#Plotting projections
unconst_proj = plot(unconst, samDF = meta, samCol = meta$strain, returnCoords = T,
     featNum = 10)
#adding link
addLink(unconst_proj, links = cbind("ALDH2", "C9H9NO6"), Views = c(1,2),
        samples = "GAERS.6")
