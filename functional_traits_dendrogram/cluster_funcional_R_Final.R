# Gower's Distance Modification (Generalization of Gower's distance) (Pavoine, 2009)
# http://rfunctions.blogspot.com/2012/07/gowers-distance-modification.html

library(raster)
library(rgdal)
library(ade4)
library(picante)
library(clue)

setwd("C:/Users/Fofis/Dropbox/Biota homogenization in NSDF/data_analysis/MARCOS")

shape <- readOGR("grids_bosque.shp")
occurrences <- read.csv("PAM_present2.csv", row.names = 1, header = T)
traits222 <- read.csv("matriz_datos_ana.csv")
traits <- read.csv("matriz_datos_ana.csv")

# Eliminate any No data in the matriz and configure the "type" for each variable used in the analyses
traits222$IUCN.category <- as.factor(traits222$IUCN.category)
traits222$Sensitivity <- as.factor(traits222$Sensitivity)
traits222$Foraging.strata <- as.factor(traits222$Foraging.strata)
traits222$Center.of.abundance <- as.factor(traits222$Center.of.abundance)
traits222$Relative.abundance <- as.factor(traits222$Relative.abundance)
traits222$Associated.to.NDF <- as.factor(traits222$Associated.to.NDF)
traits222$Number.of.ecosystems.reported <- as.factor(traits222$Number.of.ecosystems.reported)

# Nominal traits  
IUCN <- data.frame(traits222[,2])
sensitivity <- data.frame(traits222[,3])
foraging <- data.frame(traits222[,4])
c_abundance <- data.frame(traits222[,5])
r_abundance <- data.frame(traits222[,6])
ecosystems <- data.frame(traits222[,8])
associated <- data.frame(traits222[,7])
# Quantitative traits
distrib_perc <- data.frame(traits222[,9])
weight <- data.frame(log(traits222[,10]))
longitude <- data.frame(log(traits222[,11]))
sfm <- data.frame(log(traits222[,12]))
meandom <- data.frame(log(traits222[,13]))
mindom <- data.frame(log(traits222[,14]))
maxdom <- data.frame(log(traits222[,15]))
dfrange <- data.frame(log(traits222[,16]))
dfslope <- data.frame(traits222[,17])

ktab1 <- ktab.list.df(list(IUCN,sensitivity,foraging,c_abundance,r_abundance,ecosystems,associated,distrib_perc,weight,longitude,sfm,meandom,mindom,maxdom,dfrange,dfslope))

#calculate functional distance matrix
distrait <- dist.ktab(ktab1, c( "N","N","N","N","N","N","N","Q","Q", "Q", "Q", "Q","Q", "Q", "Q", "Q"), c("scaledBYrange"))

# After creating the distance matrix, we can create a functional dendrogram indicating which species are more functionally similar than others (install and load "clue" package):
disU <- ls_fit_ultrametric(distrait)

treeU <- hclust2phylog(hclust(disU, "average"))

# Transform the dendrogram into a '.phylo' file (you must install and load "picante" package):
tree<-as.phylo(treeU, use.labels=T)

# Finally, this is what we have:
taxa <- tree$tip.label
taxa <- gsub(pattern = "X", replacement = "", taxa)
taxa <- as.numeric(taxa)
taxa_labels <- colnames(occurrences)[taxa]
taxa_labels
tree$tip.label <- taxa_labels
plot(tree)

setwd("C:/Users/Fofis/Dropbox/Biota homogenization in NSDF/data_final/funcional")
write.tree(tree, file="dendro_funcional_fin")

### Finally, you can use this dendrogram to calculate (for example) functional diversity indices. Remember that Gower's method is only used to create a distance matrix. We use the UPGMA method to generate the functional dendrogram.
