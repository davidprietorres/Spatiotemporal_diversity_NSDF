library(raster)
library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
library(geiger)
library(phytools)


###Betagrid Function (available in: http://rfunctions.blogspot.com/)
betagrid<-function(gridshp, comp, xfeature, yfeature, radius, phylotree, phylobeta=F, index="sorensen"){
  data<-data.frame(gridshp[xfeature],gridshp[yfeature],comp)
  mean_turnover<-numeric(length(comp[,1]))
  mean_nestedness<-numeric(length(comp[,1]))
  mean_beta<-numeric(length(comp[,1]))
  for(i in 1:length(shape[[2]])){
    adj<-select.window(xf=data[i,1], yf=data[i,2], radius, xydata=data)[,-c(1,2)]
    if(phylobeta==F){
      res<-beta.pair(adj, index.family=index)
    }else if(phylobeta==T){
      res<-phylo.beta.pair(adj, phylotree, index.family=index)
    }
    mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE)
    mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE)
    mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE)
  }
  return(data.frame(cell=row.names(comp), mean_turnover, mean_nestedness, mean_beta))
}


#########################BETA Diversity###########################
##################################################################
#1.Load the grid (shapefile). This is a grid of 0.25 degree lat/long of STUDY AREA
shape <- readOGR(choose.files(), "grids_bosque")

#2.Load the species occurrences.
occurrences <- read.csv("E:/project_sig/heterogenization_biota_ndf/PAM_2050rcp45.csv", row.names = 1, head=T)

#3.Load the phylogeny of these species (if you want to calculate phylogenetic beta diversity):
phylo <- read.nexus("E:/project_sig/heterogenization_biota_ndf/Consenso_mayoria.nex")##Phylogenetic tree for the studied species

functional <- read.tree("E:/project_sig/heterogenization_biota_ndf/dendro_funcional_fin")#Dendrogram for the functional traits of species

# We have to know which features corresponds to the longitude (x) and latitude (y). Type the following code and observe that the second (2) feature corresponds to longitude while the third (3) corresponds to latitude. We need these numbers.
names(shape)
results <- betagrid(gridshp=shape, comp=occurrences, xfeature=2, yfeature=3, radius=1, index="sorensen")####Beta taxonomic diversity.

resultsphylo <- betagrid(gridshp=shape, comp=occurrences, xfeature=2, yfeature=3, radius=1, phylotree=phylo, phylobeta=T, index="sorensen")###Beta phylogenetic diversity

results_functional <- betagrid(gridshp=shape, comp=occurrences, xfeature=2, yfeature=3, radius=1, phylotree=functional, phylobeta=T, index="sorensen")###Beta phylogenetic diversity

#### GRAPH AND SAVE RASTER FILES (TAXONOMIC) ####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- results[,4]
# Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

# Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file.
writeRaster(rbeta, filename="204045_Taxon_beta.asc", overwrite=T, suffix='names')


#### GRAPH AND SAVE RASTER FILES (PHYLOGENETIC) ####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- resultsphylo[,4]
# Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

# Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file.
writeRaster(rbeta, filename="204045_phylo_beta.asc", overwrite=T, suffix='names')


#### GRAPH AND SAVE RASTER FILES (FUNCTIONAL) ####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- results_functional[,4]
# Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

# Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file.
writeRaster(rbeta, filename="204045_phylo_beta.asc", overwrite=T, suffix='names')


