library(raster)
library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
library(geiger)
library(phytools)

##ALPHA PHYLOGENETIC AND FUNCTIONAL DIVERSITY####################
##################################################################

#1.Load the grid (shapefile). This is a grid of 0.25 degree lat/long of STUDY AREA.
shape <- readOGR(choose.files(), "grids_bosque")

#2.Load the species occurrences. 
occurrences <- read.csv("E:/project_sig/heterogenization_biota_ndf/PAM_2050rcp45.csv", row.names = 1, head=T)

#3.Load the phylogeny of these species (if you want to calculate phylogenetic beta diversity):
phylo <- read.nexus("E:/project_sig/heterogenization_biota_ndf/Consenso_mayoria.nex")##########Arbol filogenetico de las especies
functional <- read.tree("E:/project_sig/heterogenization_biota_ndf/dendro_funcional_fin")##########arbol de distancias para los atributos funcionales de las especies

#4.We have to know which features corresponds to the longitude (x) and latitude (y). Type the following code and observe that the second (2) feature corresponds to longitude while the third (3) corresponds to latitude. We need these numbers.
names(shape)
phyloalpha<-mpd(occurrences, cophenetic(as.phylo(phylo)))####Alpha Phylogenetic
functinal_alpha<-mpd(occurrences, cophenetic(as.phylo(functional)))####Alpha Functional

#### GRAPH AND SAVE Alpha Phylogenetic DIVERSITY####
#5. Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- phyloalpha
#6. Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

#7.Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file.
writeRaster(rbeta, filename="205045_phylo_alp.asc", overwrite=T, suffix='names')#name of file.


#### GRAPH AND SAVE Alpha fUNCTIONAL DIVERSITY####
#5. Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- functinal_alpha
#6. Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

#7.Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file.
writeRaster(rbeta, filename="205045_functional_alp.asc", overwrite=T, suffix='names')#name of file.

