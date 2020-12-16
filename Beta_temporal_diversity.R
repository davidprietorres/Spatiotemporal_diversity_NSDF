library(raster)
library(rgdal)
library(rgeos)
library(picante)
library(betapart)
library(CommEcol)
library(geiger)
library(phytools)

#### FUNCTION: copy and paste into R (available in: http://rfunctions.blogspot.com/)####
tempbetagrid<-function(oc1, oc2, index="sorensen", phylotree=phylo, phylobeta=F){
  tempturn<-numeric(nrow(oc1))
  tempnest<-numeric (nrow(oc1))
  tempbeta<-numeric(nrow(oc1))
  tempturnbeta<-numeric (nrow(oc1))
  for(i in 1:nrow(oc1) ){
    namesoc1<-names(oc1)[oc1[i,]==1]
    namesoc2<-names(oc2)[oc2[i,]==1]
    both<-namesoc1[namesoc1%in%namesoc2]
    bothmat<-rbind(rep(1,length(both)),rep(1,length(both)))
    colnames(bothmat)<-both
    namoc1<-namesoc1[namesoc1%in%namesoc2==FALSE]
    nam1mat<-rbind(rep(1,length(namoc1)),rep(0,length(namoc1)))
    colnames(nam1mat)<-namoc1
    namoc2<-namesoc2[namesoc2%in%namesoc1==FALSE]
    nam2mat<-rbind(rep(0,length(namoc2)),rep(1,length(namoc2)))
    colnames(nam2mat)<-namoc2
    matcomp<-cbind(bothmat,nam1mat,nam2mat)
    forprune<-t(data.frame(names(data.frame(matcomp))))
    colnames(forprune)<-forprune
    ifelse(phylobeta==T, betas<-phylo.beta.pair(matcomp, prune.sample(forprune, phylo), index.family=index), betas<-beta.pair(matcomp, index.family=index) )
    tempturn[i]<-betas[[1]]
    tempnest[i]<-betas[[2]]
    tempbeta[i]<-betas[[3]]
    tempturnbeta[i]<-betas[[1]]/betas[[3]]}
  return(data.frame(turnover=tempturn,nestedness=tempnest,beta=tempbeta,turn.beta=tempturnbeta))}

###################################################################################################
################################# TEMPORAL BETA DIVERSITY  ########################################
###################################################################################################
#1. Then, load the grid (shapefile). This is a grid of 0.25 degree lat/long of STUDY AREA
shape <- readOGR(choose.files(), "grids_bosque")

present <- read.csv("E:/project_sig/heterogenization_biota_ndf/PAM_present2.csv", row.names = 1, head=T)#Load the species occurrences for the present time. 
future <- read.csv("E:/project_sig/heterogenization_biota_ndf/PAM_2050rcp45.csv", row.names = 1, head=T)# Now load the species occurrences for the future

#2. Load the phylogeny/functional tree of these species (if you want to calculate phylogenetic beta diversity):
phylo <- read.nexus("E:/project_sig/heterogenization_biota_ndf/Consenso_mayoria.nex")##Phylogenetic tree for the studied species

functional <- read.tree("E:/project_sig/heterogenization_biota_ndf/dendro_funcional_fin")#Dendrogram for the functional traits of species


#3. Call the function and get results for taxonomic information!

#########################
###taxonomic diversity###
#########################
results <- tempbetagrid(oc1=present, oc2=future, index="sorensen")
#### GRAPH Beta taxonomic temporal####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- results[,3]
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
writeRaster(rbeta, filename="205045_turnov_tx.asc", overwrite=T, suffix='names')#name of file.

# And now the map with "turnover/beta diversity" values
# Also, create a new layer in our grid file for the rate of "turnover/total beta diversity" in each cell.
shape$turnbeta <- results[,4]
# Also, assign values to the new raster according to the rates of "turnover/beta diversity" layer in our original grid.
rturnbeta <- rasterize(shape, field="turnbeta", emptyraster)
plot(rturnbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file
writeRaster(rbeta, filename="205045_turnov_Temp_tx.asc", overwrite=T, suffix='names')#name of file.


############################
###Phylogenetic diversity###
############################
# Call the function and get results for phylogenetic information
resultsphylo <- tempbetagrid(oc1=present, oc2=future, index="sorensen", phylobeta=T, phylotree=phylo)

#### GRAPH Beta phylogenetic temporal####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- resultsphylo[,3]
# Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

# Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file
writeRaster(rbeta, filename="205045_turnov_phy.asc", overwrite=T, suffix='names')#name file

# And now the map with "turnover/beta diversity" values
# Also, create a new layer in our grid file for the rate of "turnover/total beta diversity" in each cell.
shape$turnbeta <- resultsphylo[,4]
# Also, assign values to the new raster according to the rates of "turnover/beta diversity" layer in our original grid.
rturnbeta <- rasterize(shape, field="turnbeta", emptyraster)
plot(rturnbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file
writeRaster(rbeta, filename="205045_turnov_temp_ph.asc", overwrite=T, suffix='names')#name file.


############################
###Functional diversity###
############################
# Call the function and get results for functional information
resultsfunct <- tempbetagrid(oc1=present, oc2=future, index="sorensen", phylobeta=T, phylotree=as.phylo(phylo2))

#### GRAPH Beta functional temporal####
# Create a new layer in our grid file for the mean total beta diversity.
shape$betadiv <- resultsfunct[,3]
# Now create a raster with the same extent and resolution as our previous grid
emptyraster <- raster(extent(shape))
res(emptyraster)=0.25

# Assign values to the new raster according to the beta diversity layer in our shapefile.
rbeta <- rasterize(shape, field="betadiv", emptyraster)

# Make a cool color palette:
my.colors = colorRampPalette(c("lightgreen","green", "yellow","orangered", "red"))

# Plot the map
plot(rbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file
writeRaster(rbeta, filename="205045_turnov_fun.asc", overwrite=T, suffix='names')

# And now the map wth "turnover/beta diversity" values
# Also, create a new layer in our grid file for the rate of "turnover/total beta diversity" in each cell.
shape$turnbeta <- resultsfunct[,4]
# Also, assign values to the new raster according to the rates of "turnover/beta diversity" layer in our original grid.
rturnbeta <- rasterize(shape, field="turnbeta", emptyraster)
plot(rturnbeta, col=my.colors(255), frame.plot=F, axes=F, box=F, add=F, legend.width=0.8, legend.shrink=1)
setwd("E:/project_sig/heterogenization_biota_ndf/")#Directory to save the raster file
writeRaster(rbeta, filename="205045_turnov_temp_fun.asc", overwrite=T, suffix='names')