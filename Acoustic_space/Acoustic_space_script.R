###################################
# Acoustic hypervolume
# We based our script according to Blonder et al., 2018. New approaches for delineating n-dimensional hypervolumes. Methods in Ecology and Evolution 9(2): 305-319. https://doi.org/10.1111/2041-210X.12865.
###################################
# We used maximum, minimum and mean frequencies because they frequently correlate with vegetation type and structure (see text)
# Set your working directory and load databases
setwd("/home/marco/pCloudDrive/Maracucho_hv_cantos/")
db <- readxl::read_excel("data_matrix.xlsx")
song <- db[,c("Specie","meandom","mindom","maxdom")]
# Transform and scale variables
songs_log <- log(song[c("meandom","mindom","maxdom")])
songs_log <- data.frame(db["Specie"],songs_log)
songs_scale <- scale(songs_log[c("meandom","mindom","maxdom")])
songs_scale <- data.frame(db["Specie"],songs_scale)
species_list <-read.csv("species_list.csv")[-1]
songs_scale$Specie <- species_list$species_b
# Load species list for each regions with their corresponding measures
spp_regions <- read.csv("spp_regions.csv")

#############################
# install.packages("hypervolume")
# install.packages("alphahull")
library(hypervolume)
# Acoustic hypervolumes in the present
region_list_present = as.character(unique(spp_regions$region))
num_regions = length(region_list_present)
# compute hypervolumes for each species
hv_songlist_reg = new("HypervolumeList")
hv_songlist_reg@HVList = vector(mode="list",length=num_regions)
for (i in 1:num_regions)
{
  # make a hypervolume using auto-bandwidth
  hv_songlist_reg@HVList[[i]] <- hypervolume_gaussian(spp_regions[,2:4], name=as.character(region_list_present[i]),verbose=T)
}
# We highly recommend to save hypervolumes in a RDS file
saveRDS(hv_songlist_reg, file = "hv_songlist_reg_present.RDS")
######################################
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_regions, ncol=num_regions)
dimnames(overlap)=list(region_list_present, region_list_present)
for (i in 1:num_regions)
{
  for (j in i:num_regions)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_songlist_reg@HVList[[i]], hv_songlist_reg@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_songlist_reg)
# show pairwise overlaps
image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=TRUE,xlab='Region',ylab='Region', main="Acoustic overlap between regions", col=colorRampPalette(c("lightgray","red"))(100))

# Volume extraction
acou_space_present <- as.data.frame(get_volume(hv_songlist_reg))
acou_space_present <- data.frame(c(1:6),acou_space_present)
names(acou_space_present) <- c("region","volume")
kruskal.test(volume~region, data=acou_space_present) 
centroids_present <- as.data.frame(get_centroid(hv_songlist_reg))
###############################################################################################

# For 2050
region <- read.csv("regions_ID.csv")
PAM_2050rcp45 <- read.csv("PAM_2050rcp45.csv")
regions <- merge(region,PAM_2050rcp45, by.x = "ID", by.y = "CELLID", all = T)
regions$region <- as.factor(regions$region)
species <- as.data.frame(colnames(PAM_2050rcp45)[2:152])
colnames(species)="species"

#region 1
spp_PAM_2050rcp45_1 <- regions[regions$region==1,3:ncol(regions)]
spp_PAM_2050_1 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_1)) {ifelse(colSums(spp_PAM_2050rcp45_1[i])>0, spp_PAM_2050_1[i,1]<-colnames(spp_PAM_2050rcp45_1[i]),NA)
  
}
spp_PAM_2050_1 <- as.data.frame(spp_PAM_2050_1)
colnames(spp_PAM_2050_1) <- c("species","region")
spp_PAM_2050_1$region <-1
spp_PAM_2050rcp45_1 <- na.omit(spp_PAM_2050_1)

# region 2
spp_PAM_2050rcp45_2 <- regions[regions$region==2,3:ncol(regions)]
spp_PAM_2050_2 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_2)) {ifelse(colSums(spp_PAM_2050rcp45_2[i])>0, spp_PAM_2050_2[i,1]<-colnames(spp_PAM_2050rcp45_2[i]),NA)
  
}

spp_PAM_2050_2 <- as.data.frame(spp_PAM_2050_2)
colnames(spp_PAM_2050_2) <- c("species","region")
spp_PAM_2050_2$region <-2
spp_PAM_2050rcp45_2 <- na.omit(spp_PAM_2050_2)

# region 3
spp_PAM_2050rcp45_3 <- regions[regions$region==3,3:ncol(regions)]
spp_PAM_2050_3 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_3)) {ifelse(colSums(spp_PAM_2050rcp45_3[i])>0, spp_PAM_2050_3[i,1]<-colnames(spp_PAM_2050rcp45_3[i]),NA)
  
}

spp_PAM_2050_3 <- as.data.frame(spp_PAM_2050_3)
colnames(spp_PAM_2050_3) <- c("species","region")
spp_PAM_2050_3$region <-3
spp_PAM_2050rcp45_3 <- na.omit(spp_PAM_2050_3)

# region 4
spp_PAM_2050rcp45_4 <- regions[regions$region==4,3:ncol(regions)]
spp_PAM_2050_4 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_4)) {ifelse(colSums(spp_PAM_2050rcp45_4[i])>0, spp_PAM_2050_4[i,1]<-colnames(spp_PAM_2050rcp45_4[i]),NA)
  
}

spp_PAM_2050_4 <- as.data.frame(spp_PAM_2050_4)
colnames(spp_PAM_2050_4) <- c("species","region")
spp_PAM_2050_4$region <-4
spp_PAM_2050rcp45_4 <- na.omit(spp_PAM_2050_4)

# region 5
spp_PAM_2050rcp45_5 <- regions[regions$region==5,3:ncol(regions)]
spp_PAM_2050_5 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_5)) {ifelse(colSums(spp_PAM_2050rcp45_5[i])>0, spp_PAM_2050_5[i,1]<-colnames(spp_PAM_2050rcp45_5[i]),NA)
  
}

spp_PAM_2050_5 <- as.data.frame(spp_PAM_2050_5)
colnames(spp_PAM_2050_5) <- c("species","region")
spp_PAM_2050_5$region <-5
spp_PAM_2050rcp45_5 <- na.omit(spp_PAM_2050_5)

#region 6
spp_PAM_2050rcp45_6 <- regions[regions$region==6,3:ncol(regions)]
spp_PAM_2050_6 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp45_6)) {ifelse(colSums(spp_PAM_2050rcp45_6[i])>0, spp_PAM_2050_6[i,1]<-colnames(spp_PAM_2050rcp45_6[i]),NA)
  
}

spp_PAM_2050_6 <- as.data.frame(spp_PAM_2050_6)
colnames(spp_PAM_2050_6) <- c("species","region")
spp_PAM_2050_6$region <-6
spp_PAM_2050rcp45_6 <- na.omit(spp_PAM_2050_6)

#################################################################################
# Now lets make the new database with species list and acoustic variables for each region at 2050 scenario
spp_region1 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_1$species),]
spp_region1$region<-1
spp_region2 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_2$species),]
spp_region2$region<-2
spp_region3 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_3$species),]
spp_region3$region<-3
spp_region4 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_4$species),]
spp_region4$region<-4
spp_region5 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_5$species),]
spp_region5$region<-5
spp_region6 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp45_6$species),]
spp_region6$region<-6
# Join all regions in a data frame
spp_regions <- rbind(spp_region1,spp_region2,spp_region3,spp_region4,spp_region5,spp_region6)
#############################
# Acoustic hypervolumes for PAM_2050rcp45
region_list_PAM_2050rcp45 = as.character(unique(spp_regions$region))
num_regions = length(region_list_PAM_2050rcp45)
# compute hypervolumes for each species
hv_songlist_reg = new("HypervolumeList")
hv_songlist_reg@HVList = vector(mode="list",length=num_regions)
for (i in 1:num_regions)
{
  # make a hypervolume using auto-bandwidth
  hv_songlist_reg@HVList[[i]] <- hypervolume_gaussian(spp_regions[,2:4], name=as.character(region_list_PAM_2050rcp45[i]),verbose=T)
}
# Optionally you can save the file to avoid run it again and save time
saveRDS(hv_songlist_reg, file = "hv_songlist_reg_PAM_2050rcp45.RDS")
######################################
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_regions, ncol=num_regions)
dimnames(overlap)=list(region_list_PAM_2050rcp45, region_list_PAM_2050rcp45)
for (i in 1:num_regions)
{
  for (j in i:num_regions)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_songlist_reg@HVList[[i]], hv_songlist_reg@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_songlist_reg)

# show pairwise overlaps
image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=TRUE,xlab='Region',ylab='Region', main="Acoustic overlap between regions", col=colorRampPalette(c("lightgray","red"))(100))

###############################################################################
##################################
# For 2050_85 scenario
region <- read.csv("regions_ID.csv")
PAM_2050rcp85 <- read.csv("PAM_2050rcp85.csv")
regions <- merge(region,PAM_2050rcp85, by.x = "ID", by.y = "CELLID", all = T)
regions$region <- as.factor(regions$region)
species <- as.data.frame(colnames(PAM_2050rcp85)[2:152])
colnames(species)="species"

#region 1
spp_PAM_2050rcp85_1 <- regions[regions$region==1,3:ncol(regions)]
spp_PAM_2050_1 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_1)) {ifelse(colSums(spp_PAM_2050rcp85_1[i])>0, spp_PAM_2050_1[i,1]<-colnames(spp_PAM_2050rcp85_1[i]),NA)
  
}
spp_PAM_2050_1 <- as.data.frame(spp_PAM_2050_1)
colnames(spp_PAM_2050_1) <- c("species","region")
spp_PAM_2050_1$region <-1
spp_PAM_2050rcp85_1 <- na.omit(spp_PAM_2050_1)

# region 2
spp_PAM_2050rcp85_2 <- regions[regions$region==2,3:ncol(regions)]
spp_PAM_2050_2 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_2)) {ifelse(colSums(spp_PAM_2050rcp85_2[i])>0, spp_PAM_2050_2[i,1]<-colnames(spp_PAM_2050rcp85_2[i]),NA)
  
}

spp_PAM_2050_2 <- as.data.frame(spp_PAM_2050_2)
colnames(spp_PAM_2050_2) <- c("species","region")
spp_PAM_2050_2$region <-2
spp_PAM_2050rcp85_2 <- na.omit(spp_PAM_2050_2)

# region 3
spp_PAM_2050rcp85_3 <- regions[regions$region==3,3:ncol(regions)]
spp_PAM_2050_3 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_3)) {ifelse(colSums(spp_PAM_2050rcp85_3[i])>0, spp_PAM_2050_3[i,1]<-colnames(spp_PAM_2050rcp85_3[i]),NA)
  
}

spp_PAM_2050_3 <- as.data.frame(spp_PAM_2050_3)
colnames(spp_PAM_2050_3) <- c("species","region")
spp_PAM_2050_3$region <-3
spp_PAM_2050rcp85_3 <- na.omit(spp_PAM_2050_3)

# region 4
spp_PAM_2050rcp85_4 <- regions[regions$region==4,3:ncol(regions)]
spp_PAM_2050_4 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_4)) {ifelse(colSums(spp_PAM_2050rcp85_4[i])>0, spp_PAM_2050_4[i,1]<-colnames(spp_PAM_2050rcp85_4[i]),NA)
  
}

spp_PAM_2050_4 <- as.data.frame(spp_PAM_2050_4)
colnames(spp_PAM_2050_4) <- c("species","region")
spp_PAM_2050_4$region <-4
spp_PAM_2050rcp85_4 <- na.omit(spp_PAM_2050_4)

# region 5
spp_PAM_2050rcp85_5 <- regions[regions$region==5,3:ncol(regions)]
spp_PAM_2050_5 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_5)) {ifelse(colSums(spp_PAM_2050rcp85_5[i])>0, spp_PAM_2050_5[i,1]<-colnames(spp_PAM_2050rcp85_5[i]),NA)
  
}

spp_PAM_2050_5 <- as.data.frame(spp_PAM_2050_5)
colnames(spp_PAM_2050_5) <- c("species","region")
spp_PAM_2050_5$region <-5
spp_PAM_2050rcp85_5 <- na.omit(spp_PAM_2050_5)

#region 6
spp_PAM_2050rcp85_6 <- regions[regions$region==6,3:ncol(regions)]
spp_PAM_2050_6 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2050rcp85_6)) {ifelse(colSums(spp_PAM_2050rcp85_6[i])>0, spp_PAM_2050_6[i,1]<-colnames(spp_PAM_2050rcp85_6[i]),NA)
  
}

spp_PAM_2050_6 <- as.data.frame(spp_PAM_2050_6)
colnames(spp_PAM_2050_6) <- c("species","region")
spp_PAM_2050_6$region <-6
spp_PAM_2050rcp85_6 <- na.omit(spp_PAM_2050_6)

#################################################################################
# Now lets make the new database with species list and acoustic variables for each region at 2050_85 scenario
spp_region1 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_1$species),]
spp_region1$region<-1
spp_region2 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_2$species),]
spp_region2$region<-2
spp_region3 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_3$species),]
spp_region3$region<-3
spp_region4 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_4$species),]
spp_region4$region<-4
spp_region5 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_5$species),]
spp_region5$region<-5
spp_region6 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2050rcp85_6$species),]
spp_region6$region<-6
# Join all in a single data frame
spp_regions <- rbind(spp_region1,spp_region2,spp_region3,spp_region4,spp_region5,spp_region6)
#############################
# Acoustic hypervolumes for PAM_2050rcp85
region_list_PAM_2050rcp85 = as.character(unique(spp_regions$region))
num_regions = length(region_list_PAM_2050rcp85)
# compute hypervolumes for each species
hv_songlist_reg = new("HypervolumeList")
hv_songlist_reg@HVList = vector(mode="list",length=num_regions)
for (i in 1:num_regions)
{
  # make a hypervolume using auto-bandwidth
  hv_songlist_reg@HVList[[i]] <- hypervolume_gaussian(spp_regions[,2:4], name=as.character(region_list_PAM_2050rcp85[i]),verbose=T)
}

saveRDS(hv_songlist_reg, file = "hv_songlist_reg_PAM_2050rcp85.RDS")
######################################
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_regions, ncol=num_regions)
dimnames(overlap)=list(region_list_PAM_2050rcp85, region_list_PAM_2050rcp85)
for (i in 1:num_regions)
{
  for (j in i:num_regions)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_songlist_reg@HVList[[i]], hv_songlist_reg@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_songlist_reg)
# show pairwise overlaps
image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=TRUE,xlab='Region',ylab='Region', main="Acoustic overlap between regions", col=colorRampPalette(c("lightgray","red"))(100))

################################################################################################
##################################
# For scenario 2070_45
region <- read.csv("regions_ID.csv")
PAM_2070rcp45 <- read.csv("PAM_2070rcp45.csv")
regions <- merge(region,PAM_2070rcp45, by.x = "ID", by.y = "CELLID", all = T)
regions$region <- as.factor(regions$region)
species <- as.data.frame(colnames(PAM_2070rcp45)[2:152])
colnames(species)="species"

#region 1
spp_PAM_2070rcp45_1 <- regions[regions$region==1,3:ncol(regions)]
spp_PAM_2070_1 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_1)) {ifelse(colSums(spp_PAM_2070rcp45_1[i])>0, spp_PAM_2070_1[i,1]<-colnames(spp_PAM_2070rcp45_1[i]),NA)
  
}
spp_PAM_2070_1 <- as.data.frame(spp_PAM_2070_1)
colnames(spp_PAM_2070_1) <- c("species","region")
spp_PAM_2070_1$region <-1
spp_PAM_2070rcp45_1 <- na.omit(spp_PAM_2070_1)

# region 2
spp_PAM_2070rcp45_2 <- regions[regions$region==2,3:ncol(regions)]
spp_PAM_2070_2 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_2)) {ifelse(colSums(spp_PAM_2070rcp45_2[i])>0, spp_PAM_2070_2[i,1]<-colnames(spp_PAM_2070rcp45_2[i]),NA)
  
}

spp_PAM_2070_2 <- as.data.frame(spp_PAM_2070_2)
colnames(spp_PAM_2070_2) <- c("species","region")
spp_PAM_2070_2$region <-2
spp_PAM_2070rcp45_2 <- na.omit(spp_PAM_2070_2)

# region 3
spp_PAM_2070rcp45_3 <- regions[regions$region==3,3:ncol(regions)]
spp_PAM_2070_3 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_3)) {ifelse(colSums(spp_PAM_2070rcp45_3[i])>0, spp_PAM_2070_3[i,1]<-colnames(spp_PAM_2070rcp45_3[i]),NA)
  
}

spp_PAM_2070_3 <- as.data.frame(spp_PAM_2070_3)
colnames(spp_PAM_2070_3) <- c("species","region")
spp_PAM_2070_3$region <-3
spp_PAM_2070rcp45_3 <- na.omit(spp_PAM_2070_3)

# region 4
spp_PAM_2070rcp45_4 <- regions[regions$region==4,3:ncol(regions)]
spp_PAM_2070_4 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_4)) {ifelse(colSums(spp_PAM_2070rcp45_4[i])>0, spp_PAM_2070_4[i,1]<-colnames(spp_PAM_2070rcp45_4[i]),NA)
  
}

spp_PAM_2070_4 <- as.data.frame(spp_PAM_2070_4)
colnames(spp_PAM_2070_4) <- c("species","region")
spp_PAM_2070_4$region <-4
spp_PAM_2070rcp45_4 <- na.omit(spp_PAM_2070_4)

# region 5
spp_PAM_2070rcp45_5 <- regions[regions$region==5,3:ncol(regions)]
spp_PAM_2070_5 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_5)) {ifelse(colSums(spp_PAM_2070rcp45_5[i])>0, spp_PAM_2070_5[i,1]<-colnames(spp_PAM_2070rcp45_5[i]),NA)
  
}

spp_PAM_2070_5 <- as.data.frame(spp_PAM_2070_5)
colnames(spp_PAM_2070_5) <- c("species","region")
spp_PAM_2070_5$region <-5
spp_PAM_2070rcp45_5 <- na.omit(spp_PAM_2070_5)

#region 6
spp_PAM_2070rcp45_6 <- regions[regions$region==6,3:ncol(regions)]
spp_PAM_2070_6 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp45_6)) {ifelse(colSums(spp_PAM_2070rcp45_6[i])>0, spp_PAM_2070_6[i,1]<-colnames(spp_PAM_2070rcp45_6[i]),NA)
  
}

spp_PAM_2070_6 <- as.data.frame(spp_PAM_2070_6)
colnames(spp_PAM_2070_6) <- c("species","region")
spp_PAM_2070_6$region <-6
spp_PAM_2070rcp45_6 <- na.omit(spp_PAM_2070_6)

#################################################################################
# Now lets make the new database with species list and acoustic variables for each region at 2070_45 scenario
spp_region1 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_1$species),]
spp_region1$region<-1
spp_region2 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_2$species),]
spp_region2$region<-2
spp_region3 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_3$species),]
spp_region3$region<-3
spp_region4 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_4$species),]
spp_region4$region<-4
spp_region5 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_5$species),]
spp_region5$region<-5
spp_region6 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp45_6$species),]
spp_region6$region<-6
# Join all in a single data frame
spp_regions <- rbind(spp_region1,spp_region2,spp_region3,spp_region4,spp_region5,spp_region6)
#############################
# Acoustic hypervolumes for PAM_2070rcp45
region_list_PAM_2070rcp45 = as.character(unique(spp_regions$region))
num_regions = length(region_list_PAM_2070rcp45)
# compute hypervolumes for each species
hv_songlist_reg = new("HypervolumeList")
hv_songlist_reg@HVList = vector(mode="list",length=num_regions)
for (i in 1:num_regions)
{
  # make a hypervolume using auto-bandwidth
  hv_songlist_reg@HVList[[i]] <- hypervolume_gaussian(spp_regions[,2:4], name=as.character(region_list_PAM_2070rcp45[i]),verbose=T)
}

saveRDS(hv_songlist_reg, file = "hv_songlist_reg_PAM_2070rcp45.RDS")
######################################
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_regions, ncol=num_regions)
dimnames(overlap)=list(region_list_PAM_2070rcp45, region_list_PAM_2070rcp45)
for (i in 1:num_regions)
{
  for (j in i:num_regions)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_songlist_reg@HVList[[i]], hv_songlist_reg@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_songlist_reg)
# show pairwise overlaps
image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=TRUE,xlab='Region',ylab='Region', main="Acoustic overlap between regions", col=colorRampPalette(c("lightgray","red"))(100))

################################################################################################
##################################
# For 2070_45 scenario
region <- read.csv("identificador_regions.csv")
PAM_2070rcp85 <- read.csv("PAM_2070rcp85.csv")
regions <- merge(region,PAM_2070rcp85, by.x = "ID", by.y = "CELLID", all = T)
regions$region <- as.factor(regions$region)
# Voy a conservar el orden en que aparecen en PAM para hacer menos cambios
species <- as.data.frame(colnames(PAM_2070rcp85)[2:152])
colnames(species)="species"

#region 1
spp_PAM_2070rcp85_1 <- regions[regions$region==1,3:ncol(regions)]
spp_PAM_2070_1 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_1)) {ifelse(colSums(spp_PAM_2070rcp85_1[i])>0, spp_PAM_2070_1[i,1]<-colnames(spp_PAM_2070rcp85_1[i]),NA)
  
}
spp_PAM_2070_1 <- as.data.frame(spp_PAM_2070_1)
colnames(spp_PAM_2070_1) <- c("species","region")
spp_PAM_2070_1$region <-1
spp_PAM_2070rcp85_1 <- na.omit(spp_PAM_2070_1)

# region 2
spp_PAM_2070rcp85_2 <- regions[regions$region==2,3:ncol(regions)]
spp_PAM_2070_2 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_2)) {ifelse(colSums(spp_PAM_2070rcp85_2[i])>0, spp_PAM_2070_2[i,1]<-colnames(spp_PAM_2070rcp85_2[i]),NA)
  
}

spp_PAM_2070_2 <- as.data.frame(spp_PAM_2070_2)
colnames(spp_PAM_2070_2) <- c("species","region")
spp_PAM_2070_2$region <-2
spp_PAM_2070rcp85_2 <- na.omit(spp_PAM_2070_2)

# region 3
spp_PAM_2070rcp85_3 <- regions[regions$region==3,3:ncol(regions)]
spp_PAM_2070_3 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_3)) {ifelse(colSums(spp_PAM_2070rcp85_3[i])>0, spp_PAM_2070_3[i,1]<-colnames(spp_PAM_2070rcp85_3[i]),NA)
  
}

spp_PAM_2070_3 <- as.data.frame(spp_PAM_2070_3)
colnames(spp_PAM_2070_3) <- c("species","region")
spp_PAM_2070_3$region <-3
spp_PAM_2070rcp85_3 <- na.omit(spp_PAM_2070_3)

# region 4
spp_PAM_2070rcp85_4 <- regions[regions$region==4,3:ncol(regions)]
spp_PAM_2070_4 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_4)) {ifelse(colSums(spp_PAM_2070rcp85_4[i])>0, spp_PAM_2070_4[i,1]<-colnames(spp_PAM_2070rcp85_4[i]),NA)
  
}

spp_PAM_2070_4 <- as.data.frame(spp_PAM_2070_4)
colnames(spp_PAM_2070_4) <- c("species","region")
spp_PAM_2070_4$region <-4
spp_PAM_2070rcp85_4 <- na.omit(spp_PAM_2070_4)

# region 5
spp_PAM_2070rcp85_5 <- regions[regions$region==5,3:ncol(regions)]
spp_PAM_2070_5 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_5)) {ifelse(colSums(spp_PAM_2070rcp85_5[i])>0, spp_PAM_2070_5[i,1]<-colnames(spp_PAM_2070rcp85_5[i]),NA)
  
}

spp_PAM_2070_5 <- as.data.frame(spp_PAM_2070_5)
colnames(spp_PAM_2070_5) <- c("species","region")
spp_PAM_2070_5$region <-5
spp_PAM_2070rcp85_5 <- na.omit(spp_PAM_2070_5)

#region 6
spp_PAM_2070rcp85_6 <- regions[regions$region==6,3:ncol(regions)]
spp_PAM_2070_6 <- matrix(NA, 151,2)
for (i in 1:ncol(spp_PAM_2070rcp85_6)) {ifelse(colSums(spp_PAM_2070rcp85_6[i])>0, spp_PAM_2070_6[i,1]<-colnames(spp_PAM_2070rcp85_6[i]),NA)
  
}

spp_PAM_2070_6 <- as.data.frame(spp_PAM_2070_6)
colnames(spp_PAM_2070_6) <- c("species","region")
spp_PAM_2070_6$region <-6
spp_PAM_2070rcp85_6 <- na.omit(spp_PAM_2070_6)

#################################################################################
# Now lets make the new database with species list and acoustic variables for each region at 2070_45 scenario
spp_region1 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_1$species),]
spp_region1$region<-1
spp_region2 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_2$species),]
spp_region2$region<-2
spp_region3 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_3$species),]
spp_region3$region<-3
spp_region4 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_4$species),]
spp_region4$region<-4
spp_region5 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_5$species),]
spp_region5$region<-5
spp_region6 <- songs_scale[which(songs_scale$Specie%in%spp_PAM_2070rcp85_6$species),]
spp_region6$region<-6
# Join all in a single data frame
spp_regions <- rbind(spp_region1,spp_region2,spp_region3,spp_region4,spp_region5,spp_region6)
#############################
# Acoustic hypervolumes for PAM_2070rcp85
region_list_PAM_2070rcp85 = as.character(unique(spp_regions$region))
num_regions = length(region_list_PAM_2070rcp85)
# compute hypervolumes for each species
hv_songlist_reg = new("HypervolumeList")
hv_songlist_reg@HVList = vector(mode="list",length=num_regions)
for (i in 1:num_regions)
{
  # make a hypervolume using auto-bandwidth
  hv_songlist_reg@HVList[[i]] <- hypervolume_gaussian(spp_regions[,2:4], name=as.character(region_list_PAM_2070rcp85[i]),verbose=T)
}

saveRDS(hv_songlist_reg, file = "hv_songlist_reg_PAM_2070rcp85.RDS")
######################################
# compute all pairwise overlaps
overlap = matrix(NA, nrow=num_regions, ncol=num_regions)
dimnames(overlap)=list(region_list_PAM_2070rcp85, region_list_PAM_2070rcp85)
for (i in 1:num_regions)
{
  for (j in i:num_regions)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_songlist_reg@HVList[[i]], hv_songlist_reg@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_songlist_reg)
# show pairwise overlaps
image(x=1:nrow(overlap), y=1:nrow(overlap), z=overlap,axes=TRUE,xlab='Region',ylab='Region', main="Acoustic overlap between regions", col=colorRampPalette(c("lightgray","red"))(100))

################################################################################################
################################################################################################
# Acoustic space variation for each region among scenarios
# If you previously saved the RDS files set the working directory to the folder were they are located, otherwise you will need to run the full script in order to estimate all of them and omit lines 608-613
setwd("/home/marco/pCloudDrive/Maracucho_hv_cantos/")
hv_present <-readRDS("hv_songlist_reg_present.RDS")
hv_2050_rcp45 <- readRDS("hv_songlist_reg_PAM_2050rcp45.RDS")
hv_2050_rcp85 <- readRDS("hv_songlist_reg_PAM_2050rcp85.RDS")
hv_2070_rcp45 <- readRDS("hv_songlist_reg_PAM_2070rcp45.RDS")
hv_2070_rcp85 <- readRDS("hv_songlist_reg_PAM_2070rcp85.RDS")
# Join scenarios hypervolumes
hv_scenarios <- hypervolume_join(hv_present,hv_2050_rcp45,hv_2050_rcp85,hv_2070_rcp45,hv_2070_rcp85)

# For region 1
# compute all pairwise overlaps
overlap_1 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_1)=list(escenario,escenario)
hv_region_1 <- hv_scenarios[[c(1,7,13,19,25)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_1@HVList[[i]], hv_region_1@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_1 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_1[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_1)
# show pairwise overlaps
par(mar=c(4,6,3,1))
image(x=1:nrow(overlap_1), y=1:nrow(overlap_1), z=overlap_1,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 1",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_1)[[1]])),dimnames(overlap_1)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_1)[[2]])),dimnames(overlap_1)[[2]],las=1,cex.axis=1)

# For region 2
# compute all pairwise overlaps
overlap_2 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_2)=list(escenario,escenario)
hv_region_2 <- hv_scenarios[[c(2,8,14,20,26)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_2@HVList[[i]], hv_region_2@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_2 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_2[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_2)
# show pairwise overlaps
image(x=1:nrow(overlap_2), y=1:nrow(overlap_2), z=overlap_2,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 2",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_2)[[1]])),dimnames(overlap_2)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_2)[[2]])),dimnames(overlap_2)[[2]],las=1,cex.axis=1)

# For region 3
# compute all pairwise overlaps
overlap_3 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_3)=list(escenario,escenario)
hv_region_3 <- hv_scenarios[[c(3,9,15,21,27)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_2@HVList[[i]], hv_region_2@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_3 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_3[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_3)
# show pairwise overlaps
image(x=1:nrow(overlap_3), y=1:nrow(overlap_3), z=overlap_3,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 3",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_3)[[1]])),dimnames(overlap_3)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_3)[[2]])),dimnames(overlap_3)[[2]],las=1,cex.axis=1)

# For region 4
# compute all pairwise overlaps
overlap_4 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_4)=list(escenario,escenario)
hv_region_4 <- hv_scenarios[[c(4,10,16,22,28)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_2@HVList[[i]], hv_region_2@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_4 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_4[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_4)
# show pairwise overlaps
image(x=1:nrow(overlap_4), y=1:nrow(overlap_4), z=overlap_4,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 4",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_4)[[1]])),dimnames(overlap_4)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_4)[[2]])),dimnames(overlap_4)[[2]],las=1,cex.axis=1)

# For region 5
# compute all pairwise overlaps
overlap_5 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_5)=list(escenario,escenario)
hv_region_5 <- hv_scenarios[[c(5,11,17,23,29)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_2@HVList[[i]], hv_region_2@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_5 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_5[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_5)
# show pairwise overlaps
image(x=1:nrow(overlap_5), y=1:nrow(overlap_5), z=overlap_5,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 5",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_5)[[1]])),dimnames(overlap_5)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_5)[[2]])),dimnames(overlap_5)[[2]],las=1,cex.axis=1)

# For region 6
# compute all pairwise overlaps
overlap_6 = matrix(NA, nrow=5, ncol=5)
escenario <-as.list(c("present","2050_rcp45","2050_rcp85","2070_rcp45","2070_rcp_85"))
dimnames(overlap_6)=list(escenario,escenario)
hv_region_6 <- hv_scenarios[[c(6,12,18,24,30)]]
for (i in 1:5)
{
  for (j in i:5)
  {
    if (i!=j)
    {
      # compute set operations on each pair
      this_set = hypervolume_set(hv_region_2@HVList[[i]], hv_region_2@HVList[[j]], check.memory=FALSE)
      # calculate a Sorensen overlap_6 index (2 x shared volume / sum of |hv1| + |hv2|)
      overlap_6[i,j] = hypervolume_overlap_statistics(this_set)["sorensen"]
    }
  } 
}

# show all hypervolumes
plot(hv_region_6)
# show pairwise overlaps
image(x=1:nrow(overlap_6), y=1:nrow(overlap_6), z=overlap_6,axes=FALSE,xlab='Scenarios',ylab='', main="Acoustic overlap between scenarios of region 6",col=colorRampPalette(c("lightgray","red"))(100))
box()
axis(side=1, at=1:(length(dimnames(overlap_6)[[1]])),dimnames(overlap_6)[[1]],las=1,cex.axis=1)
axis(side=2, at=1:(length(dimnames(overlap_6)[[2]])),dimnames(overlap_6)[[2]],las=1,cex.axis=1)

### The more red, the more similar are hypervolumes
###############################################################################
# Obtain hypervolume centroid values and volume
# First we extract centroid of each scenario
centroids <- as.data.frame(get_centroid(hv_scenarios))
# install.packages("rgl")
# install.packages("car")
library(rgl)
library(car)
centroids$regions <- rownames(centroids)
centroids$regions <-as.factor(centroids$regions)
# 3Dplot of centroids for each region by scenario
# Colors for each region
# 1: Blue (present: blue4, dodgerblue, cyan, steelblue1, turquoise2)
# 2: Red (red, darkred,lightcoral, deeppink, pink)
# 3: Green (chartreuse4, chartreuse3, chartreuse, darkolivegreen3, seagreen1)
# 4: Yellow (gold, gold3, khaki, yellow, darkgoldenrod1)
# 5: Black (black, gray27, gray45, gray80, dimgray)
# 6: Purple (darkorchid4, darkmagenta, darkviolet, darkorchid1, plum)
colors <- c("blue4", "dodgerblue", "cyan", "steelblue1", "turquoise2", "red", "darkred", "lightcoral", "deeppink", "pink", "chartreuse4", "chartreuse3", "chartreuse", "darkolivegreen3", "seagreen1", "gold", "gold3", "khaki", "yellow", "darkgoldenrod1", "black", "gray27", "gray45", "gray80", "dimgray")
# Plot
par(mar=c(0,0,0,0))
plot3d( 
  x=centroids$meandom, y=centroids$mindom, z=centroids$maxdom, 
  col = colors, 
  type = 'p', 
  pch = c(15:19,8),
  size = 18,
  xlab="meandom", ylab="mindom", zlab="maxdom")

centroids$volume <- get_volume(hv_scenarios)
#######################################################
# Testing for statistical differences
acou_space_scen <- as.data.frame(get_volume(hv_scenarios))
acou_space_scen <- data.frame(rep.int(1:6,5),acou_space_scen)
names(acou_space_scen) <- c("region","volume")
acou_space_scen$region <- c("present 1","present 2","present 3","present 4","present 5","present 6","2050_rcp45 1","2050_rcp45 2","2050_rcp45 3","2050_rcp45 4","2050_rcp45 5","2050_rcp45 6","2050_rcp85 1","2050_rcp85 2","2050_rcp85 3","2050_rcp85 4","2050_rcp85 5","2050_rcp85 6","2070_rcp45 1","2070_rcp45 2","2070_rcp45 3","2070_rcp45 4","2070_rcp45 5","2070_rcp45 6","2070_rcp85 1","2070_rcp85 2","2070_rcp85 3","2070_rcp85 4","2070_rcp85 5","2070_rcp85 6")
kruskal.test(volume~region, data=acou_space_scen[c(1,7,13,19,25),]) 
kruskal.test(volume~region, data=acou_space_scen[c(2,8,14,20,26),]) 
kruskal.test(volume~region, data=acou_space_scen[c(3,9,15,21,27),]) 
kruskal.test(volume~region, data=acou_space_scen[c(4,10,16,22,28),]) 
kruskal.test(volume~region, data=acou_space_scen[c(5,11,17,23,29),]) 
kruskal.test(volume~region, data=acou_space_scen[c(6,12,18,24,30),]) 
# Get centroids by scenario for each region
centroids_scen <- as.data.frame(get_centroid(hv_scenarios))
