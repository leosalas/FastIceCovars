list.of.packages <- c('rgeos','raster','rgdal','proj4',"pacman")
# compare to existing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install missing packages
if(length(new.packages)>0) {install.packages(new.packages)}


pacman::p_load(rgeos,raster,rgdal,proj4) 

primaryproj<-projection(pts)
primaryproj<-CRS(primaryproj)


##########################################Read in all feature classes

# Read the feature class 5km grid sample points
dsn<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/fasticepoints5km.shp"
pts <- readOGR(dsn=dsn,layer="fasticepoints5km")

primaryproj<-projection(pts)
primaryproj<-CRS(primaryproj)



# Read the feature class coastal polygons
coastpoly <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/Coastline_high_res_polygon.shp",layer="Coastline_high_res_polygon")


# Read the feature class ADPE colonies points
adpe <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/adpe.shp",layer="adpe")


# Read the feature class EMPE colonies points
empe <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/empe.shp",layer="empe")


# Read the feature class Glacier points
glaciers <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/glaciers.shp",layer="glaciers")


 
 # Read the feature class fast ice merged polygon NOV
fastice <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/nic_fastice_nov_1011.shp",layer="nic_fastice_nov_1011")

 # Read the feature class fast ice merged polygon MARCH
fasticemar <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/nic_fastice_mar_1112.shp",layer="nic_fastice_mar_1112")

 # Read the feature class AOI polygons
aoi <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/TomnodSAR_AOIs_ALL.shp",layer="TomnodSAR_AOIs_ALL")

\\prbo.org\Data\Home\Petaluma\djongsomjit\Documents\projects\sealsfromspace\testing\icedissolve.shp


# dissolved <- readOGR(dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/icedissolve.shp",layer="icedissolve")
# borders = gDifference(
   # as(dissolved,"SpatialLines"),
   # as(gUnaryUnion(dissolved),"SpatialLines"),
   # byid=TRUE)
# plot(dissolved)
# plot(borders, col="red",lwd=2,add=TRUE)
# writeOGR(borders,dsn="//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/borders.shp","borders", driver="ESRI Shapefile")



#create empty base matrix to store area_by_colony_year
adpedf<-matrix(data=NA,nrow=n,ncol=3)
adpedf<-as.data.frame(empedf)
colnames(empedf)<-c("adpedist","adpesize","adpecol")


## For each point, find name of nearest polygon 
for (i in 1:nrow(pts)) {

#table of all distances
mydists<-gDistance(pts[i,], adpe, byid=TRUE)


adpedf[i,1]<-as.character(mydists[which.min(mydists)])
adpedf[i,2] <- as.character(adpe$Current_ab[which.min(mydists)])
adpedf[i,3]<-as.character(adpe$Name[which.min(mydists)])


}




##### EMPE  #need a table of each sampling location name, and distance to nearest colony based on MAPPPD data


##  Need this for ADPE!!!
#create empty base matrix to store area_by_colony_year
empedf<-matrix(data=NA,nrow=n,ncol=3)
empedf<-as.data.frame(empedf)
colnames(empedf)<-c("empedist","empesize","empecol")


## For each point, find name of nearest polygon 
for (i in 1:nrow(pts)) {

#table of all distances
mydists<-gDistance(pts[i,], empe, byid=TRUE)


empedf[i,1]<-mydists[which.min(mydists)]
empedf[i,2] <- empe$colonysize[which.min(mydists)]
empedf[i,3]<-as.character(empe$name[which.min(mydists)])

}




##### Glacier

#create empty base df to store dist
glacierdf<-as.data.frame(rep(999.999,n))
colnames(glacierdf)<-c("glacierdist")


ptm <- proc.time()

   
## For each point, find name of nearest polygon 
for (i in 1:nrow(pts)) {

#table of all distances
mydists<-gDistance(pts[i,], glaciers, byid=TRUE)


glacierdf[i,1]<-mydists[which.min(mydists)]


}


print("total time")
 proc.time() - ptm

 
 
 

##### SHORE

#create empty base df to store dist
shoredist<-as.data.frame(rep(999.999,n))
colnames(shoredist)<-c("shoredist")


ptm <- proc.time()
  
## For each point, find name of nearest polygon 
for (i in 1:n) {

#table of all distances
mydists<-gDistance(pts[i,], coastpoly, byid=TRUE)


shoredist[i,1]<-mydists[which.min(mydists)]
#shoredf[i,2] <- as.character(coastpoly$surface[which.min(mydists)])



}

print("total time")
 proc.time() - ptm
 

 ###Attribute points with fastice dates underlying and AOI
 
 #drop all columns except year/day and rename for NOVEMBER
fastice@data <- fastice@data %>% select(19)
fastice1<-subset(fastice, yrmnth %in% 1001)
names(fastice1)<-"Nov2010a"
fastice2<-subset(fastice, yrmnth %in% 1015)
names(fastice2)<-"Nov2010b"
fastice3<-subset(fastice, yrmnth %in% 1114)
names(fastice3)<-"Nov2011a"
fastice4<-subset(fastice, yrmnth %in% 1128)
names(fastice4)<-"Nov2011b"

 
 #attribute with fast ice by date 
f1<-over(pts,fastice1,returnList=FALSE)
f2<-over(pts,fastice2,returnList=FALSE)
f3<-over(pts,fastice3,returnList=FALSE)
f4<-over(pts,fastice4,returnList=FALSE)




#drop all columns except year/day for MARCH
fasticemar@data <- fasticemar@data %>% select(c(19,20,21,22))

 
 #attribute with fast ice by date
f5<-over(pts,fasticemar,returnList=FALSE)






#reduce aoi poly columns and rename
aoi@data <- aoi@data %>% select(2)
names(aoi)<-"aoi"
 pts@data$aoi<-NA
#attribute with aoi
 pts@data$aoi<-over(pts,aoi)
 
 

 
 
 
 