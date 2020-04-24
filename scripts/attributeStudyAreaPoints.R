## This script will attribute 5km spaced sampling location point shapefile with preset variables.  
## Sampling locations are defined around the Antarctic coast with a buffer using ArcMap.  
## The result can be used to subset more specific study areas - for example only points within a seasonal fast ice area.

## Results of this script can be read in by createNearestLandPoint.R script to find the coordinates of the nearest point along the shores of Antarctica or islands.

# TODO: Add comment
# 
# Author: Dennis Jongsomjit (djongsomjit@pointblue.org) and Leo Salas (lsalas@pointblue.org)
##############################################################################################


########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c('rgeos','raster','rgdal','proj4','maptools','plyr','data.table')
# compare to existing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install missing packages
if(length(new.packages)>0) {install.packages(new.packages)}

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Main Function
#########################

##Function to process all attribution of points with environmental covariates
## studyAreaPoints is the name of the shapefile with the points you want attributed
## dataPath is the full path to the location of the points and your environmental data sets
## This will create a studyarea_points spatial points data frame and save it to studyarea_points.RData in your working directory
attribute_studyAreaPoints<-function(studyAreaPoints,dataPath){
	
	
	#########################
	#Set environment variables 
	#########################
	
	
	## Set git data layers folder as the data directory
	working.dir<-paste0(dataPath,"layers/")
	setwd(working.dir)
	
	#########################
	#Read in all feature classes, rasters 
	#########################
	
	# Read the feature class 5km grid sample points - This is all points around the buffered coast study_area_points.shp
	dsn<-paste0(studyAreaPoints,".shp")
	pts <- readOGR(dsn=dsn,layer=studyAreaPoints)
	n<-nrow(pts)
	primaryproj<-projection(pts)
	primaryproj<-CRS(primaryproj)
	
	
	# Read the feature class ADPE colonies points sourced from http://www.penguinmap.com/
	adpe <- readOGR(dsn="adpe.shp",layer="adpe")
	
	
	# Read the feature class EMPE colonies points sourced from http://www.penguinmap.com/
	empe <- readOGR(dsn="empe.shp",layer="empe")
	
	# Read the feature class Glacier points provided by the Polar Geospatial Center
	glaciers <- readOGR(dsn="glaciers.shp",layer="glaciers")
	
	
	# Read the feature class 300 Contour points
	load(paste0(dataPath,"cs300.RData"))
	cs300<-data.table(cs300)
	
	# Read the feature class 800 Contour points
	load(paste0(dataPath,"cs800.RData"))
	cs800<-data.table(cs800)
	
	# Read the bathy and slope rasters.  Slope was created in degrees within ArcMap using the slope tool 
	slope<-raster("slope_deg.tif")
	bathy<-raster("ibcso_v1_bathy.tif")
	mystack<-stack(slope,bathy)
	
	
	#Read mean slope and bathymetry raster - the average slope within a given 5km cell block and only accounting for values outside of the land mask polygon
	meanslope<-raster("slope_deg_avg5km.tif")
	meansbathy<-raster("bathy_avg5km.tif")
	my5kmstack<-stack(meanslope,meansbathy)
	
	
	
	#########################
	#Attribution
	#########################
	
	##### ADPE  
	## For each point, find name of nearest polygon and attribute from that
	adpedf<-ldply(.data=1:nrow(pts),.fun=function(i,pts,adpe){
				mydists<-gDistance(pts[i,], adpe, byid=TRUE)
								
				tdf<-data.frame(adpedist=as.numeric(as.character(mydists[which.min(mydists)])),
						#adpesize=as.character(adpe$Current_ab[which.min(mydists)]),	#skip this - we will use the number of penguins within a radius
						adpecol=as.character(adpe$Name[which.min(mydists)]))
				return(tdf)
			},pts=pts,adpe=adpe)
	
		
	
	##### EMPE  
	empedf<-ldply(.data=1:nrow(pts),.fun=function(i,pts,empe){
			mydists<-gDistance(pts[i,], empe, byid=TRUE)
			
			tdf<-data.frame(empedist=mydists[which.min(mydists)],
					#empesize=empe$colonysize[which.min(mydists)],		#skip this - we will use the number of penguins within a radius
					empecol=as.character(empe$name[which.min(mydists)]))
			return(tdf)
		},pts=pts,empe=empe)


	
	##### Glacier  glacierdf
	glacierdf<-ldply(.data=1:nrow(pts), .fun=function(i,glaciers){
				mydists<-gDistance(pts[i,], glaciers, byid=TRUE)
				tdf<-data.frame(glacierdist=mydists[which.min(mydists)])
				return(tdf)
			},glaciers=glaciers)
	
	 
	##### 300M CONTOUR contdist
	cont300dist<-ldply(.data=1:nrow(pts), .fun=function(i,cs300){
				plon<-coordinates(pts[1,])[,1];plat<-coordinates(pts[i,])[,2]
				ans<-min(cs300[, sqrt(((coords.x1-plon)^2)+((coords.x2-plat)^2))])
				tdf<-data.frame(cont300dist=ans)
				return(tdf)
			},cs300=cs300)
	
	
	##### 800M CONTOUR contdist
	cont800dist<-ldply(.data=1:nrow(pts), .fun=function(i,cs800){
				plon<-coordinates(pts[1,])[,1];plat<-coordinates(pts[i,])[,2]
				ans<-min(cs800[, sqrt(((coords.x1-plon)^2)+((coords.x2-plat)^2))])
				tdf<-data.frame(cont800dist=ans)
				return(tdf)
			},cs800=cs800)
	
	 
	 ####combine data frames
	 mycovars<-cbind(adpedf, empedf, shoredist, glacierdf, cont300dist, cont800dist)
	 
	 
	#########################
	#Raster Extraction
	#########################
	
	
	##### BATHYMETRY 
	##Extract slope, bathy, and mean slope within 5km cell block 
	bathyslopemean<-extract(my5kmstack,pts,sp=TRUE)
	bathyslope<-extract(mystack,bathyslopemean,sp=TRUE)
	names(bathyslope)<-c("pointid","meanslope","meanbathy","slope","bathy")
	
	
	
	#########################
	#Combine all data and save
	#########################
	
	studyarea_points<-spCbind(bathyslope,mycovars)
	
	#Saves results to an RData file
	save(studyarea_points,file=paste0(dataPath,"studyarea_points.RData"))
	
	
	 
	return(studyarea_points)
	 
	 
 }

 
 
#########################
#Run function to Attribute data
#########################

#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/data/"


studyarea_points<-attribute_studyAreaPoints(studyAreaPoints="studyarea_points",dataPath=pathToGit)


 
 ##Results can then be read in by createNearestLandPoint.R script to find the coordinates of the nearest point along the shores of Antarctica or islands