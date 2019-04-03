# TODO: Add comment
# 
# Author: Dennis Jongsomjit (djongsomjit@pointblue.org) and Leo Salas (lsalas@pointblue.org)
##############################################################################################

## This script loads a data frame with all the points in the reference grid identified by their pointId, 
## and finds the coordinates of the nearest point along the shores of Antarctica or islands

########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","data.table","SDraw")

# see if packages are missing and install them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Set environment variables - EDIT AS NEEDED
#########################

#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/"

#########################
# Define functions 
#########################

makePointsWithNearestLand<-function(pathToGit){
	## Load the data on grid points
	load(paste0(pathToGit,"data/studyarea_points.RData"))
	
	#Get main projection that will be used 
	primaryproj<-CRS(projection(studyarea_points))
	
	## Load the data on land edges
	load(paste0(pathToGit,"data/landEdge.RData"))
	
	#find points every 50-m along the land edge
	npts<-round(lineLength(ledge)/50)
	shoredf<-spsample(ledge,n=npts,type="regular")
	shoredf<-data.table(as.data.frame(shoredf))
	
	#loop through the studyarea_points to find the nearest point in shoredf
	#how? Add each point to shoredf, calc distance, sort, take top and add to pre-dimensioned df
	studyarea_pointsdf<-as.data.frame(studyarea_points)
	studyarea_pointsdf$near.x1<-NA;studyarea_pointsdf$near.x2<-NA;studyarea_pointsdf$nearestShoreDist<-NA
	
	#time it...
	tm<-Sys.time()
	for(rr in 1:(nrow(studyarea_pointsdf))){
		shoredf$gptx<-studyarea_pointsdf[rr,"coords.x1"];shoredf$gpty<-studyarea_pointsdf[rr,"coords.x2"]
		shoredf[,dist:=sqrt(((x-gptx)^2)+((y-gpty)^2))]
		shoredf<-setorder(shoredf,dist)
		studyarea_pointsdf[rr,"near.x1"]<-shoredf[1,"x"];studyarea_pointsdf[rr,"near.x2"]<-shoredf[1,"y"];studyarea_pointsdf[rr,"nearestShoreDist"]<-shoredf[1,"dist"]
	}
	proctim<-Sys.time()-tm
	
	### Convert back to spatial points...
	studyarea_pointswLand<-studyarea_pointsdf
	coordinates(studyarea_pointswLand)<-c("coords.x1","coords.x2")
	proj4string(studyarea_pointswLand)<-primaryproj
	
	save(studyarea_pointsdf,studyarea_pointswLand,file=paste0(pathToGit,"data/studyarea_points_wNearLand.RData"))
	print("done making studyarea_points_wNearLand.RData")
	print(proctim)
	
}

####################################################
# give the points the nearest land point attribution
makePointsWithNearestLand(pathToGit=pathToGit)

## OUTPUT: saves the grid points file with nearest land point to the data folder in the git repository
