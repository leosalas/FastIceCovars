# Author: Dennis Jongsomjit [djongsomjit@pointblue.org]
# Comments:
# TO DO:
##############################################################

########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","dplyr","xml2","httr")

# see if packages are missing and install them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Set environment variables 
#########################

#SET YOUR WORKING DIRECTORY
working.dir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/"
setwd(working.dir)

#location where downloaded NIC ice layers will be saved
nicsavedir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing"
#location where analysis results will be saved
resultsdir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/"
#your local git repo path
pathToGit<-"//prbo.org/data/home/petalume/djongsomjit/documents/fasticecovars/"

#########################
# Define functions 
#########################

#FUNCTION download and unzip NIC data
# filename is the name of a single NIC dataset, retrieved with the function getNICfilename
# fileloc is the full path to the downloaded NIC data file, with the .zip file name
# dataproj is the projection of the data, here named "primaryproj"
getFastIce<-function(fileloc,dataproj){
	print("processing single fast ice date")
	
	#First function for 1 shapefile
	
	##########################Read unzipped shapefile, Select subset of points from full set.
	#use this subset to attribute with Colony and Fast ice attributes
	
	#Read shapefile
	loopdsn <-paste0(str_sub(fileloc,0,-5),".shp")
	shpname<-str_sub(fileloc,-16,-5)
	#read shapefile
	region <- readOGR(dsn=loopdsn,shpname)
	
	#Project to match full points shapefile
	region<-spTransform(region,dataproj)
	
	# Define fast ice and get region layer for creating land and ocean edges
	region@data$FP <- ifelse(region@data$FP == "08",8,0)
	#get fast ice only regions
	fast<-subset(region,region@data$FP == 8)
	
	functionList <- list("fast" = fast, "region" = region)
	
	return(functionList)
	
}

#FUNCTION to retrieve the NIC data filenames, one of which would be used to generate covariates
#getmonth is the name of a month, first three characters, first inupper case; defaults to "Nov"
#getyear is the 4-number year; defaults to 2011
getNICfilename<-function(getmonth="Nov",getyear=2011){
	#we send request here
	urlv<-"https://www.natice.noaa.gov/products/weekly_products.html"
	
	#make sure to ignore security cert error
	set_config( config( ssl_verifypeer = 0L ) )
	
	#this is what we want
	body <- list(			
			area="Antarctic",
			day0="01",
			day1="30",
			format="Shapefiles",
			month0=getmonth,
			month1=getmonth,
			oldarea="Antarctic",
			oldformat="Shapefiles",
			subareas="Hemispheric",
			year0=as.character(getyear),
			year1=as.character(getyear))
	
#request the data
	r <- POST(urlv, body = body, encode = "form")
	
#check that it worked
	r$status_code==200
	
#it returns XML, and that can be a pain to parse.
#yet, we know we want something that looks like this: ...shapefiles/hemispheric/antarc999999
#so, using regular expressions to search for it...
	nv<-as.numeric(gregexpr("shapefiles/hemispheric/antarc",content(r,encoding="UTF-8"))[[1]])
	
	filenames<-character()
	#here getting the address including year
	for(nn in nv){
		filenames<-c(filenames,substr(content(r,encoding="UTF-8"),nn-5,nn+38))
	}
	
#the above include the .html and .xml files. Filter for zips only...
	filenames<-subset(filenames,grepl("^[0-9][0-9][0-9][0-9].*hemispheric.*.zip",filenames))  
}

#FUNCTION to get edges from fast ice areas
# areas is a list of shapefiles (spatialPolygons or spatialPolygonsDataFrame) which is the output of the getFastIce function 
# areas has two named elements: region and fast. See getFastIce function for details.
getEdges<-function(areas){
	region<-myareas$region
	#fast<-myareas$fast
	
	# Dissolve regions based on fast ice or not
	unionfp <- gUnaryUnion(region, id = region@data$FP)  
	
	#create line of fast ice edge away from continent
	oceanedge = gDifference(
			as(unionfp,"SpatialLines"),
			as(gUnaryUnion(unionfp),"SpatialLines"),
			byid=TRUE)
	
	#create line of fast ice edge nearest to continent
	landedge<-as(gUnaryUnion(unionfp),"SpatialLines")
	
	return(list=c(oceanedge=oceanedge,landedge=landedge))
}

#FUNCTION to get the point nearest to the land for each survey point on fast ice
getNearLand<-function(ledge,setpoints){
	#convert land edge to df of vertices for use in nearestPointOnLine
	point_coordinates = c()
	for (i in 1: length(ledge[1,]@lines[[1]]@Lines)) {
		line1 <- ledge[1,]@lines[[1]]@Lines[[i]]
		line1coords <- line1@coords
		point_coordinates = rbind(point_coordinates, line1coords)
	}
	
	#Get coordinates of nearest point on land together with sampling location and place into data frame (1 row is a xy of the nearest location and xy of the sampling location)
	near_coordinates = c()
	for (ss in 1:10){  #this will run test on 10 sampling locations
		nearland<-nearestPointOnLine(point_coordinates, setpoints[ss,]@coords)
		near_coordinates = rbind(near_coordinates, cbind(nearland[1],nearland[2],setpoints[ss,]@coords))
	}
	
	return(nearland)	
}

#########################
# Load data 
#########################
	
# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
pointdsn<-paste0(pathToGit,"studyarea_points.shp")
pts <- readOGR(dsn=pointdsn,layer="studyarea_points")
#Get main projection that will be used 
primaryproj<-CRS(projection(pts))

#### Hook for vector loop
# To get NIC data... Remove warns of certificate at NIC
set_config( config( ssl_verifypeer = 0L ) )
# List the desired NIC filenames:
filename<-getNICfilename()	#using defaults for example

#Downloading one file using the first date listed, but alter as needed.
Filelocation<-paste0("https://www.natice.noaa.gov/pub/weekly/antarctic/",filename[1])
savename<-paste0(nicsavedir,substr(filename[1],29,46))	
download.file(Filelocation,destfile=savename,method="libcurl")
unzip(zipfile=savename,exdir=nicsavedir)	


#########################
# Workhorse
#########################

#run function to get the fast ice data ready to be processed
myareas<-getFastIce(fileloc=savename,dataproj=primaryproj)

#get the lines for land (landedge) and fast ice (ocenedge) edges: 
edges<-getEdges(areas=myareas)

#Subset ALL points to only those within fast ice
fast<-myareas$fast
subsetpoints <- pts[fast,]
 
testdf<-getNearLand(ledge=edges$landedge,setpoints=subsetpoints)
   
   




 
 ###############################
 #Calculate fast ice width
 ##############################
 
 
 
 #Section to produce a dataframe which for every sampling location will provide coordinates of nearest vertices point along the landward edge, along with coordinates of sampling location.
 

  

#Next step is to draw a line between them, continue line until it intersects seaward fast edge and cut it off here. Can use gIntersection {rgeos} to find this intersection with the fast object - this takes two sp objects as arguments {sp}. Then Calculate length of line as the width of the fast ice.



  
############END Section
   
   
   
   
   
############################ 
#Attribute subset of points with distance to nearest fast ice (ocean) edge
#############################


     #create empty base df to store dist. Same length as subset of points
     fasticeedge<-as.data.frame(rep(999.999,length(subsetpoints)))
     colnames(edgedist)<-c("fasticeedge")
     
     
     ptm <- proc.time()
       
     ## For each point, find nearest distance
     for (i in 1:length(subsetpoints)) {
     
     #table of all distances
     mydists<-gDistance(subsetpoints[i,], oceanedge, byid=TRUE)
     
     
     fasticeedge[i,1]<-mydists[which.min(mydists)]



     }
 
 
 write.csv(fasticeedge)
 
 
 
 
###################################
##################################
   