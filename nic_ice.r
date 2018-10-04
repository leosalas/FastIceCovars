
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","dplyr","xml2","pacman","httr")
# compare to existing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# install missing packages
if(length(new.packages)>0) {install.packages(new.packages)}
pacman::p_load(rgdal, proj4,rgeos,maptools,raster,stringr,dplyr,xml2,pacman,httr) 

#########################
#Set variables and define functions
#########################


#SET YOUR WORKING DIRECTORY
working.dir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/"
#location where downloaded NIC ice layers will be saved
nicsavedir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/"
#location where analysis results will be saved
resultsdir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/testing/"


setwd(working.dir)

# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
pointdsn<-"V:/Project/marine/antarctica/seals_from_space/fastice/studyarea_points.shp"
pts <- readOGR(dsn=pointdsn,layer="studyarea_points")
#Get main projection that will be used 
primaryproj<-projection(pts)
primaryproj<-CRS(primaryproj)




#FUNCTION download and unzip and merge as needed
getFastIce<-function(filenams){



print("processing single fast ice date")
Filelocation<-paste0(address,filenams[1])
savename<-paste0(nicsavedir,substr(filenams[1],29,46))

set_config( config( ssl_verifypeer = 0L ) )

download.file(Filelocation,destfile=savename,method="libcurl")


unzip(savename)	

	

#First function for 1 shapefile

##########################Read unzipped shapefile, Select subset of points from full set.
#use this subset to attribute with Colony and Fast ice attributes

#Read shapefile
loopdsn <-paste0(str_sub(savename,0,-5),".shp")
shpname<-str_sub(savename,-16,-5)
#read shapefile
region <- readOGR(dsn=loopdsn,shpname)

#Project to match full points shapefile
region<-spTransform(region,primaryproj)
  

# Define fast ice and get region layer for creating land and ocean edges
region@data$FP <- ifelse(region@data$FP == "08",8,0)

#get fast ice only regions
fast<-subset(region,region@data$FP == 8)


functionList <- list("fast" = fast, "region" = region)


return(functionList)

}





########################################
#Downloading NIC ice data sets
#########################################



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

              month0="Sep",

              month1="Sep",

              oldarea="Antarctic",

              oldformat="Shapefiles",

              subareas="Hemispheric",

              year0="2013",

              year1="2013")

 

#request the data
r <- POST(urlv, body = body, encode = "form")

 

#check that it worked
r$status_code==200

 

#it returns XML, and that can be a pain to parse.
#yet, we know we want something that looks like this: ...shapefiles/hemispheric/antarc999999
#so, using regular expressions to search for it...
nv<-as.numeric(gregexpr("shapefiles/hemispheric/antarc",content(r))[[1]])

filenams<-character()

#here getting the address including year
for(nn in nv){

    filenams<-c(filenams,substr(content(r),nn-5,nn+38))

}

 

#the above include the .html and .xml files. Filter for zips only...
filenams<-subset(filenams,grepl("^[0-9][0-9][0-9][0-9].*hemispheric.*.zip",filenams))
address<-"https://www.natice.noaa.gov/pub/weekly/antarctic/"  




#run function
myareas<-getFastIce(filenams)


###########################
 #Section to create fast ice edge and landward edge lines
###########################

region<-myareas$region
fast<-myareas$fast

# Dissolve regions based on fast ice or not
unionfp <- gUnaryUnion(region, id = region@data$FP)  #fails due to topology error


#create line of fast ice edge away from continent
oceanedge = gDifference(
   as(unionfp,"SpatialLines"),
   as(gUnaryUnion(unionfp),"SpatialLines"),
   byid=TRUE)
   

#create line of fast ice edge nearest to continent
landedge<-as(gUnaryUnion(unionfp),"SpatialLines")
   
   




#Subset ALL points to only those within fast ice
subsetpoints <- pts[fast,]
 
 
   
   




 
 ###############################
 #Calculate fast ice width
 ##############################
 
 
 
 #Section to produce a dataframe which for every sampling location will provide coordinates of nearest vertices point along the landward edge, along with coordinates of sampling location.
 

  
  #convert land edge to df of vertices for use in nearestPointOnLine
 point_coordinates = c()
   for (i in 1: length(landedge[1,]@lines[[1]]@Lines)) {
  line1 <- landedge[1,]@lines[[1]]@Lines[[i]]
  line1coords <- line1@coords
  point_coordinates = rbind(point_coordinates, line1coords)
}


 #Get coordinates of nearest point on land together with sampling location and place into data frame (1 row is a xy of the nearest location and xy of the sampling location)
 near_coordinates = c()
 
#for (ss in 1:length(subsetpoints)){ #this will run all subset points !!!

for (ss in 1:10){  #this will run test on 10 sampling locations

nearland<-nearestPointOnLine(point_coordinates, subsetpoints[ss,]@coords)
near_coordinates = rbind(near_coordinates, cbind(nearland[1],nearland[2],subsetpoints[ss,]@coords))
}



#Next step is to draw a line between them, continue line until it intersects seaward fast edge and cut it off here. Can use gIntersection {rgeos} to find this intersection with the fast object - this takes two sp objects as arguments {sp}. Then Calculate length of line as the width of the fast ice.



  
############END Section
   
   
   
   
   
############################ 
#Attribute subset of points with distance to nearest fast ice (ocean) edge
#############################


     #create empty base df to store dist. Same length as subset of points
     fasticeedge<-as.data.frame(rep(999.999,length(subsetpoints))
     colnames(edgedist)<-c("fasticeedge")
     
     
     ptm <- proc.time()
       
     ## For each point, find nearest distance
     for (i in 1:length(subsetpoints)) {
     
     #table of all distances
     mydists<-gDistance(subsetpoints[i,], oceanedge, byid=TRUE)
     
     
     fasticeedge[i,1]<-mydists[which.min(mydists)]



     }
 
 
 write.csv(fasticeedge)
 
 
 }
 
###################################
##################################
   