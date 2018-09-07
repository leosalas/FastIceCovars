library(rgdal)
library(rgeos)
library(dplyr)
library(httr)
library(xml2)
library(proj4)
library(raster)
library(maptools) #for polygon union command (depends on rgeos)
library(stringr) #for string subs

dir.create("dissolve-example")
setwd("//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/nic_ice/")


# Read the feature class 5km grid sample points This will be the pre-attributed set of points.
pointdsn<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/projected/fasticepoints5km.shp"
pts <- readOGR(dsn=pointdsn,layer="fasticepoints5km")

primaryproj<-projection(pts)
primaryproj<-CRS(primaryproj)

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
savedir<-"//prbo.org/Data/Home/Petaluma/djongsomjit/Documents/projects/sealsfromspace/nic_ice/"

#now download and unzip
for (ff in filenams){

yy<-paste0(substr(ff,7,8))

Filelocation<-paste0(address,ff)
savename<-paste0(savedir,substr(ff,29,46))

set_config( config( ssl_verifypeer = 0L ) )

download.file(Filelocation,destfile=savename,method="libcurl")


unzip(savename)	



			

#First function for 1 shapefile

##########################Read unzipped shapefile, Select subset of points from full set.
#use this subset to attribute with Colony and Fast ice attributes

loopdsn <-paste0(str_sub(savename,0,-5),".shp")
shpname<-str_sub(savename,-16,-5)
#read shapefile
region <- readOGR(dsn=loopdsn,shpname)

#Project to match full points shapefile
region<-spTransform(region,primaryproj)
  

# Define fast ice
region@data$FP <- ifelse(region@data$FP == "08",8,0)

#get fast ice only regions
fast<-subset(region,region@data$FP == 8)

#each time through this function store results in a list


}




#Here add second function if more than 1 date.  Need a if nrow>1 then run second function

####!!! If more than 1 date, Process next ice date data and Union


#now download and unzip
mylist <- list()
for (ff in 1:length(filenams)){

##########################Read unzipped shapefiles and add to list, 



Filelocation<-paste0(address,filenams[ff])
savename<-paste0(savedir,substr(filenams[ff],29,46))

set_config( config( ssl_verifypeer = 0L ) )

download.file(Filelocation,destfile=savename,method="libcurl")


unzip(savename)	


#get name of shapefile with directory location
loopdsn <-paste0(str_sub(savename,0,-5),".shp")
shpname<-str_sub(savename,-16,-5)



#read shapefile
region <- readOGR(dsn=loopdsn,shpname)

#Project to match full points shapefile
region<-spTransform(region,primaryproj)
  


#each time through this function store results in a list
mylist[[ff]] <- region

}


mergedpolygons<-do.call(bind, mylist) 


# Define fast ice
mergedpolygons@data$FP <- ifelse(mergedpolygons@data$FP == "08",8,0)

#get fast ice only regions
fast<-subset(mergedpolygons,mergedpolygons@data$FP == 8)


unionfast<- unionSpatialPolygons(fast, ID=fast@data$FP)









#subset full points within fast ice only
subsetpoints <- pts[fast,]
 
 
 
 
 #############  Section to create fast ice edge###########################
# Dissolve regions based on fast ice or not
unionfp <- gUnaryUnion(region, id = region@data$FP)


#create line of fast ice edge away from continent
oceanedge = gDifference(
   as(unionfp,"SpatialLines"),
   as(gUnaryUnion(unionfp),"SpatialLines"),
   byid=TRUE)
   

   
   #create line of fast ice edge nearest to continent
  landedge<-as(gUnaryUnion(unionfp),"SpatialLines")
   
   
   
   
   
   
   
      ############################ Attribute subset of points with nearest fast ice ocean edge
     #create empty base df to store dist
     fasticeedge<-as.data.frame(rep(999.999,length(subsetpoints))
     colnames(edgedist)<-c("fasticeedge")
     
     
     ptm <- proc.time()
       
     ## For each point, find name of nearest polygon 
     for (i in 1:length(subsetpoints)) {
     
     #table of all distances
     mydists<-gDistance(subsetpoints[i,], oceanedge, byid=TRUE)
     
     
     fasticeedge[i,1]<-mydists[which.min(mydists)]



     }
 
 
 write.csv(fasticeedge)
 
 
 }
 
 
 #One alternative method is to use rgeos dist2line function.  This gives a matrix of distance, lon, and lat
 #or maptools nearestPointOnLine(oceanedge, subsetpoint[1,])
 
#basically, find nearest point on line to sampling point.  Draw a line between them.  continue line until it reaches landward edge and cut it off here. Calculate length of line.
   
   
   
 ############################ Attribute subset of points with ADPE and EMPE dist, size, and names
   
   
   
   
   
   
   
   #writeOGR(obj=region,dsn="Q:/Petaluma/djongsomjit/Documents/projects/sealsfromspace/nic_ice",layer="region",driver="ESRI Shapefile")