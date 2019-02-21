# Author: Dennis Jongsomjit [djongsomjit@pointblue.org]
# Comments:
# TO DO:
##############################################################

########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr")

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
nicsavedir<-"c:/temp"
#location where analysis results will be saved
resultsdir<-"c:/temp/"
#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/"

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

#FUNCTION Constructing delta to calc roots
delta<-function(A,B,C){
	delval<-(B^2)-(4*A*C)
	return(delval)
}

#FUNCTION to calculate the quadratic roots
getQuadRoots<-function(x,cdf){
	A<-cdf[x,"qA"];B<-cdf[x,"qB"];C<-cdf[x,"qC"];m<-cdf[x,"m"];b<-cdf[x,"b"]
	px<-cdf[x,"px"];py<-cdf[x,"py"]
	if(delta(A,B,C) > 0){ # first case D>0
		x_1<-(-B+sqrt(delta(A,B,C)))/(2*A);y_1<-(m*x_1)+b 
		x_2<-(-B-sqrt(delta(A,B,C)))/(2*A);y_2<-(m*x_2)+b 
		#take the closest to px,py
		d1<-((px-x_1)^2)+((py-y_1)^2)
		d2<-((px-x_2)^2)+((py-y_2)^2)
		if(d1<d2){
			rtx<-x_1;rty<-y_1
		}else{
			rtx<-x_2;rty<-y_2
		}
	}
	else if(delta(a,b,c) == 0){ # second case D=0
		rtx<--B/(2*a);rty<-(m*rtx)+b
		
	}else {rtx=rty=NA} # third case D<0
	tdf<-cdf[x,]
	tdf$rtx<-rtx;tdf$rty<-rty
	return(tdf)
}


#FUNCTION to get the point nearest to the land for each survey point on fast ice
getNearLand<-function(ledge,setpoints,dist=100000){
	#convert land edge to df of vertices for use in nearestPointOnLine
	pcoordf<-ldply(.data=1:length(ledge@lines[[1]]@Lines),.fun=function(i,ledge){
				line1 <- ledge[1,]@lines[[1]]@Lines[[i]];
				line1coords <- line1@coords;
				line1coords <- as.data.frame(line1coords)
				line1coords$lineId<-i
				return(line1coords)
			},ledge=ledge)
	
	pcoordf$pid<-1:nrow(pcoordf)
	nrec<-nrow(pcoordf)
	near_coordinates<-ldply(.data=1:nrow(setpoints),.fun=function(ss,setpoints,pcoordf){
				pco<-setpoints[ss,]@coords;
				pcoordf$dist<-apply(pcoordf[,c("x","y")],1,FUN=function(x,pco){
							ed=sqrt(((x[1]-pco[1])^2)+((x[2]-pco[2])^2));
							return(ed)},pco=pco);
				pcoordsort<-pcoordf[order(pcoordf$dist),];
				topLineId<-pcoordsort[1,"lineId"];
				qq<-as.matrix(subset(pcoordf,lineId==topLineId,select=c("x","y")));
				nearland<-nearestPointOnLine(qq, pco);
				neardf<-data.frame(nx=nearland[1],ny=nearland[2],px=pco[1],py=pco[2],pointId=as.character(setpoints$pointid[ss]),lineId=topLineId);
				return(neardf)
			},setpoints=setpoints,pcoordf=pcoordf)
	
	#We now have two points to convert to a polyLine to intersect with the fastice edges: the grid point and its nearest coastline point
	#WE SHOULD generate the nearest coastline points for all grid points 
	#Next we find the nearest fastice edge segment, then intersect the polylines, retrieve the point, and calculate distance to nearest respective coastline point
	#ASK DJ about nearest fastice edge point first; do not re-create the wheel
	
	# Calculate locations where spatial lines intersect
	int.pts <- gIntersection(sl1, sl2, byid = TRUE)
	int.coords <- int.pts@coords
	
	
	
	#DO THIS BY HAND First
	#x2 = p, x1=n
	#convert the x coordinates to positives (longitudes are negative in the East)
	#near_coordinates$nx<--1*near_coordinates$nx;near_coordinates$px<--1*near_coordinates$px
	near_coordinates$m<-(near_coordinates$py-near_coordinates$ny)/(near_coordinates$px-near_coordinates$nx)
	near_coordinates$b<-((near_coordinates$ny*near_coordinates$px)-(near_coordinates$py*near_coordinates$nx))/(near_coordinates$px-near_coordinates$nx);
	near_coordinates$Z<-near_coordinates$b-near_coordinates$ny
	
	near_coordinates$qA<-(near_coordinates$m^2)+1
	near_coordinates$qB<-2*(((near_coordinates$m*near_coordinates$Z)-near_coordinates$nx)-near_coordinates$nx)
	near_coordinates$qC<-(near_coordinates$nx^2)+(near_coordinates$Z^2)-(dist^2)
	rootdf<-ldply(.data=1:nrow(near_coordinates),.fun=getQuadRoots,cdf=near_coordinates)
	return(rootdf)	
}

#now we have the table with all the data, plus the points that make the segment: nx,ny,px,py,rtx,rty
#need to intersect each segment with the set of segments in edge nearest to px,py
#fasticewidth is then the distance between nx,ny and edgx,edgy


#########################
# Load data 
#########################
	
# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
load(paste0(pathToGit,"studyarea_points.RData"))
#Get main projection that will be used 
primaryproj<-CRS(projection(studyarea_points))

#### Hook for vector loop
# To get NIC data... Remove warns of certificate at NIC
set_config( config( ssl_verifypeer = 0L ) )
# List the desired NIC filenames:
filename<-getNICfilename()	#using defaults for example

#Downloading one file using the first date listed, but alter as needed.
Filelocation<-paste0("https://www.natice.noaa.gov/pub/weekly/antarctic/",filename[1])
savename<-paste0(nicsavedir,"/",substr(filename[1],29,46))	
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
subsetpoints <- studyarea_points[fast,]

#debugging
ledge=edges$landedge
#hrr.shp.2 <- spTransform(hrr.shp, CRS("+init=epsg:26978"))
setpoints=subsetpoints[1:10,]
dist=100000


tm<-Sys.time()
testdf<-getNearLand(ledge=edges$landedge,setpoints=subsetpoints[1:10,]) #example using first 10 points
Sys.time()-tm

tm<-Sys.time()
testdf<-getNearLand(ledge=edges$landedge,setpoints=subsetpoints[1:100,]) #example using first 10 points
Sys.time()-tm

   




 
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
   