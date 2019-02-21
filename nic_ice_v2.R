# TODO: Add comment
# 
# Author: djongsomjit & lsalas
###############################################################################


########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","ggplot2")

# see if packages are missing and install them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Set environment variables 
#########################

#location where downloaded NIC ice layers will be saved
nicsavedir<-"c:/temp"
#location where analysis results will be saved
resultsdir<-"c:/temp/"
#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/"

#########################
# Define functions 
#########################

## FUNCTION download and unzip NIC data
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

## FUNCTION to retrieve the NIC data filenames, one of which would be used to generate covariates
# getmonth is the name of a month, first three characters, first inupper case; defaults to "Nov"
# getyear is the 4-number year; defaults to 2011
getNICfilename<-function(getmonth="Nov",getyear=2011){
	#we send request here
	urlv<-"https://www.natice.noaa.gov/products/weekly_products.html"
	
	#make sure to ignore security cert error
	set_config( config( ssl_verifypeer = 0L ) )
	
	dd<-ifelse(getmonth=="Feb","28",
			ifelse(getmonth %in% c("Apr","Jun","Sep","Nov"),"30","31"))
	#this is what we want
	body <- list(			
			area="Antarctic",
			day0="01",
			day1=dd,
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

## FUNCTION to get edges from fast ice areas
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

## FUNCTION to make segments in 8 directions from each grid point
# pt is the point to use
# dst is the length in m for the segment
makeSegments<-function(pnt,dst,rosprj){
	px<-as.integer(pnt$coords.x1);py<-as.integer(pnt$coords.x2)
	seqx<-round(c(cos(0),cos(pi/4),0,cos(3*pi/4),cos(pi),cos(3*pi/4),0,cos(pi/4)) *dst)
	seqy<-round(c(0,sin(pi/4),sin(pi/2),sin(3*pi/4),0,sin(5*pi/4),sin(3*pi/2),sin(7*pi/4)) *dst)
	linlst<-list()
	for(ss in c(1:8)){
		ll<-cbind(c(px,px+seqx[ss]), c(py,py+seqy[ss]))
		Sl <- Line(ll)
		Sll <- Lines(list(Sl), ID = letters[ss])
		linlst[[ss]]<-Sll
	}
	ptlines<-SpatialLines(linlst)
	projection(ptlines)<-CRS(rosprj)
	return(ptlines)
}

##FUNCTION to intersect segments with ice polygons and determine shortest distance to edge, fast ice width, total ice distance
# ptcoord are the xy coordinates of the grid point
# gplines is the spatial lines data frame with the 8 line segments
# fip is the fast ice polygons data frame
# neland are the xy coordinates of the nearland point
# edgl is the spatial lines for the land's edge
getNearestIceEdge<-function(ptcoord,gplines,fip,neland,edgl){
	xld<-as.numeric(neland[1]);yld<-as.numeric(neland[2])
	xgp<-as.numeric(ptcoord[1]);ygp<-as.numeric(ptcoord[2])
	distdf<-data.frame()
	for(dd in 1:8){
		tst<-gIntersection(gplines[dd],edgl)
		if(!is.null(tst)){
			nlp<-nrow(coordinates(tst))
			landint<-TRUE
		}else{
			nlp=0;landint<-FALSE
		}
		aa<-gIntersection(gplines[dd],fip)
		outpt<-getIntersectPoints(aa,xgp,ygp)
		xout<-as.numeric(outpt[1]);yout<-as.numeric(outpt[2])
		dld<-round(sqrt(((xld-xout)^2)+((yld-yout)^2)))
		dgp<-round(sqrt(((xgp-xout)^2)+((ygp-yout)^2)))
		tdf<-data.frame(dir=dd,xout=xout,yout=yout,distland=dld,distedge=dgp,nlp=nlp,landint=landint)
		distdf<-rbind(distdf,tdf)
	}
	distdf$diffdist<-distdf$distland-distdf$distedge
	
	#retrieve the edgepoint with smallest positive diffdist that either has no landint, or even-numbered landint
	totID<-sum(distdf$distedge)
	if(FALSE %in% landint){
		distdf<-subset(distdf,landint==FALSE & diffdist>0)
	}else{
		distdf<-subset(distdf, (nlp %% 2 == 0) & diffdist>0)
	}
	if(nrow(distdf)>0){
		distdf<-distdf[order(distdf$distedge,decreasing=FALSE),]
	}else{
		distdf$xout<-NA;distdf$yout<-NA;distdf$distedge<-NA
	}
	
	resdf<-data.frame(edge_x=distdf[1,"xout"],edge_y=distdf[1,"yout"],distEdge=distdf[1,"distedge"],totalIceDist=totID)

	return(resdf)
}

## FUNCTION to extract the coordinates from the intersection points of the lines and the fast ice polygon
# aa is the intersection object
# xgp and ygp are the coordinates of the grid point
getIntersectPoints<-function(aa,xgp,ygp){
	intdf<-data.frame()
	for(ii in 1:length(aa)){
		intc<-coordinates(aa)[[ii]][[1]][2,]
		intct<-data.frame(intx=intc[1],inty=intc[2])
		intdf<-rbind(intdf,intct)
	}
	if(nrow(intdf)>1){
		intdf$disttmp<-sqrt(((intdf$intx-xgp)^2)+((intdf$inty-ygp)^2))
		intdf<-intdf[order(intdf$disttmp),]
	}
	return(intdf[1,c("intx","inty")])
}

## FUNCTION to get fast ice width for each grid point
# gp is the grid points data frame (not spatial)
# fip is the fast ice spatial polygons data frame
# edgp is the spatial lines data on land edges
# dist is the distance away in each direction to extend from the point and intersect with the ice polygon
calcFasIceWidth<-function(gp,fip,edgp,dist=100000){
	#for each point, set the 4 lines at distance 100 km
	fastdf<-ldply(.data=c(1:nrow(gp)),.fun=function(pp,gpnt,distnc,fipol,edgpt){
				spdf<-gpnt[pp,]; lid<-spdf$nearLineId
				#create the 8 segments
				gplines<-makeSegments(pnt=spdf,dst=distnc,rosprj=projection(fipol));
				#overlap with ice polygons and return closest distance furthest from nearland
				edglin<-edgpt@lines[[1]]@Lines[lid];edglins<-Lines(edglin,ID=lid);edglinslst<-list(edglins)
				edgl<-SpatialLines(edglinslst);	projection(edgl)<-projection(fipol)
				icepoint<-getNearestIceEdge(ptcoord=spdf[,c("coords.x1","coords.x2")],gplines=gplines,fip=fipol,neland=spdf[,c("near_x","near_y")],edgl=edgl);
				icepoint$iceWidth<-icepoint$distEdge+spdf$distToShore;
				return(icepoint)
			},gpnt=gp,distnc=dist,fipol=fip,edgpt=edgp)
	
	res<-cbind(gp,fastdf)
	return(res)
}

#########################
# Load data 
#########################

# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
load(paste0(pathToGit,"studyarea_points_wNearLand.RData"))
#Get main projection that will be used 
primaryproj<-CRS(projection(studyarea_pointswLand))

#### Hook for vector loop
# To get NIC data... Remove warns of certificate at NIC
set_config( config( ssl_verifypeer = 0L ) )
# List the desired NIC filenames:
filename<-getNICfilename()	#using defaults for example,Nov 2011

#Downloading one file using the first date listed, but alter as needed.
Filelocation<-paste0("https://www.natice.noaa.gov/pub/weekly/antarctic/",filename[1])
savename<-paste0(nicsavedir,"/",substr(filename[1],29,46))	
download.file(Filelocation,destfile=savename,method="libcurl")
unzip(zipfile=savename,exdir=nicsavedir)	


#########################
# Process the data
#########################

#run function to get the fast ice data ready to be processed
myareas<-getFastIce(fileloc=savename,dataproj=primaryproj)

#get the lines for land (landedge) and fast ice (ocenedge) edges: 
edges<-getEdges(areas=myareas)

#Subset ALL points to only those within fast ice
fast<-myareas$fast
subsetpoints <- studyarea_pointswLand[fast,]

#see what we got:
subdf<-as.data.frame(subsetpoints)
edgesldf<-SpatialLinesDataFrame(edges$landedge, data=data.frame(ID=1))
edgedf<-fortify(edgesldf)
p<-ggplot(data=edgedf, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
p<-p + geom_point(data=subdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subdf,aes(x=near_x,y=near_y),color="red",size=0.1)

#zoom to near Erebus bay to check - add fast ice polygon
icedf<-fortify(myareas$fast)

subedge<-subset(edgedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
subsubdf<-subset(subdf,coords.x1>100000 & coords.x1<1000000 & coords.x2< -500000 & coords.x2 > -1500000)
subicedf<-subset(icedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
p<-ggplot(data=subedge, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
p<-p + geom_point(data=subsubdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subsubdf,aes(x=near_x,y=near_y),color="red",size=0.1)
p<- p + geom_path(data=subicedf,aes(x=long, y=lat,group=group),color="gray")


###############################
#Calculate fast ice width
##############################
edgp<-edges$landedge
tm<-Sys.time()
subsetpoints_wIceWidth<-calcFasIceWidth(gp=subdf[1:10,],fip=myareas$fast,edgp=edgp,dist=100000)
Sys.time()-tm


subedge<-subset(edgedf,long< -100000 & long> -250000 & lat< 2500000 & lat > 2000000)
subsubdf<-subset(subdf,coords.x1>-250000 & coords.x1< -100000 & coords.x2< 2500000 & coords.x2 > 2000000)
subicedf<-subset(icedf,long< -100000 & long> -250000 & lat< 2500000 & lat > 2000000)
p<-ggplot(data=subedge, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
p<-p + geom_point(data=subsubdf,aes(x=coords.x1,y=coords.x2),color="blue",size=1) + geom_point(data=subsubdf,aes(x=near_x,y=near_y),color="red",size=1)
p<- p + geom_path(data=subicedf,aes(x=long, y=lat,group=group),color="gray")
p<- p + geom_point(x=subdf[1,"coords.x1"],y=subdf[1,"coords.x2"],color="green",size=2)

p<- p + geom_point(data=subsetpoints_wIceWidth[1,],aes(x=edge_x,y=edge_y),color="black",size=2) +
		geom_point(data=subsetpoints_wIceWidth[1,],aes(x=near_x,y=near_y),color="red",size=2)


