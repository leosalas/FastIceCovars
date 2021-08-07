# TODO: Add comment
# 
# Author: lsalas
###############################################################################


libs <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","data.table","SDraw")
lapply(libs, require, character.only = TRUE)

pathToGit<-"c:/users/lsalas/git/fasticecovars/data/"

source("c:/users/lsalas/git/fasticecovars/scripts/fastIceCovars_utils.R")

#load the problem points
load(paste0(pathToGit,"weseNoIce.RData"))

#Function to generate sampling points up to N kms away (the until parameter) from the problem point, every X meters (the by parameter)
extendGrid<-function(df,by=1000,until=5000){
	stps=seq(-1*until,until,by=by)
	edf<-ldply(1:nrow(df),function(rr,df,stps){
				pid<-df[rr,"pointid"]
				plon<-df[rr,"coords.x1"];plat<-df[rr,"coords.x2"]
				gridpts.lon<-sapply(stps,function(x,plon){x+plon},plon=plon)
				gridpts.lat<-sapply(stps,function(x,plat){x+plat},plat=plat)
				gg<-expand.grid(gridpts.lon,gridpts.lat);names(gg)<-c("coords.x1","coords.x2")
				gg$pointid<-pid
				gg$gpid<-paste0(pid,"_",1:nrow(gg))
				return(gg)
			},df=df,stps=stps)
}

attribute_noIceGrid<-function(pts,pathToGit,glaciers,cs300,cs800,my5kmstack){
	
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
	mycovars<-cbind(glacierdf, cont300dist, cont800dist)
	
	
	##Extract slope, bathy, and mean slope within 5km cell block 
	bathyslopemean<-extract(my5kmstack,pts,sp=TRUE)
	names(bathyslopemean)<-c("pointid","gpid","meanslope","meanbathy")
	
	noIceGrid_points<-spCbind(bathyslopemean,mycovars)
	
	return(noIceGrid_points)
	
	
}

makePointsWithNearestLand<-function(pts,pathToGit){
	
	#Get main projection that will be used 
	primaryproj<-CRS(projection(pts))
	
	## Load the data on land edges
	load(paste0(pathToGit,"landEdge.RData"))
	
	#find points every 50-m along the land edge
	npts<-round(lineLength(ledge)/50)
	shoredf<-spsample(ledge,n=npts,type="regular")
	shoredf<-data.table(as.data.frame(shoredf))
	
	#loop through the studyarea_points to find the nearest point in shoredf
	ptsdf<-as.data.frame(pts)
	ptsdf$near.x1<-NA;ptsdf$near.x2<-NA;ptsdf$nearestShoreDist<-NA
	for(rr in 1:(nrow(ptsdf))){
		shoredf$gptx<-ptsdf[rr,"coords.x1"];shoredf$gpty<-ptsdf[rr,"coords.x2"]
		shoredf[,dist:=sqrt(((x-gptx)^2)+((y-gpty)^2))]
		shoredf<-setorder(shoredf,dist)
		ptsdf[rr,"near.x1"]<-shoredf[1,"x"];ptsdf[rr,"near.x2"]<-shoredf[1,"y"];ptsdf[rr,"nearestShoreDist"]<-shoredf[1,"dist"]
	}
	
	### Convert back to spatial points...
	pts_wLand<-ptsdf
	coordinates(pts_wLand)<-c("coords.x1","coords.x2")
	proj4string(pts_wLand)<-primaryproj
	
	return(pts_wLand)
	
}


##############################################################
# make the grid
noIcePts<-weseNoIce[,c("pointid","coords.x1","coords.x2")]
noIceGrid<-extendGrid(df=noIcePts,by=1000,until=3000)
nrow(noIceGrid)==nrow(noIcePts)*((((3000/1000)*2)+1)^2)

###############################################################
#attribute to obtain: meanslope, meanbathy, cont300m, and cont800m
#first convert grid to spatial points
coordinates(noIceGrid)<-c("coords.x1","coords.x2")
dataproj<-CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(noIceGrid)<-dataproj

# Read the feature class Glacier points provided by the Polar Geospatial Center
glaciers <- readOGR(dsn=paste0(pathToGit,"layers"),"glaciers")

layerspth<-paste0(pathToGit,"layers/")
# Read the feature class 300 Contour points
load(paste0(pathToGit,"cs300.RData"))
cs300<-data.table(cs300)

# Read the feature class 800 Contour points
load(paste0(pathToGit,"cs800.RData"))
cs800<-data.table(cs800)

#Read mean slope and bathymetry raster - the average slope within a given 5km cell block and only accounting for values outside of the land mask polygon
meanslope<-raster(paste0(layerspth,"slope_deg_avg5km.tif"))
meansbathy<-raster(paste0(layerspth,"bathy_avg5km.tif"))
my5kmstack<-stack(meanslope,meansbathy)

attr_noIceGrid<-attribute_noIceGrid(pts=noIceGrid,pathToGit==pathToGit,glaciers=glaciers,cs300=cs300,cs800=cs800,my5kmstack=my5kmstack)


###############################################################
#attribute to obtain: nearest distance to land
attr_noIceGridwLand<-makePointsWithNearestLand(pts=attr_noIceGrid,pathToGit=pathToGit)
save(attr_noIceGridwLand,file="c:/users/lsalas/desktop/tmpSOS/attr_noIceGridwLand.RData")

###############################################################
load(file="c:/users/lsalas/desktop/tmpSOS/attr_noIceGridwLand.RData")
dataproj<-CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

#attribute to obtain: ice covariates
#### CHOOSE NIC DATE
# List the desired NIC filenames:
icemonth="Nov";iceyear=2011
filename<-getNICfilename(getmonth=icemonth,getyear=iceyear)	#   ********************* User specifies date + point buffer and the function does the rest
nicdtdf<-data.frame(NICdate=1:NROW(filename),FileName=filename);print(nicdtdf)

fn<-1 #  ******************************************** This is the date we chose

#nf<-readline(paste0("Which NICdate to use? (1 to ",NROW(filename),"): "))
#Downloading one file using the date selected, and unzipping.
if(fn<1 | fn>NROW(filename)){fn<-1}
nicsavedir<-"c:/temp"
savename<-nicDownload(filename=filename[fn],nicsavedir=nicsavedir)	

#run function to get the fast ice data ready to be processed
myareas<-getFastIce(fileloc=savename,dataproj=dataproj,nicsavedir=nicsavedir)


###############################
#Calculate fast ice width
##############################
names(attr_noIceGridwLand)<-gsub("nearestShoreDist","distToShore",names(attr_noIceGridwLand))
FastIceGridPoints<-calcFasIceWidth(savename=savename,primaryproj=dataproj,myareas=myareas,studyarea_pointswLand=attr_noIceGridwLand,buffwidth=20000,plotit=FALSE)


save(FastIceGridPoints,file=paste0(pathToGit,"FastIceGridPoints_weseNoIce_",icemonth,iceyear,".RData"))
#save as shapefile too
writeOGR(obj=FastIcePoints, dsn=paste0(resultsdir,"FastIcePoints_wFastIceCovars_",icemonth,iceyear), layer=paste0("FastIcePoints",icemonth,iceyear), driver="ESRI Shapefile")




