# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## The following are utility functions to be used with the file getNICice.R


## FUNCTION to download and unzip the target NIC data
# filename is the name of the NIC file to retrieve - obtained with function getNICfilename
# returns the name of the unpacked file, including its path
nicDownload<-function(filename){
	Filelocation<-paste0("https://www.natice.noaa.gov/pub/weekly/antarctic/",filename[1])
	savename<-paste0(nicsavedir,"/",substr(filename[1],29,46))	
	download.file(Filelocation,destfile=savename,method="libcurl")
	unzip(zipfile=savename,exdir=nicsavedir)
	return(savename)
}

## FUNCTION to access the downloaded NIC data
# filename is the name of a single NIC dataset, retrieved with the function getNICfilename
# fileloc is the full path to the downloaded and unzipped NIC data file (shapefile), with the .shp file name
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
	# To get NIC data... Remove warns of certificate at NIC
	set_config( config( ssl_verifypeer = 0L ) )
	
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

## FUNCTION to plot grid points, and their nearest land points
# subdf is the data.frame obtained from a spatial_points data.frame with all points on fast ice for the NIC date chosen
# edges is a spatial_lines data.frame outlining the continent's edge
# myareas is a list with the spatial polygons of fast ice for the NIC date chosen
showPoints<-function(subdf,edges,myareas){
	
	edgesldf<-SpatialLinesDataFrame(edges$landedge, data=data.frame(ID=1))
	edgedf<-fortify(edgesldf)
	pc<-ggplot(data=edgedf, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
	pc<-pc + geom_point(data=subdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subdf,aes(x=near_x,y=near_y),color="red",size=0.1)
	
	icedf<-fortify(myareas$fast)
	subedge<-subset(edgedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
	subsubdf<-subset(subdf,coords.x1>100000 & coords.x1<1000000 & coords.x2< -500000 & coords.x2 > -1500000)
	subicedf<-subset(icedf,long>100000 & long<1000000 & lat< -500000 & lat > -1500000)
	pe<-ggplot(data=subedge, aes(x=long, y=lat)) + geom_path(aes(group=group)) + theme_bw()
	pe<-pe + geom_point(data=subsubdf,aes(x=coords.x1,y=coords.x2),color="blue",size=0.2) + geom_point(data=subsubdf,aes(x=near_x,y=near_y),color="red",size=0.1)
	pe<- pe + geom_path(data=subicedf,aes(x=long, y=lat,group=group),color="gray")
	
	print(pc)
	dev.new();print(pe)
	
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
# savename is the name of the downloaded NIC ice data, including its path
# primaryproj is the projection used by the spatial points file of grid points
# myareas is the selected fast ice data (the areas with fast ice for the date chosen) - a list of spatial polygons 
# studyarea_pointswLand is the spatial points table of grid points, already attributed with nearest land point
# buffwidth is the size (in meters) of a circular buffer to count how many 50-m edgepoints are within it - a metric of fast ice stability and abundance
# plotit indicates if plots of the continent-wide and Erebus-only grid points and their nearest land points be made
calcFasIceWidth<-function(savename,primaryproj,myareas,studyarea_pointswLand,buffwidth=20000,plotit=TRUE){
	
	#get the lines for land (landedge) and fast ice (ocenedge) edges: 
	edges<-getEdges(areas=myareas)
	
	#Subset ALL points to only those within fast ice
	fast<-myareas$fast
	subsetpoints <- studyarea_pointswLand[fast,]
	subdf<-as.data.frame(subsetpoints)
	
	if(plotit==TRUE){
		#see what we got:				************************************************ Optional
		showPoints(subdf=as.data.frame(subsetpoints),edges=edges,myareas=myareas)
	}
	
	#get fast ice edge and add points 50m along it
	oedge<-edges$oceanedge
	npts<-round(lineLength(oedge)/50)	
	icedgedf<-spsample(oedge,n=npts,type="regular")
	icedgedf<-data.table(as.data.frame(icedgedf))
	
	#now loop through the subsetpoints to find the nearest point in icedgedf
	subdf$iceedge.x1<-NA;subdf$iceedge.x2<-NA;subdf$distNearestIceEdge<-NA;subdf$fastIceAbund<-NA
	edgpts<-icedgedf[,c("x","y")]
	coordinates(edgpts)<-c("x","y")
	projection(edgpts)<-CRS(projection(subsetpoints))
	
	#time it...
	tm<-Sys.time()
	for(rr in 1:(nrow(subdf))){
		icedgedf$gptx<-subdf[rr,"coords.x1"];icedgedf$gpty<-subdf[rr,"coords.x2"]
		icedgedf[,dist:=sqrt(((x-gptx)^2)+((y-gpty)^2))]
		icedgedf<-setorder(icedgedf,dist)
		
		tdf<-data.frame(gptx=icedgedf[1,"gptx"],gpty=icedgedf[1,"gpty"])
		coordinates(tdf)<-c("gptx","gpty")
		projection(tdf)<-CRS(projection(subsetpoints))
		pntbuff<-gBuffer(tdf,width=20000)
		qq<-over(tdf,pntbuff)
		subdf[rr,"iceedge.x1"]<-icedgedf[1,"x"];subdf[rr,"iceedge.x2"]<-icedgedf[1,"y"];
		subdf[rr,"distNearestIceEdge"]<-icedgedf[1,"dist"];subdf$fastIceAbund<-sum(!is.na(qq))
	}
	proctim<-Sys.time()-tm
	print(paste("Processing time:",proctim))
	
	#ice width is the sum of dist to edge and distance to land
	subdf$fastIceWidth<-subdf$distToShore+subdf$distNearestIceEdge
	
	### Convert back to spatial points...
	FastIcePoints<-subdf
	coordinates(FastIcePoints)<-c("coords.x1","coords.x2")
	proj4string(FastIcePoints)<-primaryproj
	
	return(FastIcePoints)
}

## Function to get the data for fast ice stability, persistence, and predictability
# iceyear the ice year being sought from NIC
# iceareas the set of fast ice polygons for the year and date requested
# primaryproj is the projection used by the spatial points file of grid points
# studyarea_pointswLand is the spatial points table of grid points, already attributed with nearest land point
getFastIceSPPdata<-function(iceyear,iceareas,primaryproj,studyarea_pointswLand){
	#Subset ALL points to only those within fast ice
	targetpoints <- studyarea_pointswLand[iceareas,]
	stabilitydf<-data.frame(pointId=targetpoints$pointid)
	
	## We want some specific metrics of fast ice:
	# Three values: 
	#	stability (is there fast ice at the end of December?)
	filename<-getNICfilename(getmonth="Dec",getyear=iceyear)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,DecemberIcePresence=1)
	stabilitydf<-merge(stabilitydf,lastMdf,by="pointId",all.x=TRUE)
	stabilitydf$DecemberIcePresence<-ifelse(is.na(stabilitydf$DecemberIcePresence),0,1)
	
	#	persistence (is there fast ice throughout February? Derived metric: is persistent if it’s there the past 2 years)
	persistencedf<-data.frame(pointId=targetpoints$pointid)
	year1prior<-as.character(as.numeric(iceyear)-1);year2prior<-as.character(as.numeric(iceyear)-2)
	
	filename<-getNICfilename(getmonth="Feb",getyear=year1prior)	
	lastFeb<-NROW(filename);persistFeb1y<-nicDownload(filename[lastFeb])	
	lastM_areas<-getFastIce(fileloc=persistFeb1y,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,February1Prior=1)
	persistencedf<-merge(persistencedf,lastMdf,by="pointId",all.x=TRUE)
	persistencedf$February1Prior<-ifelse(is.na(persistencedf$February1Prior),0,1)
	
	filename<-getNICfilename(getmonth="Feb",getyear=year2prior)	
	lastFeb<-NROW(filename);persistFeb2y<-nicDownload(filename[lastFeb])	
	lastM_areas<-getFastIce(fileloc=persistFeb2y,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,February2Prior=1)
	persistencedf<-merge(persistencedf,lastMdf,by="pointId",all.x=TRUE)
	persistencedf$February2Prior<-ifelse(is.na(persistencedf$February2Prior),0,1)
	
	persistencedf$Persistence2Years<-apply(persistencedf[,2:3],1,sum)
	persistencedf<-persistencedf[,c("pointId","Persistence2Years")]
	
	# 	predictability (derived from stability: stable in the past 5 years, but use the number of stable ice years).
	#	Must be predictable by end of December
	predictabilitydf<-data.frame(pointId=targetpoints$pointid)
	year3prior<-as.character(as.numeric(iceyear)-3)
	year4prior<-as.character(as.numeric(iceyear)-4)
	year5prior<-as.character(as.numeric(iceyear)-5)
	
	#Y1
	filename<-getNICfilename(getmonth="Dec",getyear=year1prior)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,Dec1yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="pointId",all.x=TRUE)
	predictabilitydf$Dec1yP<-ifelse(is.na(predictabilitydf$Dec1yP),0,1)
	
	#Y2
	filename<-getNICfilename(getmonth="Dec",getyear=year2prior)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,Dec2yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="pointId",all.x=TRUE)
	predictabilitydf$Dec2yP<-ifelse(is.na(predictabilitydf$Dec2yP),0,1)
	
	#Y3
	filename<-getNICfilename(getmonth="Dec",getyear=year3prior)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,Dec3yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="pointId",all.x=TRUE)
	predictabilitydf$Dec3yP<-ifelse(is.na(predictabilitydf$Dec3yP),0,1)
	
	#Y4
	filename<-getNICfilename(getmonth="Dec",getyear=year4prior)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,Dec4yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="pointId",all.x=TRUE)
	predictabilitydf$Dec4yP<-ifelse(is.na(predictabilitydf$Dec4yP),0,1)
	
	#Y5
	filename<-getNICfilename(getmonth="Dec",getyear=year5prior)	
	lastDec<-NROW(filename);stabilityDec<-nicDownload(filename[lastDec])
	lastM_areas<-getFastIce(fileloc=stabilityDec,dataproj=primaryproj)
	lastMPoints <- studyarea_pointswLand[lastM_areas$fast,]
	lastMdf<-data.frame(pointId=lastMPoints$pointid,Dec5yP=1)
	predictabilitydf<-merge(predictabilitydf,lastMdf,by="pointId",all.x=TRUE)
	predictabilitydf$Dec5yP<-ifelse(is.na(predictabilitydf$Dec5yP),0,1)
	
	predictdf<-predictabilitydf
	predictdf$PredictabilityDec5Years<-apply(predictDecdf[,2:6],1,sum)
	
	
	## Prepare and return output df
	outdf<-merge(stabilitydf,persistencedf,by="pointId")
	outdf<-merge(outdf,predictdf,by="pointId")
	
	return(outdf)
}

