# TODO: Add comment
# 
# Author: lsalas
###############################################################################


## This script generates a data frame with all the points in the reference grid identified by their pointId, 
## and the coordinates of the nearest point along the shores of Antarctica or islands

########################
# Dependencies
########################
# Load packages
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","data.table","SDraw")

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Set environment variables 
#########################

#location where the resulting data file will be saved
resultsdir<-"/home/lsalas/nicICE/"

#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/"

#FUNCTION Constructing delta to calc roots
delta<-function(A,B,C){
	delval<-((B)^2)-(4*A*C)
	return(delval)
}

#FUNCTION to calculate the quadratic roots
getQuadRoots<-function(x,cdf,resmx,dist){
	A<-cdf[x,"qA"];B<-cdf[x,"qB"];C<-cdf[x,"qC"];m<-cdf[x,"m"];b<-cdf[x,"b"]
	px<-cdf[x,"coords.x1"];py<-cdf[x,"coords.x2"];nrvx<-cdf[x,"near_x"];nrvy<-cdf[x,"near_y"]
	#first handle separately cases with infinite slopes, where near_x and coords.x1 are equal
	#or where m is near 0
	if(m!=Inf && m!=-Inf && all.equal(px,nrvx)==FALSE && m>0.01 && mm<100){	#for flat slopes or very near 90deg, this is a 1/100th of 1% error (50m in 500KM)
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
		else if(delta(A,B,C) == 0){ # second case D=0
			rtx<--B/(2*a);rty<-(m*rtx)+b
			
		}else {rtx=rty=NA} # third case D<0
		resmx[x,1]<-rtx;resmx[x,2]<-rty
	}else if(m<=0.01){ #flat on x, add/take on X
		resmx[x,2]<-nrvy
		resmx[x,1]<-ifelse(nrvx>px,nrvx-dist,nrvx+dist)
	}else if(m>=100){ #slope close to 90 degrees, flat on y, add/take on Y
		resmx[x,1]<-nrvx
		resmx[x,2]<-ifelse(nrvy>py,nrvy-dist,nrvy+dist)
	}else if(m==Inf){ #m is Inf, meaning same x and nrvy<py, so add on Y??		#py-nrvy
		resmx[x,1]<-nrvx;resmx[x,2]<-nrvy+dist
	}else{	#m is -Inf, meaning same x and nrvy>py, so take on Y??
		resmx[x,1]<-nrvx;resmx[x,2]<-nrvy-dist
	}
	return(resmx)
	#tdf<-cdf[x,]
	#tdf$rtx<-rtx;tdf$rty<-rty
	#return(tdf)
}


####################################################
## Process all points to find nearest land edge point

## Load the data on grid points
# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
load(paste0(resultsdir,"studyarea_points.RData"))
#load(paste0(pathToGit,"studyarea_points.RData"))

#Get main projection that will be used 
primaryproj<-CRS(projection(studyarea_points))

## Load the data on land edges
#load(paste0(resultsdir,"landEdge.RData"))
load(paste0(pathToGit,"landEdge.RData"))


### NEW from here to...
npts<-round(lineLength(ledge)/50)

shoredf<-spsample(ledge,n=npts,type="regular")
shoredf<-data.table(as.data.frame(shoredf))

studyarea_pointsdf<-as.data.frame(studyarea_points)

#now loop through the studyarea_points to find the nearest point in shoredf
#how? Add each point to shoredf, calc distance, sort, take top and add to pre-dimensioned df.
studyarea_pointsdf$near.x1<-NA;studyarea_pointsdf$near.x2<-NA;studyarea_pointsdf$nearestShoreDist<-NA

#time it...
tm<-Sys.time()
for(rr in 1:(nrow(studyarea_pointsdf))){
	shoredf$gptx<-studyarea_pointsdf[rr,"coords.x1"];shoredf$gpty<-studyarea_pointsdf[rr,"coords.x2"]
	shoredf[,dist:=sqrt(((x-gptx)^2)+((y-gpty)^2))]
	shoredf<-setorder(shoredf,dist)
	studyarea_pointsdf[rr,"near.x1"]<-shoredf[1,"x"];studyarea_pointsdf[rr,"near.x2"]<-shoredf[1,"y"];studyarea_pointsdf[rr,"nearestShoreDist"]<-shoredf[1,"dist"]
}
Sys.time()-tm

### Just in case, save up to this point:
save(studyarea_pointsdf,file=paste0(resultsdir,"studyarea_points_wNearLand_v2.RData"))

### Convert back to spatial points...
studyarea_pointswLand<-studyarea_pointsdf
coordinates(studyarea_pointswLand)<-c("coords.x1","coords.x2")
proj4string(studyarea_pointswLand)<-primaryproj

save(studyarea_pointsdf,studyarea_pointswLand,file=paste0(resultsdir,"studyarea_points_wNearLand_v2.RData"))

### NEW till here

#################################
## IGNORE code below - it's the one using the geospatial intersections.



#construct data.frame of vertices
pcoordf<-ldply(.data=1:length(ledge@lines[[1]]@Lines),.fun=function(i,ledge){
			line1 <- ledge[1,]@lines[[1]]@Lines[[i]];
			line1coords <- line1@coords;
			line1coords <- as.data.frame(line1coords)
			line1coords$lineId<-i
			return(line1coords)
		},ledge=ledge)

pcoordf$edge_vertId<-1:nrow(pcoordf)

#always use a long distance for the point to draw a line and intersect with the edge of the fast ice
dist<-500000	#500 km
studyarea_pointsdf<-as.data.frame(studyarea_points)
npts<-nrow(studyarea_pointsdf)
near_lineId<-rep(9999,times=npts)
near_mx<-matrix(rep(-9999999.99999,times=npts*2),ncol=2)
#time it...
tm<-Sys.time()
for(pp in 1:npts){ #npts
	pco<-studyarea_points[pp,]@coords;
	pcoordf$dist<-apply(pcoordf[,c("x","y")],1,FUN=function(x,pco){
				ed=sqrt(((x[1]-pco[1])^2)+((x[2]-pco[2])^2));
				return(ed)},pco=pco);
	pcoordsort<-pcoordf[order(pcoordf$dist),];
	topLineId<-pcoordsort[1,"lineId"];
	qq<-as.matrix(subset(pcoordf,lineId==topLineId,select=c("x","y")));
	nearland<-nearestPointOnLine(qq, pco);
	near_lineId[pp]<-topLineId
	near_mx[pp,]<-c(nearland[1],nearland[2])
	print(pp)
}
Sys.time()-tm

if(NROW(subset(near_lineId,near_lineId<9999))==nrow(studyarea_pointsdf)){
	print("adding near line ID")
	studyarea_pointsdf$nearLineId<-near_lineId
}

near_mx<-as.data.frame(near_mx);names(near_mx)<-c("near_x","near_y")
if(nrow(subset(near_mx,near_x > -9999990))==nrow(studyarea_pointsdf)){
	print("adding nearest land point data")
	studyarea_pointsdf<-cbind(studyarea_pointsdf,near_mx)
}

### Just in case, save up to this point:
save(studyarea_pointsdf,file=paste0(resultsdir,"studyarea_points_wNearLand.RData"))

#######

## Then do the calculations for distance and 300km line...
tm<-Sys.time()
studyarea_pointsdf$distToShore<-sqrt(((studyarea_pointsdf$coords.x1 - studyarea_pointsdf$near_x)^2)+((studyarea_pointsdf$coords.x2 - studyarea_pointsdf$near_y)^2))

coordinates(studyarea_pointsdf)<-c("coords.x1","coords.x2")
proj4string(studyarea_pointsdf)<-primaryproj

studyarea_pointswLand<-studyarea_pointsdf

Sys.time()-tm

save(studyarea_pointswLand,file=paste0(resultsdir,"studyarea_points_wNearLand.RData"))

##This approach did not work
#studyarea_pointsdf$m<-(studyarea_pointsdf$coords.x2-studyarea_pointsdf$near_y)/(studyarea_pointsdf$coords.x1-studyarea_pointsdf$near_x)
#studyarea_pointsdf$b<-((studyarea_pointsdf$near_y*studyarea_pointsdf$coords.x1)-(studyarea_pointsdf$coords.x2*studyarea_pointsdf$near_x))/(studyarea_pointsdf$coords.x1-studyarea_pointsdf$near_x);
#studyarea_pointsdf$Z<-studyarea_pointsdf$b-studyarea_pointsdf$near_y
#studyarea_pointsdf$qA<-((studyarea_pointsdf$m)^2)+1
#studyarea_pointsdf$qB<-2*((studyarea_pointsdf$m*studyarea_pointsdf$Z)-studyarea_pointsdf$near_x)
#studyarea_pointsdf$qC<-((studyarea_pointsdf$near_x)^2)+((studyarea_pointsdf$Z)^2)-((dist)^2)

#rtmx<-matrix(rep(-9999999.99999,times=npts*2),ncol=2)
#for(x in 1:(npts)){
#	rtmx<-getQuadRoots(x,cdf=studyarea_pointsdf,resmx=rtmx,dist=dist)
#}
#rootdf<-ldply(.data=1:nrow(gridpoints),.fun=getQuadRoots,cdf=gridpoints)

#studyarea_pointsdf$far_x<-rtmx[,1];studyarea_pointsdf$far_y<-rtmx[,2]

### Next: create batch version of code to attribute studyareadf with some ice date data
## include description of fast ice metrics: number of intersects, fast ice width, 
## some metric of the complexity of the hosting ice polygon? some metrics of complexity of overall ice dataset?





