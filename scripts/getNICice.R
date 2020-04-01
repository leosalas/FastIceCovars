# TODO: Think of ways to re-write loop processing subdf with plyr, parallelization
# 
# Author: Dennis Jongsomjit (djongsomjit@pointblue.org) and Leo Salas (lsalas@pointblue.org)
##############################################################################################

## NOTE: you may want to use the jupyter notebook instead.

########################
# Dependencies
########################
# list packages required
list.of.packages <- c("rgdal", "proj4","rgeos","maptools","raster","stringr","plyr","dplyr","xml2","httr","ggplot2","data.table","SDraw")

# see if packages are missing and install them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) {install.packages(new.packages)}

# load the packages
lapply(list.of.packages, require, character.only = TRUE)


#########################
#Set environment variables - EDIT AS NEEDED
#########################

#location where downloaded NIC ice layers will be saved
nicsavedir<-"c:/temp"
#location where analysis results will be saved
resultsdir<-"c:/temp/"
#your local git repo path
pathToGit<-"c:/users/lsalas/git/fasticecovars/"

#########################
# Source functions 
#########################

source(paste0(pathToGit,"scripts/fasticeCovars_utils.R"))

#########################
# Load the NIC data 
#########################

#### CHOOSE NIC DATE
# List the desired NIC filenames:
icemonth="Nov";iceyear=2011
filename<-getNICfilename(getmonth=icemonth,getyear=iceyear)	#   ********************* User specifies date + point buffer and the function does the rest
nicdtdf<-data.frame(NICdate=1:NROW(filename),FileName=filename);print(nicdtdf)

print("******************************************************************************************************************")
print("***** STOP: review the NIC dates printed above; choose one date by entering the line number below  ***************")
print("******************************************************************************************************************")

fn<-1 #  ******************************************** User specifies the desired available NIC date

#nf<-readline(paste0("Which NICdate to use? (1 to ",NROW(filename),"): "))
#Downloading one file using the date selected, and unzipping.
if(fn<1 | fn>NROW(filename)){fn<-1}
savename<-nicDownload(filename[fn])	


#########################
# Load the grid point data 
#########################
# Read the feature class 5km grid sample points This will be the pre-attributed set of points with 227507
load(paste0(pathToGit,"data/studyarea_points_wNearLand.RData"))
#Get main projection that will be used 
primaryproj<-CRS(projection(studyarea_pointswLand))
#run function to get the fast ice data ready to be processed
myareas<-getFastIce(fileloc=savename,dataproj=primaryproj)


###############################
#Calculate fast ice width
##############################
FastIcePoints<-calcFasIceWidth(savename=savename,primaryproj=primaryproj,myareas=myareas,studyarea_pointswLand=studyarea_pointswLand,buffwidth=20000,plotit=TRUE)

sppdata<-getFastIceSPPdata(iceyear=iceyear,iceareas=myareas$fast,primaryproj=primaryproj,studyarea_pointswLand=studyarea_pointswLand)

FastIcePoints<-merge(FastIcePoints,sppdata,by.x="pointid",by.y="pointId",all.x=T)
fip<-merge(fip,sppdata,by.x="pointid",by.y="pointId",all.x=T)

save(FastIcePoints,file=paste0(resultsdir,"FastIcePoints_wFastIceCovars_",icemonth,iceyear,".RData"))
#save as shapefile too
writeOGR(obj=FastIcePoints, dsn=paste0(resultsdir,"FastIcePoints_wFastIceCovars_",icemonth,iceyear), layer=paste0("FastIcePoints",icemonth,iceyear), driver="ESRI Shapefile")


