# FastIceCovars
Scripts to generate the fast ice covariate data for the Seals from Space project.  
For more information, contact Dennis Jongsomjit [djongsomjit@pointblue.org] or  
Leo Salas [lsalas@p0intblue.org]

## Overall approach  
We start with a raster dataset of bathymetric data. From this raster we construct a grid of the desired resolution  
(in our case 5 km cells). This grid is then used to sample NIC data. All but the ice data are treated here as immutable:  
mean depth, distance to shore, distance to the 300 m depth isobath, mean slope, distance to ADPE and EMPE colonies,  
distance to glacier.

We provide a function to query the NIC database and then choose a specific date. Once chosen, the function filters  
the grid for only points on fast ice, and retrieves distance to the edge, fast ice width and a simple metric of  
how good the fast ice is around the point.

### Preliminary work in ArcGIS
A bathymetric grid of the southern ocean was obtained at a 500m resolution (IBCSO v1.0; Arndt et al. 2013).  
Within ArcMap a land and ice shelf layer (ADD v2.0; British Antarctic Survey 1998) was used to mask out these  
areas from the bathymetric grid.

From this layer then:
* A bathymetry line shapefile at -300 m was created using the Contour tool of ArcGIS.  
* Slope was calculated in degrees using the Slope tool of ArcGIS. 
* To derive the mean depth, and mean slope, the bathymetric grid and slope grid were averaged across 
   10 x 10 500-m cells using the Aggregate tool of ArcGIS. The aggregate tool allows for the creation of grids at  
   different resolutions. Because we aggregated by 10 x 10 cells, the resulting grid is of 5-km resolution. 
* The 5-km grids was then used to create the 5-km sampling location shapefile (i.e. the cell centroids from the raster)  
   using the Raster to Point tool of ArcGIS.  These sampling locations can be filtered to locations of interest.  
   For the purpose of this work, we used areas where fast ice is potentially present. 
* We limited the points to a minimum latitude of 60 degees South, and a maximum determined by the land-plus-ice shelf mask.
* The points were thus attributed with mean slope and bathymetry. We then calculated distance to the 300-m isobath  
   using the XX tool of ArcGIS. (Dennis to complete)
     
The resulting shapefile is ready for use in the R environment.

### Filtering and attribution of points using R 
The attributed shapefile is read into the R environment. Within R, each sampling location was further attributed with  
spatially overlapping grid values or distances to shapefile features.

NEED to explain what each file does... (Leo to do)

## Projection
For processing, all objects were in or re-projected to the Antarctic Polar Stereographic WGS84 projection  
  
Proj4 = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m  
            +no_defs +ellps=WGS84 +towgs84=0,0,0"
 
## Citations
- Arndt, J.E., H. W. Schenke, M. Jakobsson, F. Nitsche, G. Buys, B. Goleby, M. Rebesco, F. Bohoyo, J.K. Hong, J. Black,  
   R. Greku, G. Udintsev, F. Barrios, W. Reynoso-Peralta, T. Morishita, R. Wigley, "The International Bathymetric Chart  
   of the Southern Ocean (IBCSO) Version 1.0 - A new bathymetric compilation covering circum-Antarctic waters",  
   Geophysical Research Letters, doi: 10.1002/grl.50413

- British Antarctic Survey. 1998. Antarctic Digital Database, Version 2.0. Manual and bibliography. Scientific Committee  
   on Antarctic Research, Cambridge. 74 pp. https://www.add.scar.org/
