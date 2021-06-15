##First order adjustment function function 'coord_bump' ---------------------
coord_bump <- function(x,ras,radius = 10, printmessage=F){
  
  require(dplyr)
  require(sf)
  require(raster)
  require(nngeo)
  
  # x = dataframe with a column for lat and long that has the same projection as ras
  # ras = bathymetry raster
  # radius = search radius in km (default is 10 or 10000 m)
  # printmessage = logical whether you want the current port to be printed. 
    #this comes before the calculations so is good for tracking down bugs
  
  if(printmessage){print(paste0("Working on ",as.character(x$portname)))}
  
  rasproj <- proj4string(ras)
  
  x <- st_as_sf(x,coords=c("long","lat"),crs=rasproj) #convert to sf
  x2 <- x%>%st_transform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")
  
  #set up a dynamic azimutal equidistance projection for each point 
  #https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world/121539#121539
  
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  st_coordinates(x2)[2], st_coordinates(x2)[1])
  
  cr <- x%>%
        st_transform(aeqd)%>%
        st_buffer(radius*1000)%>% # km search radius for water
        st_as_sf()%>%
        st_transform(proj4string(ras))%>% #Back to the projection of the raster
        extent()

  bathy <- crop(ras,cr)
  
  xy <- st_coordinates(x)%>%SpatialPoints()
  
  depth <- raster::extract(bathy,xy)

  if(depth>0 | is.na(depth)){
    
   bathypoints <- coordinates(bathy)[values(bathy)<0,]%>%
                  data.frame()%>%
                  st_as_sf(coords=c("x","y"),crs=rasproj)
                  
  nearest_point <- bathypoints[suppressMessages(st_nn(x,bathypoints,progress = FALSE))[[1]],]
  
  return(st_coordinates(nearest_point)%>%data.frame()%>%rename(lon_a=X,lat_a=Y))
  
  } else {return(st_coordinates(x)%>%data.frame()%>%rename(lon_a=X,lat_a=Y))}
      
}

##Secondary adjustment function function 'coord_bump2' ---------------------
coord_bump2 <- function(x,ras,radius = 10, printmessage=F){
  
  #This will work if the initial cost-distance returns an INF meaning the coarse resoluiton
  #transition layer landlocked the destination coordinate which itself was assigned a coordinate using
  #coord_bump but was not in the 'open' water. This process is unsupervised so could cause issues depending on the 
  #topography of the seascape tested and, in particular, for those locations near an narrow peice of land
  #whereby the selection could go to the opposite side. These may need to be manually checked. 
  
  ##Variables --- 
  # x = dataframe with a column for adjusted 'lon_a' and 'lat_a" after coord_bump -- projection latlong
  # ras = bathymetry raster if null one will be created using the r natural earth shapefile
  # radius = search radius in km (default is 10 or 10000 m)
  # printmessage = logical whether you want the current port to be printed. 
  #this comes before the calculations so is good for tracking down bugs
  
  #returned nearest variable in the projection of the raster. 
  
  #load libraries -----
  require(dplyr)
  require(sf)
  require(raster)
  require(nngeo)
  require(units)
  
  #progress message ----
  if(printmessage){print(paste0("Working on ",as.character(x$portname)))}
  
  #set up a dynamic projection ----
  
  #latitude and longitude projection 
  latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  
  x = st_as_sf(x,coords=c("lon_a","lat_a"),crs=latlong) #convert to sf
  
  #set up a dynamic azimutal equidistance projection for each point 
  #https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world/121539#121539
  
  aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                  st_coordinates(x)[2], st_coordinates(x)[1])
  
  
  #find the geometric centre of water in a radius from the point
  rasproj <- proj4string(ras)
  
  cr <- x%>%
    st_transform(aeqd)%>%
    st_buffer(.,dist=units::set_units(radius,"km"))%>% # km search radius for water
    st_as_sf()%>%
    st_transform(proj4string(ras))%>%
    extent()
  
  #create a small bathy
  bathy <- crop(ras,cr)
  
  #use clumping to identify the cluster of contiguous points using raster::climp
  #Here I take advantage of the cluster seperator. I assign all non water values as NA
  #The cluster function will identify contiguous clusters. Because of the small spatial scale
  #I can pick the cluster with the most observations (using table) as the largest and then find the 
  #nearest coordinate assocaited with that cluster
  
  bathy2 <- bathy
  bathy2[bathy2==10] <- NA
  
  clumpras <- raster::clump(bathy2)
  
  #identify the biggest cluster
  bc <- clumpras%>%
    values(.)%>%
    table()%>%
    data.frame()%>%
    rename(cluster=1)%>%
    arrange(-Freq)%>%
    .[1,"cluster"]%>%
    as.numeric()
  
  #identify the points of interest that correspond to that cluster
  bathypoints <- coordinates(bathy)[which(values(clumpras)==bc),]%>%
    data.frame()%>%
    st_as_sf(coords=c("x","y"),crs=rasproj)
  
  #identify the nearest point within the large water cluster
  nearest_point <- bathypoints[suppressMessages(st_nn(x,bathypoints,progress = FALSE))[[1]],]
  
  #Note that this is unsupervised so there could be situations where narrow land points may push the 
  #coordinate to the other side of a land barrier. This is why you should start with small search radii
  
  return(st_coordinates(nearest_point)%>%data.frame()%>%rename(lon_a=X,lat_a=Y))
  
}   

coord_bump_trans <- function(x,xproj,ras,radius = 10, printmessage=F){
  
  #This will work if the initial cost-distance returns an INF meaning the coarse resoluiton
  #transition layer (converted to a raster) landlocked the destination coordinate which itself was assigned a coordinate using
  #coord_bump but was not in the 'open' water. This process is unsupervised so could cause issues depending on the 
  #topography of the seascape tested and, in particular, for those locations near an narrow peice of land
  #whereby the selection could go to the opposite side. These may need to be manually checked. 
  
  ##Variables --- 
  # x = dataframe with lat and long variables. This can come after the initial coord_bump
  # xproj - is the projection of the raster that was used to create the trans object and of the x coordinats (these must be the same. )
  # trans = is the transition object
  # radius = search radius in km (default is 10 or 10000 m)
  # printmessage = logical whether you want the current port to be printed. 
  # this comes before the calculations so is good for tracking down bugs
  
  #returned nearest variable in the projection of the raster. 
  
  #load libraries -----
  require(dplyr)
  require(sf)
  require(raster)
  require(nngeo)
  require(units)
  
  #progress message ----
  if(printmessage){print(paste0("Working on ",as.character(x$name)))}
  
  x <- st_as_sf(x,coords=c("long","lat"),crs=xproj)
  
  #check if a adjustement is needed
  xval <- raster::extract(ras,as_Spatial(x))
  
  if(is.nan(xval)){
    
    #set up a dynamic projection ----
    x1 = x%>%st_transform("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0") #convert to sf
    
    #set up a dynamic azimutal equidistance projection for each point 
    #https://gis.stackexchange.com/questions/121489/1km-circles-around-lat-long-points-in-many-places-in-world/121539#121539
    
    aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    st_coordinates(x1)[2], st_coordinates(x1)[1])
    
    #crop the raster to speed up the search
    cr <- x%>%
      st_transform(aeqd)%>%
      st_buffer(.,dist=units::set_units(radius,"km"))%>% # km search radius for water
      st_as_sf()%>%
      st_transform(xproj)%>%
      extent()
    
    #create a small bathy
    bathy <- crop(ras,cr)
    
    #use clumping to identify the cluster of contiguous points using raster::climp
    #Here I take advantage of the cluster seperator. I assign all non water values as NA
    #The cluster function will identify contiguous clusters. Because of the small spatial scale
    #I can pick the cluster with the most observations (using table) as the largest and then find the 
    #nearest coordinate assocaited with that cluster
    
    bathy2 <- bathy
    bathy2[is.nan(bathy[])] <- NA
    
    clumpras <- raster::clump(bathy2)
    
    #identify the biggest cluster
    bc <- clumpras%>%
      values(.)%>%
      table()%>%
      data.frame()%>%
      rename(cluster=1)%>%
      arrange(-Freq)%>%
      .[1,"cluster"]%>%
      as.numeric()
    
    #identify the points of interest that correspond to that cluster
    bathypoints <- coordinates(bathy)[which(values(clumpras)==bc),]%>%
      data.frame()%>%
      st_as_sf(coords=c("x","y"),crs=xproj)
    
    #identify the nearest point within the large water cluster
    nearest_point <- bathypoints[suppressMessages(st_nn(x,bathypoints,progress = FALSE))[[1]],]
    
    #Note that this is unsupervised so there could be situations where narrow land points may push the 
    #coordinate to the other side of a land barrier. This is why you should start with small search radii
    
    #return output with the euclidean distance between the start point and the offset noted. 
    
    return(st_coordinates(nearest_point)%>%
             data.frame()%>%
             rename(lon_a=X,lat_a=Y)%>%
             mutate(dist_offset=as.numeric(st_distance(x,nearest_point))/1000))
    
  } #end first if   
  
  if(!is.nan(xval)){
    
    return(st_coordinates(x)%>%
             data.frame()%>%
             rename(lon_a=X,lat_a=Y)%>%
             mutate(dist_offset=0))
    
  }# end no adjustment needed if
  
}#end funciton
