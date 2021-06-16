#load libraries --------
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(tidyr)
library(dplyr)
library(raster)
library(ggplot2)
library(gdistance)
library(nngeo)
library(fasterize)
library(reshape2)

#Projections -----
  latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Source functions for analysis ------------
  source("code/coord_bump.R")
  source("code/lcp_function.R")

#Source basemap -------------- 
  basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                     dplyr::select(name_en,geometry)%>%
                     st_as_sf()%>%
                     st_union()%>%
                     st_transform(latlong)%>%
                     st_as_sf()%>%
                     mutate(country="Canada"),
                   ne_states(country = "United States of America",returnclass = "sf")%>%
                     dplyr::select(name_en,geometry)%>%
                     st_as_sf()%>%
                     st_union()%>%
                     st_transform(latlong)%>%
                     st_as_sf()%>%
                     mutate(country="USA"))%>%st_union()%>%st_as_sf()%>%suppressMessages()

#load geographic coordinates for sites and match up format between east and west coast ---------
  
  east <- read.csv("data/Site coordinates/EastCoast.csv")%>%
          mutate(Longitude=Longitude*-1)%>%
          st_as_sf(coords=c("Longitude","Latitude"),crs=latlong,remove=FALSE)%>%
          rename(lat=Latitude,
                 long=Longitude,
                 site_id=Site_Code)%>%
          mutate(coast="east")
    
  west <- read.csv("data/Site coordinates/WestCoast.csv")%>%
          st_as_sf(coords=c("CTD_Long","CTD_Lat"),crs=latlong,remove=FALSE)%>%
          rename(lat=CTD_Lat,
                 long=CTD_Long,
                 site_id=Site.ID)%>%
        mutate(coast="west")

##Group coordinates into a single 'sf' data.frame
  sample_sites <- rbind(east%>%dplyr::select(coast,site_id,geometry),
                      west%>%dplyr::select(coast,site_id,geometry))%>%
                mutate(long=st_coordinates(.)[,1],
                       lat=st_coordinates(.)[,2])%>%
                st_transform(latlong)

##least cost path analysis using custom function--------------
  east_lcp <- lcp_function(sample_sites%>%filter(coast=="east"),basemap=basemap)
  west_lcp <- lcp_function(sample_sites%>%filter(coast=="west"),basemap=basemap)
    
#Great Circle Distance (as the bird flies) -------------------
  east_gcd <- as.matrix(st_distance(east)/1000)%>%
              matrix(., ncol=nrow(east), nrow=nrow(east))

      rownames(east_gcd) <- east%>%data.frame()%>%pull(site_id) #adjust column/row names
      colnames(east_gcd) <- east%>%data.frame()%>%pull(site_id) 

    west_gcd <- as.matrix(st_distance(west)/1000)%>%
        matrix(., ncol=nrow(west), nrow=nrow(west))
      
    rownames(west_gcd) <- west%>%data.frame()%>%pull(site_id) #adjust column names
    colnames(west_gcd) <- west%>%data.frame()%>%pull(site_id) 
    
#model comparisons
    
    plot_comparison_df <- rbind(data.frame(gcd=east_gcd[lower.tri(east_gcd)],lcd=east_lcp[lower.tri(east_lcp)],coast="East coast"),
                                data.frame(gcd=west_gcd[lower.tri(west_gcd)],lcd=west_lcp[lower.tri(west_lcp)],coast="West coast"))
    
    ggplot(data=plot_comparison_df,aes(x=gcd,y=lcd))+
      geom_point()+
      geom_abline(slope=1,intercept=0,lty=2)+
      geom_smooth(method="lm")+
      facet_grid(~coast)+
      labs(x="Great cirle Distance",y="Least-cost path")+
      theme_bw()
    
    #East coast linear model
    lm(east_gcd[lower.tri(east_gcd)]~east_lcp[lower.tri(east_lcp)])%>%summary()
    
    #west coast linear model
    lm(west_gcd[lower.tri(west_gcd)]~west_lcp[lower.tri(west_lcp)])%>%summary()

#Write the outputs ------------
    write.csv(east_lcp,"output/east_lcp_distances.csv",row.names = TRUE)
    write.csv(east_gcd,"output/east_gcd_distances.csv",row.names = TRUE)
    
    write.csv(west_lcp,"output/west_lcp_distances.csv",row.names = TRUE)
    write.csv(west_gcd,"output/west_gcd_distances.csv",row.names = TRUE)
