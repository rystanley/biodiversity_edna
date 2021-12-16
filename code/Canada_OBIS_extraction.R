## Extraction of Canadian species from OBIS

#load libraries
library(robis)
library(sf)
library(dplyr)
library(ggplot2)
library(taxize)
library(tidyr)

#function for extracting species taxonomic information
classify_row <- function(x,db="worms"){
  
  require(dplyr)
  require(tidyr)
  require(taxize)
  
  output <- taxize::classification(x,db=db,verbose=FALSE)%>%
            .[[1]]%>%
            as.data.frame(.,stringsAsFactors=FALSE)%>%
            dplyr::select(-id)%>%
            dplyr::filter(rank%in%c("Kingdom","Phylum","Class","Order",
                                     "Family","Genus","Species"))%>%
            tidyr::spread(rank,name)
    
  return(output)
  
}


#set a standard map projection
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

##load the regions - this requires a shape file that I can provide. 
bioregions <- st_read("data/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")%>%
  st_transform(latlong)%>%
  rowwise()%>%
  mutate(region=as.character(Legend),
         region=strsplit(region,"/",fixed=TRUE)[[1]][1],
         region=gsub("\\.","",region),
         region=gsub("[0-9]+","",region),
         region=trimws(region))%>%
  st_as_sf()%>%
  filter(REGION %in% c("Arctic","Pacific","Atlantic"))

#make bounding boxes
bb_arctic <- bioregions%>%
              filter(REGION == "Arctic")%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_as_sf()

bb_atlantic <- bioregions%>%
              filter(REGION == "Atlantic")%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_as_sf()

bb_pacific <- bioregions%>%
              filter(REGION == "Pacific")%>%
              st_bbox()%>%
              st_as_sfc()%>%
              st_as_sf()

#plot it to see what they look lik

bb_boxes <- rbind(bb_arctic%>%mutate(REGION="Arctic"),
                  bb_atlantic%>%mutate(REGION="Atlantic"),
                  bb_pacific%>%mutate(REGION="Pacific"))

ggplot()+geom_sf(data=bb_boxes,aes(fill=REGION))

##there is overlap to be removed
bb_arctic <- st_difference(bb_arctic,bb_pacific)
bb_arctic <- st_difference(bb_arctic,bb_atlantic)

ggplot()+geom_sf(data=bb_arctic%>%st_as_sf())


#Run OBIS extractions ------

fieldcodes <- c("year","date_year","country","decimalLongitude","decimalLatitude","aphiaID", 
                "kingdom","phylum","class","order","family","genus","species","depth","minimumDepthInMeters",
                "shoredistance")

if(!file.exists("data/OBIS/Arctic_obis.RData")){ #unless you want to re-run these extractions

  arctic_obis <- occurrence(geometry=st_as_text(bb_arctic$x),
                          fields = fieldcodes)
  
  if(!dir.exists("data/OBIS/")){dir.create("data/OBIS/")}
  
  save(arctic_obis,file="data/OBIS/Arctic_obis.RData") #note these are not tracked with github
  
}else(load("data/OBIS/Arctic_obis.RData"))

if(!file.exists("data/OBIS/Atlantic_obis.RData")){

  atlantic_obis <- occurrence(geometry=st_as_text(bb_atlantic$x),
                              fields = fieldcodes)
  
  if(!dir.exists("data/OBIS/")){dir.create("data/OBIS/")}
  
  save(atlantic_obis,file="data/OBIS/Atlantic_obis.RData")

}else(load("data/OBIS/Atlantic_obis.RData"))

if(!file.exists("data/OBIS/Pacific_obis.RData")){
  
  pacific_obis <- occurrence(geometry=st_as_text(bb_pacific$x),
                             fields = fieldcodes)
  
  if(!dir.exists("data/OBIS/")){dir.create("data/OBIS/")}
  
  save(pacific_obis,file="data/OBIS/Pacific_obis.RData")
  
}else(load("data/OBIS/Pacific_obis.RData"))


#Clean OBIS data -------

  #Arctic
      arctic_obis_sf <- arctic_obis%>%
                        st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
                        st_join(.,filter(bioregions,REGION == "Arctic")%>%
                                  dplyr::select(REGION,geometry)%>%
                                  st_buffer(0.5), #half degree buffer for the near shore critters
                                join = st_intersects)%>%
                        filter(!is.na(REGION))%>% #get rid of those outside the region
                        mutate(indspecies = paste(species,aphiaID,sep="_")) #for later matching
      
      arcticplot <- ggplot()+geom_sf(data=filter(bioregions,REGION == "Arctic"))+
        geom_sf(data=arctic_obis_sf%>%sample_frac(0.25),size=0.5)+ #plot a fraction of them just for visual inspection
        theme_bw()
      
      ggsave("data/OBIS/Arcticplot.png",arcticplot,height = 6,width=6,units="in",dpi=600)
      
      arctic_obis_taxonomy <- arctic_obis_sf%>%
        data.frame(.,stringsAsFactors = FALSE)%>%
        dplyr::select(-geometry)%>% #this gets rid of the spatial elements we don't need
        filter(!is.na(species) & !is.na(aphiaID))%>%
        distinct(aphiaID,.keep_all=TRUE)%>%
        arrange(species)%>%
        mutate(region="Arctic")%>%
        dplyr::select(region,kingdom,phylum,class,order,family,genus,species,aphiaID)
    
    #Atlantic
        atlantic_obis_sf <- atlantic_obis%>%
          st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
          st_join(.,filter(bioregions,REGION == "Atlantic")%>%
                    dplyr::select(REGION,geometry)%>%
                    st_buffer(0.5), #half degree buffer for the near shore critters
                  join = st_intersects)%>%
          filter(!is.na(REGION))%>% #get rid of those outside the region
          mutate(indspecies = paste(species,aphiaID,sep="_")) #for later matching
        
        atlanticplot <- ggplot()+geom_sf(data=filter(bioregions,REGION == "Atlantic"))+
          geom_sf(data=atlantic_obis_sf%>%sample_frac(0.25),size=0.5)+ #plot a fraction of them just for visual inspection
          theme_bw()
        
        ggsave("data/OBIS/Atlanticplot.png",atlanticplot,height = 6,width=6,units="in",dpi=600)
        
        atlantic_obis_taxonomy <- atlantic_obis_sf%>%
          data.frame(.,stringsAsFactors = FALSE)%>%
          dplyr::select(-geometry)%>% #this gets rid of the spatial elements we don't need
          filter(!is.na(species) & !is.na(aphiaID))%>%
          distinct(aphiaID,.keep_all=TRUE)%>%
          arrange(species)%>%
          mutate(region="Atlantic")%>%
          dplyr::select(region,kingdom,phylum,class,order,family,genus,species,aphiaID)
    
    #Pacific
        pacific_obis_sf <- pacific_obis%>%
          st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)%>%
          st_join(.,filter(bioregions,REGION == "Pacific")%>%
                    dplyr::select(REGION,geometry)%>%
                    st_buffer(0.5), #half degree buffer for the near shore critters
                  join = st_intersects)%>%
          filter(!is.na(REGION))%>% #get rid of those outside the region
          mutate(indspecies = paste(species,aphiaID,sep="_")) #for later matching
        
        pacificplot <- ggplot()+geom_sf(data=filter(bioregions,REGION == "Pacific"))+
          geom_sf(data=pacific_obis_sf%>%sample_frac(0.25),size=0.5)+ #plot a fraction of them just for visual inspection
          theme_bw()
        
        ggsave("data/OBIS/Pacificplot.png",pacificplot,height = 6,width=6,units="in",dpi=600)
        
        pacific_obis_taxonomy <- pacific_obis_sf%>%
          data.frame(.,stringsAsFactors = FALSE)%>%
          dplyr::select(-geometry)%>% #this gets rid of the spatial elements we don't need
          filter(!is.na(species) & !is.na(aphiaID))%>%
          distinct(aphiaID,.keep_all=TRUE)%>%
          arrange(species)%>%
          mutate(region="Pacific")%>%
          dplyr::select(region,kingdom,phylum,class,order,family,genus,species,aphiaID)

 ## assemble it at the end
        Canada_obis <- rbind(arctic_obis_taxonomy,atlantic_obis_taxonomy,pacific_obis_taxonomy)
        write.csv(Canada_obis,"data/OBIS/Canada_obis_taxa.csv",row.names = FALSE)
        