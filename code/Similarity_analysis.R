## (Dis)Similarity analysis dataprep

#load libraries ---
library(tidyverse)
library(dplyr)
library(sf)
library(vegan)
library(tidyr)
library(readxl)



### Data Prep ----------------

#East coast data prep -----

  #read in the data with the official nomencalutre to be matched with the site level data
  
  east_siening_nomenclature <- read_excel("data/Biodiveristy_Coastal.xlsx",sheet="east_seining")%>% #from seines
    data.frame()%>%
    rename(count = Number.of.fish)
  
  colnames(east_siening_nomenclature) <- gsub('\\.', '_',colnames(east_siening_nomenclature))%>%tolower()
  
  east_edna_nomenclature <- read_excel("data/Biodiveristy_Coastal.xlsx",sheet="east_edna")%>% #from eDNA
    data.frame()
  
  colnames(east_edna_nomenclature) <- gsub('\\.', '_',colnames(east_edna_nomenclature))%>%tolower()


  #East coast catch by site

    east_12S_step1 <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="12S_east")%>% #10 sites in the east
                dplyr::select(1:11)%>%
                mutate(count = rowSums(.[,2:11]),
                       marker="12S")%>%
                suppressMessages()%>%
                left_join(.,east_edna_nomenclature%>%select(-e_16s_reads)%>%rename(count=e_12s_reads)%>%filter(!is.na(count)))%>%
                data.frame()%>%
                rename(spec = site) # better naming
                
      #summed counts for  "Chaetodontidae_Chaetodon_ocellatus_12S" "Clupeidae_Alosa_aestivalis_12S" are both 419 so we need a filter
      repeat_sp <- east_12S_step1%>%group_by(spec)%>%summarise(count=n())%>%filter(count>1)%>%pull(spec)
      
      #create a logical index to remove rows that were merged twice due tot he joint count of 419
        logical_index=rep(TRUE,nrow(east_12S_step1))  
        for(i in 1:nrow(east_12S_step1)){
          temp=east_12S_step1[i,]
              if(temp$spec %in% repeat_sp & grepl(temp$genus,temp$spec)){logical_index[i] <- FALSE} # mismatches can be found by id'ing those whose genus does not match the full name
        }
        
      #apply index to get rid of issue      
        east_12S <- east_12S_step1[logical_index,]
        
      #clean up the variables that can be used again later.   
        rm(east_12S_step1,temp,i,logical_index,repeat_sp)
    
    #16S has no duplicate counts so it merges without issue
   east_16S <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="16S_east")%>% 
                dplyr::select(1:11)%>%
                mutate(count = rowSums(.[,2:11]),
                       marker="16S")%>%
                suppressMessages()%>%
                left_join(.,east_edna_nomenclature%>%select(-e_12s_reads)%>%rename(count=e_16s_reads)%>%filter(!is.na(count)))%>%
                data.frame()%>%
                rename(spec = site)%>% # better naming
                mutate(scientific_name = ifelse(scientific_name == "Ammodytes hexapterus | Ammodytes personatus",genus,scientific_name)) #two groups have the same grouping but a different genus name (AM2 and AM3)
    
   #Seining has 4 species that have the count of 1
    east_seining_step1 <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="seining_east")%>% 
                    dplyr::select(1:11)%>%
                    mutate(count = rowSums(.[,2:11]),
                           marker="Seining")%>%
                    suppressMessages()%>%
                    rename(WRK02 = Wreck,
                           CON = Conrod,
                           TAE = TaylorEast,
                           TAW = TaylorWest,
                           MOS = Moosehead,
                           FAI = FiftyAcre, 
                           RSE = RoseBay,
                           CAB = Cable,
                           FRK = FranksGeorge,
                           PLS = PointPleasant)%>%
                    left_join(.,east_siening_nomenclature)%>%
                    rename(spec = site)%>% # better naming
                    dplyr::select(names(east_16S))%>%
                    data.frame() 
    
    #summed counts for "AmericanEel"  "Lumpfish" "MackerelScad" "RoughScad" are all 1 so we need a filter
    repeat_sp <- east_seining_step1%>%group_by(spec)%>%summarise(count=n())%>%filter(count>1)%>%pull(spec)
    
    #create a logical index to remove rows that were merged twice due tot he joint count of 419
    logical_index=rep(TRUE,nrow(east_seining_step1))  
    
    for(i in 1:nrow(east_seining_step1)){
      
      temp=east_seining_step1[i,]%>%
           mutate(common_name = tolower(gsub(" ","",common_name)),
                  spec2=tolower(spec))
      
      if(temp$spec %in% repeat_sp & temp$common_name != temp$spec2){logical_index[i] <- FALSE} # mismatches can be found by id'ing those whose genus does not match the full name
    }
    
    #apply index to get rid of issue      
    east_seining <- east_seining_step1[logical_index,]
    
    #clean up the variables that can be used again later.   
    rm(east_seining_step1,temp,i,logical_index,repeat_sp)
      
    
  east_samples <- rbind(east_12S,east_16S,east_seining)%>%
                    rename(species_code=spec)%>%
                    mutate( coast="east")%>%
                    data.frame()%>%
                    select(marker,index,family,genus,scientific_name,common_name,species_code,
                           TAW,PLS,TAE,FRK,RSE,MOS,FAI,WRK02,CAB,CON)
                    

#West coast data prep ---------
    west_siening_nomenclature <- read_excel("data/Biodiveristy_Coastal.xlsx",sheet="west_seining")%>%
      data.frame()%>%
      rename(count = Number.of.fish)
    
    colnames(west_siening_nomenclature) <- gsub('\\.', '_',colnames(west_siening_nomenclature))%>%tolower()
    
    west_edna_nomenclature <- read_excel("data/Biodiveristy_Coastal.xlsx",sheet="west_edna")%>%
      data.frame()
    
    colnames(west_edna_nomenclature) <- gsub('\\.', '_',colnames(west_edna_nomenclature))%>%tolower()
    
    
    #west coast catch by site
    west_12S <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="12S_west")%>% #10 sites in the west
      dplyr::select(1:10)%>%
      mutate(count = rowSums(.[,2:10]),
             marker="12S")%>%
      suppressMessages()%>%
      left_join(.,west_edna_nomenclature%>%select(-w_16s_reads)%>%rename(count=w_12s_reads)%>%filter(!is.na(count)))%>%
      data.frame()%>%
      rename(spec = site, # better naming
             PB = PB_e) # not sure why this is named PB_e
    
    west_16S <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="16S_west")%>% 
      dplyr::select(1:10)%>%
      mutate(count = rowSums(.[,2:10]),
             marker="16S")%>%
      suppressMessages()%>%
      left_join(.,west_edna_nomenclature%>%select(-w_12s_reads)%>%rename(count=w_16s_reads)%>%filter(!is.na(count)))%>%
      data.frame()%>%
      rename(spec = site, # better naming
             PB = PB_e) #not sure why this is named PB_e
    
    west_seining_step1 <- read_excel("data/iNEXT/data.iNEXT.xlsx",sheet="seining_west")%>% 
      dplyr::select(1:10)%>%
      mutate(count = rowSums(.[,2:10]),
             marker="Seining")%>%
      suppressMessages()%>%
      left_join(.,west_siening_nomenclature)%>%
      mutate(common_name = gsub("/"," ",common_name))%>%
      rename(spec = site)%>% # better naming
      dplyr::select(names(west_16S))%>%
      data.frame() 
    
    #summed counts for various species cause some issues. 
    repeat_sp <- west_seining_step1%>%group_by(spec)%>%summarise(count=n())%>%filter(count>1)%>%pull(spec)
    
    #create a logical index to remove rows that were merged twice due to the joint count of 419
    logical_index=rep(TRUE,nrow(west_seining_step1))  
    
    for(i in 1:nrow(west_seining_step1)){
      #step that will make the common_name match up to the common_name of the seine. 
      temp=west_seining_step1[i,]%>%
        mutate(common_name = tolower(gsub(" ","",common_name)), #fixes to some naming inconsistencies
               common_name = gsub("\\.","",common_name),
               common_name = gsub("-","",common_name),
               spec2=gsub("_"," ",spec)%>%tolower(),
               spec2=gsub(" ","",spec2),
               spec2=gsub("kelpsurfperch","kelpperch",spec2))#this was a change in the common name I noticed. 
      
      if(temp$spec %in% repeat_sp & temp$common_name != temp$spec2){logical_index[i] <- FALSE} # mismatches can be found by id'ing those whose genus does not match the full name
    }
    
    #apply index to get rid of issue      
    west_seining <- west_seining_step1[logical_index,]%>%
                    mutate(scientific_name=ifelse(spec =="Copper_QuillbackRockfish","Sebastes caurinus",scientific_name))#for some reason the full name of copper quillback rockfish was not spelled out. 
    
    #clean up the variables that can be used again later.   
    rm(west_seining_step1,temp,i,logical_index,repeat_sp)
    
    west_samples <- rbind(west_12S,west_16S,west_seining)%>%
      rename(species_code=spec)%>%
      mutate(coast="west")%>%
      data.frame()%>%
      select(marker,index,family,genus,scientific_name,common_name,species_code,
             PB,SS3,SS4,MB,SB,BB,PL,WS,MuB) # order of stations taken from the distance calculations. 
    
### Similarity analysis -----------------
    
   #read in the distances calculations
    west_gc <- read.csv("output/west_gcd_distances.csv")
    west_siteorder <- west_gc$X## this will be used to better line up the comparisons 
    rownames(west_gc) <- west_gc$X
    west_gc <- west_gc%>%dplyr::select(-X)
    
    west_lcd <- read.csv("output/west_lcp_distances.csv")
    rownames(west_lcd) <- west_lcd$X
    west_lcd <- west_lcd%>%dplyr::select(-X)
    
    east_gc <- read.csv("output/east_gcd_distances.csv")
    east_siteorder <- east_gc$X
    rownames(east_gc) <- east_gc$X
    east_gc <- east_gc%>%dplyr::select(-X)
    
    east_lcd <- read.csv("output/east_lcp_distances.csv")
    rownames(east_lcd) <- east_lcd $X
    east_lcd  <- east_lcd %>%dplyr::select(-X)
    
    
    #function to create the sample by species dataframes required to do the analysis
    simmilarity_format <- function(x,siteorder,PA=FALSE){
      
      x <- x%>%
        dplyr::select(c("scientific_name",siteorder))%>%
        mutate(scientific_name = gsub(" ","_",scientific_name))%>%
        gather("site","val",-scientific_name)%>%
        pivot_wider(names_from = scientific_name,values_from = val)%>%
        data.frame()
      
      rownames(x) <- x$site
      
      x <- x%>%dplyr::select(-site)
      
      #convert to binary (presence/absence)
      if(PA){x[x>0] <- 1} 
      
      return(x)
    }
  
  ## West Coast analysis
    
    #Jaccard 
    west_12S_sim_jaccard <- filter(west_samples,marker=="12S")%>%
                    simmilarity_format(.,siteorder = west_siteorder,PA=TRUE)%>%
                    vegdist(.,method="jaccard")
    
    west_16S_sim_jaccard <- filter(west_samples,marker=="16S")%>%
                    simmilarity_format(.,siteorder = west_siteorder,PA=TRUE)%>%
                    vegdist(.,method="jaccard")
    
    west_seining_sim_jaccard <- filter(west_samples,marker=="Seining")%>%
                        simmilarity_format(.,siteorder = west_siteorder,PA=TRUE)%>%
                        vegdist(.,method="jaccard")
    
    #Bray-Curtis
    west_12S_sim_bray <- filter(west_samples,marker=="12S")%>%
      simmilarity_format(.,siteorder = west_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    west_16S_sim_bray <- filter(west_samples,marker=="16S")%>%
      simmilarity_format(.,siteorder = west_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    west_seining_sim_bray <- filter(west_samples,marker=="Seining")%>%
      simmilarity_format(.,siteorder = west_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    
 ## East Coast analysis
    
    #Jaccard
    east_12S_sim_jaccard <- filter(east_samples,marker=="12S")%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=TRUE)%>%
      vegdist(.,method="jaccard")
    
    east_16S_sim_jaccard <- filter(east_samples,marker=="16S")%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=TRUE)%>%
      vegdist(.,method="jaccard")
    
    east_seining_sim_jaccard <- filter(east_samples,marker=="Seining")%>%
      mutate(scientific_name = ifelse(is.na(scientific_name),paste(family,genus,sep=""),scientific_name), #some aren't to the species level but they are the only ones wihtin that family/genus so we can count as an ~OTU
             scientific_name = gsub("daeNA","dae",scientific_name))%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=TRUE)%>%
      vegdist(.,method="jaccard")
    
    #Bray-Curtis
    east_12S_sim_bray <- filter(east_samples,marker=="12S")%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    east_16S_sim_bray <- filter(east_samples,marker=="16S")%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    east_seining_sim_bray <- filter(east_samples,marker=="Seining")%>%
      mutate(scientific_name = ifelse(is.na(scientific_name),paste(family,genus,sep=""),scientific_name), #some aren't to the species level but they are the only ones wihtin that family/genus so we can count as an ~OTU
             scientific_name = gsub("daeNA","dae",scientific_name))%>%
      simmilarity_format(.,siteorder = east_siteorder,PA=FALSE)%>%
      decostand(.,method="hellinger")%>%
      vegdist(.,method="bray")
    
    
    #Stitch it together
    west_site_combo<- t(combn(colnames(as.matrix(west_seining_sim_jaccard)), 2))
    east_site_combo<- t(combn(colnames(as.matrix(east_seining_sim_jaccard)), 2))
    
    west_sim <- rbind(west_site_combo%>%
                   data.frame()%>%
                   rename(start=1,end=2)%>%
                   mutate(primer_12S = as.matrix(west_12S_sim_jaccard)[west_site_combo],
                          primer_16S = as.matrix(west_16S_sim_jaccard)[west_site_combo],
                          Seining = as.matrix(west_seining_sim_jaccard)[west_site_combo],
                          great_circle = as.matrix(west_gc)[west_site_combo],
                          least_cost_distance = as.matrix(west_lcd)[west_site_combo],
                          coast="Pacific",
                          dissim = "Jaccard"),
                   west_site_combo%>%
                     data.frame()%>%
                     rename(start=1,end=2)%>%
                     mutate(primer_12S = as.matrix(west_12S_sim_bray)[west_site_combo],
                            primer_16S = as.matrix(west_16S_sim_bray)[west_site_combo],
                            Seining = as.matrix(west_seining_sim_bray)[west_site_combo],
                            great_circle = as.matrix(west_gc)[west_site_combo],
                            least_cost_distance = as.matrix(west_lcd)[west_site_combo],
                            coast="Pacific",
                            dissim = "Bray-Curtis"))
    
    east_sim <- rbind(east_site_combo%>%
                        data.frame()%>%
                        rename(start=1,end=2)%>%
                        mutate(primer_12S = as.matrix(east_12S_sim_jaccard)[east_site_combo],
                               primer_16S = as.matrix(east_16S_sim_jaccard)[east_site_combo],
                               Seining = as.matrix(east_seining_sim_jaccard)[east_site_combo],
                               great_circle = as.matrix(east_gc)[east_site_combo],
                               least_cost_distance = as.matrix(east_lcd)[east_site_combo],
                               coast="Atlantic",
                               dissim = "Jaccard"),
                      east_site_combo%>%
                        data.frame()%>%
                        rename(start=1,end=2)%>%
                        mutate(primer_12S = as.matrix(east_12S_sim_bray)[east_site_combo],
                               primer_16S = as.matrix(east_16S_sim_bray)[east_site_combo],
                               Seining = as.matrix(east_seining_sim_bray)[east_site_combo],
                               great_circle = as.matrix(east_gc)[east_site_combo],
                               least_cost_distance = as.matrix(east_lcd)[east_site_combo],
                               coast="Atlantic",
                               dissim = "Bray-Curtis"))
    
    #one unified dataframe
    simmilarity_df <- rbind(west_sim,east_sim)
    
    
### Plot the data -----------
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    
    ggplot(data=simmilarity_df%>%
             gather("marker","val",primer_12S,primer_16S,Seining),
           aes(x=great_circle,y=val,col=marker,group=marker))+
      geom_point()+
      facet_grid(coast~dissim)+
      theme_bw()+
      stat_smooth(method="lm")
    
    ggplot(data=simmilarity_df%>%gather("marker","val",primer_12S,primer_16S,Seining),
           aes(x=least_cost_distance,y=val,col=marker,group=marker))+
      geom_point()+
      facet_grid(coast~dissim)+
      theme_bw()+
      stat_smooth(method="lm")
    
    mod_lcp <- glm(great_circle~marker*coast,data=simmilarity_df%>%gather("marker","val",primer_12S,primer_16S,Seining)%>%filter(dissim=="Jaccard"))
    
    summary(mod_lcp)
    summary(mod_lcp)$coefficient
    