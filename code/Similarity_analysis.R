## (Dis)Similarity analysis dataprep

#load libraries ---
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
                rename(spec = site) # better naming
    
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
    
    #create a logical index to remove rows that were merged twice due tot he joint count of 419
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
    west_seining <- west_seining_step1[logical_index,]
    
    #clean up the variables that can be used again later.   
    rm(west_seining_step1,temp,i,logical_index,repeat_sp)
    
    west_samples <- rbind(west_12S,west_16S,west_seining)%>%
      rename(species_code=spec)%>%
      mutate(coast="west")%>%
      data.frame()%>%
      select(marker,index,family,genus,scientific_name,common_name,species_code,
             BB,MB,SS4,SS3,WS,MuB,SB,PB,PL)
    
    