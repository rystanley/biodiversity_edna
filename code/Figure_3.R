## This code generates Figure 3

#load libraries -------------
    #install.packages('devtools')
    #library(devtools)
    #install_github('JohnsonHsieh/iNEXT')
    
    #devtools::install_github("vmikk/metagMisc")
    #library(metagMisc)
    
    library(iNEXT)
    library(ggplot2)
    library(gridExtra)
    library(patchwork)
    library(dplyr)

## Richness analysis ----------

  #read in data
    east1 <-read.table("data/iNEXT/data.iNEXT.east.txt",header=T)
    west1 <-read.table("data/iNEXT/data.iNEXT.west.txt",header=T)
  
  #set number of samples 
    east.num <- 20
    west.num <- 18
  
  #Richness analysis
    east.inc <- iNEXT(east1, q=0, datatype="incidence_freq", endpoint=east.num, knots = 20, se = TRUE, conf = 0.95,
                     nboot = 1000)
    
    west.inc <- iNEXT(west1, q=0, datatype="incidence_freq", endpoint=west.num, knots = 18, se = TRUE, conf = 0.95,
                      nboot = 1000)

## Plotting----------

          #assemble data from the iNEXT function using the handy formatting function 'fortify' modified by the iNEXT package
            plotdata <- rbind(fortify.iNEXT(east.inc,type=1),fortify.iNEXT(east.inc,type=3),
                              fortify.iNEXT(west.inc,type=1),fortify.iNEXT(west.inc,type=3))%>%
              mutate(region=ifelse(grepl("east",site),"Atlantic","Pacific"),
                     type=gsub("east_","",site),
                     type=gsub("west_","",type),
                     type=gsub("seining","Seining",type),
                     type=factor(type,levels=c("Seining","12S","16S","eDNA")),
                     denominator = ifelse(plottype==1,"Number of sampling sites","Sample coverage"),
                     region=factor(region,levels=c("Atlantic","Pacific")))%>%
              arrange(region,type,x)%>%
              dplyr::select(region,type,method,x,y,y.lwr,y.upr,denominator)
            
            #assemble data for the asymptotic diversity estimates
            assym_div <- rbind(east.inc$AsyEst,west.inc$AsyEst)%>%
              mutate(region=ifelse(grepl("east",Site),"Atlantic","Pacific"),
                     type=gsub("east_","",Site),
                     type=gsub("west_","",type),
                     type=gsub("seining","Seining",type),
                     type=factor(type,levels=c("Seining","12S","16S","eDNA")))%>%
              dplyr::select(-Site)
      
            #Rarefaction plots
            p1 <- ggplot()+
              geom_ribbon(data=plotdata,aes(x=x,ymin=y.lwr,ymax=y.upr,fill=type),alpha=0.5)+
              geom_line(data=plotdata,aes(x=x,y=y,col=type),lty=2,lwd=1.25)+
              geom_line(data=filter(plotdata,method!="extrapolated"),aes(x=x,y=y,col=type),lwd=1.25)+
              geom_point(data=filter(plotdata,method=="observed"),aes(x=x,y=y,col=type),size=4)+
              theme_bw()+
              facet_grid(region~denominator,scales="free_x")+
              theme(strip.background = element_rect(, colour = "black", fill = "white"),
                    strip.text.x = element_text(colour = "black",size=14), 
                    strip.text.y = element_text(colour = "black",size=14),
                    axis.text = element_text(colour = "black",size=12),
                    axis.title = element_text(colour = "black",size=12),
                    legend.position="bottom")+
              scale_y_continuous(expand=c(0,0.02))+
              labs(x="",y="Species Richness",fill="",col="");p1

            ggsave("output/Figure3.png",p1,width=12,height=12,units="in",dpi=600)
      
