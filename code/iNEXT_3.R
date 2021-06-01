#New code has been generated

# the intuitive part of working with github and RStudio projects is that all users start in the same root directory. 
# so there is no need to actually set a working directory because the paths will work on whichever computer it is on provided the project
# is initialized. In the first version of the code, for example you are writing plots to the data folder, whereas in the file structure those should really
# go to the outputs so that raw data cannot be accidentally saved over. This is the 

#setwd("~/Documents/GitHub/biodiversity_edna/data/iNEXT") 


#load libaries -------------
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
  #Original method using grob objects - XP 
      east.plot <- ggiNEXT(east.inc, type=1) + theme_bw(base_size = 14) + 
        theme(legend.position="none") +
        theme(panel.grid=element_blank(),
              #panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              #axis.text.x = element_blank(),
              axis.title.x = element_blank())+ 
        labs(x="",y="Species richness")+
       geom_text(x=5, y=75, label="(a) Atlantic samples", color="grey30")
      
      east.plot 
      
      east.plot3 <- ggiNEXT(east.inc, type=3) + theme_bw(base_size = 14) + xlim(c(0.8,1)) +
        theme(panel.grid=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),legend.title=element_blank())+ 
        xlab("Sample coverage") + ylab("Species richness") +
        geom_text(x=0.9, y=75, label="(b) Atlantic samples", color="grey30")
      
      east.plot3+facet_wrap("Atlantic"~.,strip.position = "right")
      
      west1 <-read.table("data/iNEXT/data.iNEXT.west.txt",header=T)
      west.inc <- iNEXT(west1, q=0, datatype="incidence_freq", endpoint=18, knots = 18, se = TRUE, conf = 0.95,
                        nboot = 1000)
      west.plot <- ggiNEXT(west.inc, type=1) + theme_bw(base_size = 14) + 
        theme(legend.position="none") 
        theme(panel.grid=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black")) + 
        xlab("Number of sampling sites") + ylab("Species richness") +
        geom_text(x=5, y=90, label="(c) Pacific samples", color="grey30")
      
      west.plot 
      
      west.plot3 <- ggiNEXT(west.inc, type=3) + theme_bw(base_size = 14) + xlim(c(0.8,1)) +
        theme(panel.grid=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),legend.title=element_blank()) + 
        xlab("Sample coverage") + ylab("Species richness") +
        geom_text(x=0.9, y=90, label="(d) Pacific samples", color="grey30")
      
      west.plot3 
      
      plots <-arrangeGrob(east.plot,east.plot3, west.plot, west.plot3, ncol=2)
      ggsave(plot=plots,file="species richness extrapolation 2021May27.pdf", width=12, height=12)

      #iNEXT modified facet method RS
      
          #assemble data from the iNEXT function using the handy formatting function 'fortify' modified by the iNEXT package
            plotdata <- rbind(fortify.iNEXT(east.inc,type=1),fortify.iNEXT(east.inc,type=3),
                              fortify.iNEXT(west.inc,type=1),fortify.iNEXT(west.inc,type=3))%>%
              mutate(region=ifelse(grepl("east",site),"Atlantic","Pacific"),
                     type=gsub("east_","",site),
                     type=gsub("west_","",type),
                     type=gsub("seining","Seining",type),
                     type=factor(type,levels=c("Seining","12S","16S","eDNA")),
                     denominator = ifelse(plottype==1,"Number of sampling sites","Sample coverage"),
                     region=factor(region,levels=c("Pacific","Atlantic")))%>%
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

            ggsave("output/RichnessAnalysis.png",p1,width=12,height=12,units="in",dpi=600)
      
          #Diversity plots      
            p2 <- ggplot()+
              geom_point(data=assym_div,aes(x=type,y=Estimator),size=2)+
              geom_errorbar(data=assym_div,aes(x=type,y=Estimator,ymin=LCL,ymax=UCL),width=0.25)+
              facet_grid(region~Diversity)+
              theme_bw()+
              labs(#y="Asymptotic diversity estimate ± 95% CI", #note if the code shows this plus minus as jibberish re-insert using alt+241
                y="Asymptotic diversity estimate", 
                  x="")+
              theme(strip.background = element_rect(, colour = "black", fill = "white"),
                    strip.text.x = element_text(colour = "black",size=14), 
                    strip.text.y = element_text(colour = "black",size=14),
                    axis.text = element_text(colour = "black",size=12),
                    axis.title = element_text(colour = "black",size=12));p2
            
            ggsave("output/AsymptoticDiversityEstimates.png",p2,width=12,height=9,units="in",dpi=600)
         