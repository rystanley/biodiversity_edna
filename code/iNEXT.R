#New code has been generated


setwd("~/Documents/GitHub/biodiversity_edna/data/iNEXT")

#install.packages('devtools')
#library(devtools)
#install_github('JohnsonHsieh/iNEXT')

#devtools::install_github("vmikk/metagMisc")
#library(metagMisc)

library(iNEXT)
library(ggplot2)
library(gridExtra)
main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                     axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    text=element_text(family="sans"))

east1 <-read.table("data.iNEXT.east.txt",header=T)
east.inc <- iNEXT(east1, q=0, datatype="incidence_freq", endpoint=30, knots = 30, se = TRUE, conf = 0.95,
                 nboot = 1000)
east.plot <- ggiNEXT(east.inc, type=1) + theme_bw(base_size = 12) +
  theme(legend.title=element_blank(), text=element_text(size=12))+ 
  main_theme + xlab("Number of sampling sites") + ylab("Species richness") +
 geom_text(x=5, y=80, label="Atlantic samples", color="black")

east.plot 


west1 <-read.table("data.iNEXT.west.txt",header=T)
west.inc <- iNEXT(west1, q=0, datatype="incidence_freq", endpoint=30, knots = 30, se = TRUE, conf = 0.95,
                  nboot = 1000)
west.plot <- ggiNEXT(west.inc, type=1) + theme_bw(base_size = 12) + 
  theme(legend.title=element_blank(), text=element_text(size=12)) + 
  main_theme + xlab("Number of sampling sites") + ylab("Species richness") +
  geom_text(x=5, y=90, label="Pacific samples", color="black")

west.plot 

plots <-arrangeGrob(east.plot,west.plot, ncol=1)
ggsave(plot=plots,file="species richness extrapolation 2021Apr16.pdf", width=12, height=14)



