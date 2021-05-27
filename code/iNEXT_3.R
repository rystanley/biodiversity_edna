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

east1 <-read.table("data.iNEXT.east.txt",header=T)
east.inc <- iNEXT(east1, q=0, datatype="incidence_freq", endpoint=20, knots = 20, se = TRUE, conf = 0.95,
                 nboot = 1000)

east.plot <- ggiNEXT(east.inc, type=1) + theme_bw(base_size = 14) + 
  theme(legend.position="none") +
  theme(panel.grid=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"))+ 
  xlab("Number of sampling sites") + ylab("Species richness") +
 geom_text(x=5, y=75, label="(a) Atlantic samples", color="grey30")

east.plot 

east.plot3 <- ggiNEXT(east.inc, type=3) + theme_bw(base_size = 14) + xlim(c(0.8,1)) +
  theme(panel.grid=element_blank(),panel.border = element_blank(),axis.line = element_line(colour = "black"),legend.title=element_blank())+ 
  xlab("Sample coverage") + ylab("Species richness") +
  geom_text(x=0.9, y=75, label="(b) Atlantic samples", color="grey30")
east.plot3

west1 <-read.table("data.iNEXT.west.txt",header=T)
west.inc <- iNEXT(west1, q=0, datatype="incidence_freq", endpoint=18, knots = 18, se = TRUE, conf = 0.95,
                  nboot = 1000)
west.plot <- ggiNEXT(west.inc, type=1) + theme_bw(base_size = 14) + 
  theme(legend.position="none") +
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



