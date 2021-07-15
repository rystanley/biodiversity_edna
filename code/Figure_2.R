## This code generates Figure 2

# Comparisons of fish taxa detected by eDNA and seining at the family and species level. For the lower panel, only taxa at the species level were considered. 

#load libraries -------------
library(VennDiagram)
library(gridExtra)
library(ggplot2)

#set up venn contrasts using the results outlined in Results section of manuscript 

grid.newpage() #set up clean plotting window so the plot assigned to the object can be inspected
species_com<-draw.pairwise.venn(area1 = 103, area2 = 52, cross.area = 45, category = c("eDNA-species", "seining-species"),lty = rep("solid", 2),
                                fill = c("light blue", "pink"), alpha = rep(0.7, 2), scaled = TRUE,cat.dist = c(0.040, 0.05),cat.pos = c(-155,155),
                                cex = 2.5, cat.cex = 2)

grid.newpage()
family_com<-draw.pairwise.venn(area1 = 44, area2 = 28, cross.area = 27, category = c("eDNA-family", "seining-family"),
                               fill = c("light blue", "pink"), alpha = rep(0.7, 2), cat.dist = c(0.03, 0.055),cat.pos = c(-155,155),
                               cex = 2.5, cat.cex = 2, scaled = TRUE)

#Combine plots into a single object and save -------------

Venn_com<-grid.arrange(grobTree(family_com),grobTree(species_com),ncol=1)
ggsave(plot=Venn_com, file="output/Figure_2.png", width = 7, height = 14,units="in",dpi=600)
