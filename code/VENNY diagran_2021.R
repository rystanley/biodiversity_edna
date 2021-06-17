grid.newpage()
setwd("~/GRDI_R/Venn diagram")
#install.packages("VennDiagram")
library(VennDiagram)
########################################
species_com<-draw.pairwise.venn(area1 = 103, area2 = 52, cross.area = 45, category = c("eDNA-species", "seining-species"),lty = rep("solid", 2),
                   fill = c("light blue", "pink"), alpha = rep(0.7, 2), scaled = TRUE,cat.dist = c(0.040, 0.05),cat.pos = c(-155,155),
                   cex = 2.5, cat.cex = 2)

family_com<-draw.pairwise.venn(area1 = 44, area2 = 28, cross.area = 27, category = c("eDNA-family", "seining-family"),
                              fill = c("light blue", "pink"), alpha = rep(0.7, 2), cat.dist = c(0.03, 0.055),cat.pos = c(-155,155),
                              cex = 2.5, cat.cex = 2, scaled = TRUE)
#library(grid)
#library(gridBase)
#install.packages("gridBase")

library(ggplot2)
library(gridExtra)
#Venn_com<-grid.arrange(grobTree(species_com),grobTree(genus_com),grobTree(family_com),ncol=1)

Venn_com<-grid.arrange(grobTree(family_com),grobTree(species_com),ncol=1)
ggsave(plot=Venn_com, file="Venn_comparison_2021May14.pdf", width = 7, height = 14)