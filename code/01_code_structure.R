## Code can all go in this folder. I would advise using headers to show what is happening in each code chunk
## All data should be stored and accessed through the common 'data' folder and outputs should all be directed
## to the 'output' folder. This way it will work seamlessly on each of our computers. 

#for example

#load libraries
library(dplyr)
library(ggplot2)
library(viridis)

#load data
dat <- read.csv("data/eDNA-Seine_Data_AugSep2019-2.csv")

#make plot of species by count

plotdata <- dat%>%
            rename(count=3)%>%#the name of the third column has a '-' which generally is not a good idea because it is a wildcard in R
            group_by(Species)%>%
            summarise(count=sum(count))%>%
            ungroup()%>%
            arrange(count)%>% #arrange in descending order
            mutate(Species=factor(Species,levels=Species)) #make a factor according to that order


p1 <- ggplot(data=plotdata)+
  geom_bar(aes(x=Species,y=count),stat="identity")+
  coord_flip()+
  theme_bw()+
  scale_y_log10(expand=c(0,0))
 

#save plot to the 'outputs folder' using a common png format. Note that these will not be tracked in Github
#becuase they are really changing in a way that git can track (think track changes in word). I will direct .gitignore to ignore all 
#pngs 

ggsave("output/testplot.png",p1,height=8,width=4,units="in",dpi=600)
