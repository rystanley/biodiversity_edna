library(vegan)
setwd("~/Documents/GitHub/biodiversity_edna/data/Vegan")

east <-read.table("FishCatchDiversityData3.txt", header=TRUE)
west <-read.table("West_Coast_seining_data2.txt", header=TRUE)

east_12S <-read.table("12S-final-fish-table-site_collapsed_East2.txt", header=TRUE)
west_12S <-read.table("12S-final-fish-table-site_collapsed_West2.txt", header=TRUE)

east_16S <-read.table("16S-final-fish-table-site-collapsed_East2.txt", header=TRUE)
west_16S <-read.table("16S-final-fish-table-site-collapsed_West2.txt", header=TRUE)

dim(east)
#10, 23
dim(west)
#9, 39
dim(east_12S)
#10, 43
dim(west_12S)
#9, 59
dim(east_16S)
#10, 41
dim(west_16S)
#9, 55

curve_east = specaccum(east[2:23], method = "random", 
                      permutations = 1000)
curve_east_12S = specaccum(east_12S[2:43], method = "random", 
                       permutations = 1000)
curve_east_16S = specaccum(east_16S[2:41], method = "random", 
                           permutations = 1000)

curve_west = specaccum(west[2:39], method = "random", 
                       permutations = 1000)
curve_west_12S = specaccum(west_12S[2:59], method = "random", 
                       permutations = 1000)
curve_west_16S = specaccum(west_16S[2:55], method = "random", 
                           permutations = 1000)

plot(curve_east_12S,col = "black",xlab="Sites", ylab="Number of species")
plot(curve_east, add = TRUE, col="red")
plot(curve_east_16S, add = TRUE, col="blue")
legend("bottomright", legend=c("seining", "12S", "16S"),
       col=c("red", "black", "blue"), lty=1, cex=0.8)
text(5, 5, "East coast")
     
plot(curve_west_12S,col = "black",xlab="Sites", ylab="Number of species")
plot(curve_west, add = TRUE, col="red")
plot(curve_west_16S, add = TRUE, col="blue")
legend("bottomright", legend=c("seining", "12S", "16S"),
       col=c("red", "black", "blue"), lty=1, cex=0.8)
text(5, 7, "West coast")

specpool(east[2:23]);specpool(east_12S[2:43]);specpool(east_16S[2:41])
specpool(west[2:39]);specpool(west_12S[2:59]);specpool(west_16S[2:55])

