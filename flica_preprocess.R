
#This script is for making a matrix with [ID, time, xy] as [1:11, 1:end of time, 1:2] for the high res periods

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_highres_level0.RData"
id_file <- 'coati_ids.RData'


#-------SETUP-------

library(fields)
library(viridis)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

dat <- array(NA, dim = c(nrow(xs), ncol(xs), 2))
dat[,,1] <- xs
dat[,,2] <- ys



#save(dat, file = paste0(outdir,'galaxy_xy_highres_array_flica.RData'))  

#check it opens correctly
#rm(list=ls())
#load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy_xy_highres_array_flica.RData")
#it opens correctly woo


#filter to fission time for Roi - chosen 27.12.21 12:10:00 - 12:40:00
#27.12.21 12:10:00
#ts[47401]
#27.12.21 12:40:00
#ts[49201]

#fissiontime <- ts[47401:49201]
#save(fissiontime, file = paste0(data_dir,'galaxy_xy_highres_fissiontimes_flica.RData'))  

#fission_event <- dat[,47401:49201,]
#save(fission_event, file = paste0(data_dir,'galaxy_xy_highres_fissionarray_flica.RData'))  


#now want to plot the times there are NA's in the dataset

#make data frame for first coati
Quasar <- as.data.frame(xs[1,])
colnames(Quasar)[colnames(Quasar) == "xs[1, ]"] <- "x"
Quasar$is_na[ is.na (Quasar$x) ] <- 0 
Quasar$is_na[ Quasar$x > 1] <- 1
Quasar$time <- ts

i=1
for (i in 1:11){
 
  xs_i<-data.frame("x"=xs[i,])
  print(str(xs_i))
  xs_i$is_na[ is.na (xs_i$x) ] <- i - 0.3 
  xs_i$is_na[ xs_i$x > 1] <- i
  xs_i$time <- ts
  assign(paste0("ID_",i),xs_i)
  
  
}

Sys.setenv(TZ='UTC')
day1_start <- which(ts == as.POSIXct("2021-12-24 11:00:00 UTC"))
day1_end <-   which(ts == as.POSIXct("2021-12-24 13:59:59 UTC"))
day2_start <- which(ts == as.POSIXct("2021-12-25 11:00:00 UTC"))
day2_end <-   which(ts == as.POSIXct("2021-12-25 13:59:59 UTC"))
day3_start <- which(ts == as.POSIXct("2021-12-26 11:00:00 UTC"))
day3_end <-   which(ts == as.POSIXct("2021-12-26 13:59:59 UTC"))
day4_start <- which(ts == as.POSIXct("2021-12-27 11:00:00 UTC"))
day4_end <-   which(ts == as.POSIXct("2021-12-27 13:59:59 UTC"))
day5_start <- which(ts == as.POSIXct("2021-12-28 11:00:00 UTC"))
day5_end <-   which(ts == as.POSIXct("2021-12-28 13:59:59 UTC"))
day6_start <- which(ts == as.POSIXct("2021-12-29 11:00:00 UTC"))
day6_end <-   which(ts == as.POSIXct("2021-12-29 13:59:59 UTC"))
day7_start <- which(ts == as.POSIXct("2021-12-30 11:00:00 UTC"))
day7_end <-   which(ts == as.POSIXct("2021-12-30 13:59:59 UTC"))
day8_start <- which(ts == as.POSIXct("2021-12-31 11:00:00 UTC"))
day8_end <-   which(ts == as.POSIXct("2021-12-31 13:59:59 UTC"))
day9_start <- which(ts == as.POSIXct("2022-01-01 11:00:00 UTC"))
day9_end <-   which(ts == as.POSIXct("2022-01-01 13:59:59 UTC"))
day10_start <- which(ts == as.POSIXct("2022-01-02 11:00:00 UTC"))
day10_end <-   which(ts == as.POSIXct("2022-01-02 13:59:59 UTC"))
day11_start <- which(ts == as.POSIXct("2022-01-03 11:00:00 UTC"))
day11_end <-   which(ts == as.POSIXct("2022-01-03 13:59:59 UTC"))
day12_start <- which(ts == as.POSIXct("2022-01-04 11:00:00 UTC"))
day12_end <-   which(ts == as.POSIXct("2022-01-04 13:59:59 UTC"))
day13_start <- which(ts == as.POSIXct("2022-01-05 11:00:00 UTC"))
day13_end <-   which(ts == as.POSIXct("2022-01-05 13:59:59 UTC"))
day14_start <- which(ts == as.POSIXct("2022-01-06 11:00:00 UTC"))
day14_end <-   which(ts == as.POSIXct("2022-01-06 13:59:59 UTC"))
day15_start <- which(ts == as.POSIXct("2022-01-07 11:00:00 UTC"))
day15_end <-   which(ts == as.POSIXct("2022-01-07 13:59:59 UTC"))
day16_start <- which(ts == as.POSIXct("2022-01-08 11:00:00 UTC"))
day16_end <-   which(ts == as.POSIXct("2022-01-08 13:59:59 UTC"))


str(ts)
day1 <- day1_start:day1_end
day2 <- day2_start:day2_end
day3 <- day3_start:day3_end
day4 <- day4_start:day4_end
day5 <- day5_start:day5_end
day6 <- day6_start:day6_end
day7 <- day7_start:day7_end
day8 <- day8_start:day8_end
day9 <- day9_start:day9_end
day10 <- day10_start:day10_end
day11 <- day11_start:day11_end
day12 <- day12_start:day12_end
day13 <- day13_start:day13_end
day14 <- day14_start:day14_end
day15 <- day15_start:day15_end
day16 <- day16_start:day16_end

test<-list(day1,day2,day3, day4, day5, day6, day7, day8,day9, day10, day11,
          day12, day13, day14, day15, day16)
i=1

png(height = 1080, width = 1080, units = 'px', filename = paste0(plot_dir,'na_highres.png'))
par(mfrow=c(4,4), mar = c(2,1,1,1))
j<-0
for(i in test){
  j<-j+1
#day <- paste0(day, i)

cex = 0.5

#xlim = c(as.POSIXct("2022-01-08 11:00:00 UTC"),as.POSIXct("2022-01-08 13:59:59 UTC") 
#choose times to plot
plot(ID_1[rownames(ID_1)%in%i,]$time, ID_1[rownames(ID_1)%in%i,]$is_na, col= "#9e0142", ylim = c(0,11),  pch = 20, cex = cex, xlab = "Time", ylab = "ID", main = paste0("Day ",j))
points(ID_2[rownames(ID_2)%in%i,]$time, ID_2[rownames(ID_2)%in%i,]$is_na, col= "#d53e4f", pch = 20, cex = cex)
points(ID_3[rownames(ID_3)%in%i,]$time, ID_3[rownames(ID_3)%in%i,]$is_na, col= "#f46d43", pch = 20, cex = cex)
points(ID_4[rownames(ID_4)%in%i,]$time, ID_4[rownames(ID_4)%in%i,]$is_na, col= "#fdae61", pch = 20, cex = cex)
points(ID_5[rownames(ID_5)%in%i,]$time, ID_5[rownames(ID_5)%in%i,]$is_na, col= "#fee08b", pch = 20, cex = cex)
points(ID_6[rownames(ID_6)%in%i,]$time, ID_6[rownames(ID_6)%in%i,]$is_na, col= "#ffffbf",pch = 20, cex = cex)
points(ID_7[rownames(ID_7)%in%i,]$time, ID_7[rownames(ID_7)%in%i,]$is_na, col= "#e6f598",pch = 20, cex = cex)
points(ID_8[rownames(ID_8)%in%i,]$time, ID_8[rownames(ID_8)%in%i,]$is_na, col= "#abdda4", pch = 20, cex = cex)
points(ID_9[rownames(ID_9)%in%i,]$time, ID_9[rownames(ID_9)%in%i,]$is_na, col= "#66c2a5",pch = 20, cex = cex)
points(ID_10[rownames(ID_10)%in%i,]$time, ID_10[rownames(ID_10)%in%i,]$is_na, col= "#3288bd", pch = 20, cex = cex)
points(ID_11[rownames(ID_11)%in%i,]$time, ID_11[rownames(ID_11)%in%i,]$is_na, col= "#5e4fa2",pch = 20, cex = cex)


}

dev.off()

#need to add other group members to this, have their is_na value from 1 to 11 then put them onto same plot 


