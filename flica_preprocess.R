
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
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy_xy_highres_array_flica.RData")
#it opens correctly woo


#filter to fission time for Roi - chosen 27.12.21 12:10:00 - 12:40:00
#27.12.21 12:10:00
ts[47401]
#27.12.21 12:40:00
ts[49201]

fissiontime <- ts[47401:49201]
#save(fissiontime, file = paste0(data_dir,'galaxy_xy_highres_fissiontimes_flica.RData'))  

fission_event <- dat[,47401:49201,]
#save(fission_event, file = paste0(data_dir,'galaxy_xy_highres_fissionarray_flica.RData'))  


#now want to plot the times there are NA's in the dataset

#make data frame for first coati
Quasar <- as.data.frame(xs[1,])
colnames(Quasar)[colnames(Quasar) == "xs[1, ]"] <- "x"
Quasar$is_na[ is.na (Quasar$x) ] <- 0 
Quasar$is_na[ Quasar$x > 1] <- 1
Quasar$time <- ts

Estrella <- as.data.frame(xs[2,])
colnames(Estrella)[colnames(Estrella) == "xs[2, ]"] <- "x"
Estrella$is_na[ is.na (Estrella$x) ] <- 0 
Estrella$is_na[ Estrella$x > 1] <- 2
Estrella$time <- ts

#perhaps this graph can be done as a forloop through the ID's

#choose times to plot
plot(Quasar$time[27401:49201], Quasar$is_na[27401:49201], ylim = c(0,3))
points(Estrella$time[27401:49201], Estrella$is_na[27401:49201])

#need to add other group members to this, have their is_na value from 1 to 11 then put them onto same plot 


