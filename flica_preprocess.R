
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



save(dat, file = paste0(outdir,'galaxy_xy_highres_array_flica.RData'))  

#check it opens correctly

rm(list=ls())
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy_xy_highres_array_flica.RData")
#it opens correctly woo


