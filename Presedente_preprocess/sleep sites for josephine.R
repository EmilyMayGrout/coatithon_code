#getting sleep site info for Josephine

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level0_sleep.RData"
id_file <- 'presedente_coati_ids.RData'

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#filter to Gendry, Kenyatta, Lula, Mandela, Sam, Truss and Ardern
inds <- c(1,5,8,9,10,18,20)
Xs <- xs[inds,]
Ys <- ys[inds,]
coati_ids <- coati_ids[inds,]

xs <- data.frame(t(Xs))
name_x <- paste(coati_ids$name, "_x")
colnames(xs)<- c(name_x)

ys <- data.frame(t(Ys))
name_y <- paste(coati_ids$name, "_y")
colnames(ys) <- c(name_y)

df <- cbind(xs, ys, ts)

df$na_count <- rowSums(is.na(df))
#removing rows where all inds have NA
df <- df[!df$na_count == 14,]

save(df, file = "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/sleep_site_UTM.RDa")
load("C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/sleep_site_UTM.RDa")









