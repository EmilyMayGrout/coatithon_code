#making plots for the matrix per day
#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level1.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)

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

#removing males and wildflower
xs <- xs[c(1:4,6,7,11:17,19,20,22),]
ys <- ys[c(1:4,6,7,11:17,19,20,22),]
coati_ids <- coati_ids[-c(5,8,9,10,18,21),]

R = 50

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

first_time <- c("1", "79", "157", "235", "313", "391", "469", "547", "625", "703", "781", "859", "937", "1015", "1093")
last_time <- c("78", "156", "234", "312", "390", "468", "546", "624", "702", "780", "858", "936", "1014", "1092", "1165")

times <- data.frame(first_time, last_time)

png(height = 2000, width = 2000, units = 'px', filename = paste0(plot_dir,'perday_allinds_level1.png'))
par(mfrow=c(4,4),mgp=c(3, 1, 0), mar=c(11,10,4,2))


for (m in 1:15){
  subgroup_data <- get_subgroup_data(xs[,times$first_time[m]:times$last_time[m]], ys[,times$first_time[m]:times$last_time[m]], R=R)
  ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
      sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
      ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
    }
  }
  diag(ff_net) <- NA
  #new_order <- c(13,1,9,10,3,4,12,2,11,7,5,16,8,15,6,14)
  new_order <- c(21,1,13,20,7,14,17,6,10,19,3,4,11,12,16,2,15,22,5,8,9,18)
  ffnet_reorder <- ff_net[new_order, new_order]
  
  visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
  mtext(paste("Day", m))

}

dev.off()








