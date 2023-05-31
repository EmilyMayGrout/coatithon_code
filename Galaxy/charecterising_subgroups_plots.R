#plots for the low res fission fusion paper galaxy

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
gps_file <- "galaxy_xy_10min_level1.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 50

#-------SETUP-------
#making plots for characterizing the subgroups

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

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#make 3 plots next to each other, first is hist of number of subgroups at 50m, next it distribution of individuals in each group size with 2 subgroups, last is at 3 subgroups

png(height = 400, width = 1280, units = 'px', filename = paste0(plot_dir,'50m_charecterisations.png'))

par(mfrow=c(1,3), mar = c(6,5,2,1)) #(bottom, left, top, right)
subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs],main = "50m radii", xlab =  'Number of subgroups', col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 6))

subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), xaxp = c(1,11, 10))
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), xaxp = c(1,11, 10))

dev.off()


#------------------------------------------------------------------
#different days for the matrix

png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_28thDec.png'))
par(mfrow=c(1,1), mar = c(6,7,3,3))#(bottom, left, top, right)

subgroup_data <- get_subgroup_data(xs[,313:390], ys[,313:390], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("28.12.2021", cex = 3)
dev.off()

#-----------------------------------------------------------
png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_29thDec.png'))
subgroup_data <- get_subgroup_data(xs[,391:468], ys[,391:468], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5, 1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("29.12.21", cex = 3)
dev.off()
#-------------------------------------------------------------

png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_30thDec.png'))

subgroup_data <- get_subgroup_data(xs[,469:546], ys[,469:546], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5, 1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("30.12.21", cex = 3)
dev.off()
#----------------------------------------------------------
png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_31stDec.png'))


subgroup_data <- get_subgroup_data(xs[,547:624], ys[,547:624], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5, 1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("31.12.21", cex = 3)
dev.off()
#----------------------------------------------------------
png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_1stJan.png'))

subgroup_data <- get_subgroup_data(xs[,625:702], ys[,625:702], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5, 1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("01.01.22", cex = 3)
dev.off()


#-------------------------------------------------------------------
png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_2ndJan.png'))
par(mfrow=c(1,1), mar = c(6,7,3,3))#(bottom, left, top, right)

subgroup_data <- get_subgroup_data(xs[,703:780], ys[,703:780], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("02.01.2022", cex = 3)

dev.off()
#----------------------------------------------------------------------
png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_3rdJan.png'))
par(mfrow=c(1,1), mar = c(6,7,3,3))#(bottom, left, top, right)

subgroup_data <- get_subgroup_data(xs[,781:858], ys[,781:858], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5, 1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("03.01.2022", cex = 3)

dev.off()


#vertical random days


png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_4thJan.png'))
par(mfrow=c(1,1), mar = c(6,7,3,3))#(bottom, left, top, right)


subgroup_data <- get_subgroup_data(xs[,859:936], ys[,859:936], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("04.01.2022",  cex = 3)

dev.off()
#--------------------------------------------------------------

png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'rando_matrix_6thJan.png'))
par(mfrow=c(1,1), mar = c(6,7,3,3))#(bottom, left, top, right)

subgroup_data <- get_subgroup_data(xs[,1015:1092], ys[,1015:1092], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("06.01.2022", cex = 3)

dev.off()
