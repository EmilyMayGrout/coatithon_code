#This script is for looking at which individuals sleep close together
#need to fix bugs in a couple of the graphs



#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_night_level0.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)

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

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


#who is in the same sub-group in total data

subgroup_data <- get_subgroup_data(xs, ys, R=50)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network_night.png'))
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
dev.off()



#do same for each day:
night1 <- ts[1:66] #go from 00:00 to 10:50 which is 7pm to 5.50am 
night2 <- ts[67:132]
night3 <- ts[133:198]
night4 <- ts[199:264]
night5 <- ts[265:330]
night6 <- ts[331:396]
night7 <- ts[397:462]
night8 <- ts[463:528]
night9 <- ts[529:594]
night10<- ts[595:660]
night11<- ts[661:726]
night12<- ts[727:792]
night13<- ts[793:858]
night14<- ts[859:924]
night15<- ts[925:990]
night16<- ts[991:1056]

par(mfrow=c(4,4), mar = c(5,5,5,5))

#day1:
subgroup_data <- get_subgroup_data(xs[,1:66], ys[,1:66], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 1")
#-----------------------------------------------------------------
#day2:
subgroup_data <- get_subgroup_data(xs[,67:132], ys[,67:132], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 2")
#-----------------------------------------------------------------
#day3:
subgroup_data <- get_subgroup_data(xs[,133:198], ys[,133:198], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 3")
#-----------------------------------------------------------------
#day4:
subgroup_data <- get_subgroup_data(xs[,199:264], ys[,199:264], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 4")
#-----------------------------------------------------------------
#day5:
subgroup_data <- get_subgroup_data(xs[,265:330], ys[,265:330], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 5")
#-----------------------------------------------------------------
#day6:
subgroup_data <- get_subgroup_data(xs[,331:396], ys[,331:396], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 6")
#-----------------------------------------------------------------
#day7:
subgroup_data <- get_subgroup_data(xs[,397:462], ys[,397:462], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 7")
#-----------------------------------------------------------------
#day8:
subgroup_data <- get_subgroup_data(xs[,463:528], ys[,463:528], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 8")
#-----------------------------------------------------------------
#day9:
subgroup_data <- get_subgroup_data(xs[,529:594], ys[,529:594], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 9")
#-----------------------------------------------------------------
#day10:
#subgroup_data <- get_subgroup_data(xs[,595:660], ys[,595:660], R=50)
#ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
#for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
#diag(ff_net) <- NA
#new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
#ffnet_reorder <- ff_net[new_order, new_order]
#visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
#mtext("Night 10") #dodge
#-----------------------------------------------------------------
#day11:
subgroup_data <- get_subgroup_data(xs[,661:726], ys[,661:726], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 11") 
#-----------------------------------------------------------------
#day12:
subgroup_data <- get_subgroup_data(xs[,727:792], ys[,727:792], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 12")
#-----------------------------------------------------------------
#day13:
subgroup_data <- get_subgroup_data(xs[,793:858], ys[,793:858], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 13")
#-----------------------------------------------------------------
#day14:
#subgroup_data <- get_subgroup_data(xs[,859:924], ys[,859:924], R=50)
#ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
#for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
#diag(ff_net) <- NA
#new_order <- c(1,11,4,10,2,3,6,7,8,9,5)
#ffnet_reorder <- ff_net[new_order, new_order]
#visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
#mtext("Night 14") #dodge
#-----------------------------------------------------------------
#day15:
subgroup_data <- get_subgroup_data(xs[,925:990], ys[,925:990], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 15")
#-----------------------------------------------------------------
#day16:
subgroup_data <- get_subgroup_data(xs[,991:1056], ys[,991:1056], R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Night 16")

dev.off()

