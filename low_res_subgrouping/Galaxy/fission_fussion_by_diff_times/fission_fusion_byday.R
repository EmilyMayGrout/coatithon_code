#this code is for making an association matrix for each day with the 1/10 gps for the whole day

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
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

first24 <- as.POSIXct('2021-12-26 11:00', tz = 'UTC')
last24 <-  as.POSIXct('2022-12-26 22:50', tz = 'UTC')

day1 <- ts[1:78]
day2 <- ts[79:156]
day3 <- ts[157:234]
day4 <- ts[235:312]
day5 <- ts[313:390]
day6 <- ts[391:468]
day7 <- ts[469:546]
day8 <- ts[547:624]
day9 <- ts[625:702]
day10 <- ts[703:780]
day11 <- ts[781:858]
day12 <- ts[859:936]
day13 <- ts[937:1014]
day14 <- ts[1015:1092]
day15 <- ts[1093:1170]
day16 <- ts[1171:1248]


#---------make matrix-------------------------------------------
png(height = 1050, width = 1050, units = 'px', filename = paste0(plot_dir,'subgroup_network_perday.png'))

par(mfrow=c(4,4), mar = c(5,5,5,5))

#day1:
subgroup_data <- get_subgroup_data(xs[,1:78], ys[,1:78], R=50)
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
mtext("Day 1")
#-----------------------------------------------------------------
#day2:
subgroup_data <- get_subgroup_data(xs[,79:156], ys[,79:156], R=50)
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
mtext("Day 2")
#-----------------------------------------------------------------
#day3:
subgroup_data <- get_subgroup_data(xs[,157:234], ys[,157:234], R=50)
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
mtext("Day 3")
#-----------------------------------------------------------------
#day4:
subgroup_data <- get_subgroup_data(xs[,235:312], ys[,235:312], R=50)
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
mtext("Day 4")
#-----------------------------------------------------------------
#day5:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 5")
#-----------------------------------------------------------------
#day6:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 6")
#-----------------------------------------------------------------
#day7:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 7")
#-----------------------------------------------------------------
#day8:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 8")
#-----------------------------------------------------------------
#day9:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 9")
#-----------------------------------------------------------------
#day10:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 10")
#-----------------------------------------------------------------
#day11:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 11")
#-----------------------------------------------------------------
#day12:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 12")
#-----------------------------------------------------------------
#day13:
subgroup_data <- get_subgroup_data(xs[,937:1014], ys[,937:1014], R=50)
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
mtext("Day 13")
#-----------------------------------------------------------------
#day14:
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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
mtext("Day 14")
#-----------------------------------------------------------------
#day15:
subgroup_data <- get_subgroup_data(xs[,1093:1170], ys[,1093:1170], R=50)
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
mtext("Day 15")
#-----------------------------------------------------------------
#day16:
subgroup_data <- get_subgroup_data(xs[,1171:1248], ys[,1171:1248], R=50)
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
mtext("Day 16")

dev.off()


#comparing morning vs afternoon
#find times for morning (6am-12pm) and afternoons (12:10pm-6pm)
#UTC times 11:00-17:00 morning
#17:10-23:00 evening
mornings <- ts[c(1:37,79:115, 157:193,236:271,313:349,391:427,469:505,547:583,625:661,
                 703:739,781:817,859:895, 937:973,1015:1051,1093:1129, 1171:1207)]

afternoons <- ts[c(38:78,116:156,194:235, 272:312, 350:390, 428:468, 506:546, 584:624,
                   662:702,740:780,818:858,896:936,974:1014,1052:1092,1130:1170,1208:1248)]

#these are the timestamps for mornings and afternoons which I made manually
mor1 <- ts[1:37]
aft1 <- ts[38:78]
mor2 <- ts[79:115]
aft2 <- ts[116:156]
mor3 <- ts[157:193]
aft3 <- ts[194:235]
mor4 <- ts[236:271]
aft4 <- ts[272:312]
mor5 <- ts[313:349]
aft5 <- ts[350:390]
mor6 <- ts[391:427]
aft6 <- ts[428:468]
mor7 <- ts[469:505]
aft7 <- ts[506:546]
mor8 <- ts[547:583]
aft8 <- ts[584:624]
mor9 <- ts[625:661]
aft9 <- ts[662:702]
mor10 <-ts[703:739]
aft10 <-ts[740:780]
mor11 <-ts[781:817]
aft11 <-ts[818:858]
mor12 <-ts[859:895]
aft12 <-ts[896:936]
mor13 <-ts[937:973]
aft13 <-ts[974:1014]
mor14 <-ts[1015:1051]
aft14 <-ts[1052:1092]
mor15 <-ts[1093:1129]
aft15 <-ts[1130:1170]
mor16 <-ts[1171:1207]
aft16 <-ts[1208:1248]

png(height = 500, width = 1000, units = 'px', filename = paste0(plot_dir,'subgroup_network_morning_afternoon.png'))

par(mfrow=c(1,2), mar = c(5,5,5,5))


#mornings
subgroup_data <- get_subgroup_data(xs[,c(1:37,79:115, 157:193,236:271,313:349,391:427,469:505,547:583,625:661,
                                         703:739,781:817,859:895, 937:973,1015:1051,1093:1129, 1171:1207)], 
                                   ys[,c(1:37,79:115, 157:193,236:271,313:349,391:427,469:505,547:583,625:661,
                                         703:739,781:817,859:895, 937:973,1015:1051,1093:1129, 1171:1207)], R=50)
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
mtext("Mornings")
#-----------------------------------------------------------
#afternoons
subgroup_data <- get_subgroup_data(xs[,c(38:78,116:156,194:235, 272:312, 350:390, 428:468, 506:546, 584:624,
                                         662:702,740:780,818:858,896:936,974:1014,1052:1092,1130:1170,1208:1248)], 
                                   ys[,c(38:78,116:156,194:235, 272:312, 350:390, 428:468, 506:546, 584:624,
                                         662:702,740:780,818:858,896:936,974:1014,1052:1092,1130:1170,1208:1248)], R=50)
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
mtext("Afternoons")

dev.off()

