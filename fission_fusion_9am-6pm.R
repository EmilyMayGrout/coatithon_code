#this script is to build a social network with the low res data from 9am - 6pm


#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
#looking at data from 9am-6pm
gps_file <- "galaxy_xy_10min_9am_level0.RData"
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

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)



#-----------plot 1: social network for high res when group is split-----------------------
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


png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network_9am-6pm.png'))

visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
dev.off()


#-----------plot 2: social network for high res data when group is together --------------
R = 50
subgroup_data <- get_subgroup_data(xs, ys, R)

#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)
#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs[,full_group_index]
subset_y <- ys[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 10)


new_order <- c(1,11,4,10,2,3,6,7,8,9,5)

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'withingroup_network_withgus_9am-6pm.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids[new_order,])
dev.off()
