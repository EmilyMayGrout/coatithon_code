
#this script is to build a social network for the times when the group is together
#need to find times when all individuals are together with and without singletons

#------DIRECTORIES AND FILES------

#EMILY DIRECTORIES
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'

#ARI DIRECTORIES
#data_dir <- "/Users/Ari/Dropbox/coatithon/processed/2022/"
#code_dir <- '/Users/Ari/Dropbox/code_ari/coatithon_code/'
#plot_dir <- '/Users/Ari/Dropbox/coatithon/results/'

gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData' 

#-----PARAMETERS-----
min_tracked <- 8 #the minimum number of individuals tracked at any point to consider a time point in this analysis
R <- 50 # the radius used to define subgroups
R_within <- 10 #the radius used to define subgroups within the larger group

#-----LIBRARIES AND SOURCED FUNCTIONS-----
library(fields)
library(viridis)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#get who is tracked at each time point
tracked <- !is.na(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(tracked)

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is R (default 50m)
subgroup_data <- get_subgroup_data(xs, ys, R)

#also get the subgroup data for a smaller radius, R_within (default = 10)
subgroup_data_within <- get_subgroup_data(xs, ys, R_within)

#run through time by time and check if each time is when the group is together
#want to ignore times when there are singletons
together <- c()

for(t in 1:(n_times-1)){
  
  #get subgroup membership now and later
  group_membership <- subgroup_data$ind_subgroup_membership[,t]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(group_membership))==0){
    next
  }
  
  #get number of subgroups in group_membership
  n_subgroups <- length(unique(group_membership[!is.na(group_membership)]))
  #get number of singletons in group
  n_singletons <- sum(table(group_membership)==1)
  
  #determine if this time step is when the full group is together
  if(n_subgroups==1 
  ){
    together <- c(together, t)
  }
  #if we have more than one group, but rest are singletons, still class group as together
  if(n_subgroups > 1 
     & (n_singletons + 1) == n_subgroups
  ){
    together <- c(together, t)
  }
}

#together is the indx of time points when the full group is together or when the full group is together and there's singletons


#NOTE: I (Ari) removed this bit because you can just index in to the original matrices at the correct time indices (in other words it gets filtered in the loop below) 
#filter the xs and ys to the indx of times when the group is together
# 
# xs_fullgroup <- xs[1:11, together]
# ys_fullgroup <- ys[1:11, together]
# 
# fullgroup_data <- get_subgroup_data(xs_fullgroup, ys_fullgroup, R=50)
# 
# n_inds <- nrow(xs_fullgroup)
# n_times <- ncol(xs_fullgroup)
# 
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)

i = 5
j = 6
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #only time steps where:
    #1) at least min_tracked individuals are tracked from the entire group
    #2) i and j are both tracked 
    #3) i and j are in the same subgroup (at a scale of R=50)
    usable_times <- which(n_tracked >= min_tracked 
                           & !is.na(subgroup_data$ind_subgroup_membership[i,]) 
                           & !is.na(subgroup_data$ind_subgroup_membership[j,])
                           & subgroup_data$ind_subgroup_membership[i,] == subgroup_data$ind_subgroup_membership[j,])
    
    #find the usable time steps where the group is also together
    usable_times_together <- intersect(together, usable_times)
    
    #getting subgroup id for individual i and j, only at times when they were together and that are usable according to definition above
    #Note that this uses the smaller radius, because we want to look at within-group structure
    sub_ids_i <- subgroup_data_within$ind_subgroup_membership[i,usable_times_together]
    sub_ids_j <- subgroup_data_within$ind_subgroup_membership[j,usable_times_together]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]

visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])


