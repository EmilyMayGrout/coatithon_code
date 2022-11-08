#this script is to make a dataframe for every time step for which individual is in which subgroup, whether its in the full group and whether all individuals were tracked at that time. Getting speed for each individual at each time step.

#---PARAMS----
R <- 50
dt <- 10 #time interval between points (=10 for coati low res). 10 is metres per minute.

#---DIRECTORIES----

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

#-----LIBRARIES-----

library(fields)
library(viridis)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')


#----LOAD DATA----
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

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, R)

#------------------------

#make a dataframe with time and individual

speed_df <- data.frame(t = rep(1:(n_times-1), n_inds), UTC_time = rep(ts[1:(n_times-1)], n_inds) ,ind = rep(1:n_inds, each = (n_times-1)), n_tracked = rep(n_tracked[1:(n_times-1)], n_inds), subgroup_size = NA, split = NA, speed = NA)

#loop through to calculate subgroup size, whether there is a split, and ind speed
i=1
for (i in 1:nrow(speed_df)){
  
  #get the current individual and time for that row
  time <- speed_df$t[i]
  ind <- speed_df$ind[i]
  
  #get current subgroups for all individuals at that time
  sub_data <- subgroup_data$ind_subgroup_membership[,time]
  
  #get current subgroup for the focal individual at that time
  ind_subgroup <- subgroup_data$ind_subgroup_membership[ind,time]
  
  #if that individual is not tracked, leave the NA and go to next row
  if(is.na(ind_subgroup)){
    next
  }
  
  #get subgroup size for the focal individual
  subgroup_size <- sum(sub_data == ind_subgroup, na.rm=T)
  
  #store subgroup size in df
  speed_df$subgroup_size[i] <- subgroup_size
  
  #determine whether group is "split" i.e. has at least 2 subgroups with at least 2 members each
  sub_counts <- subgroup_data$subgroup_counts[,time]
  n_subgroups_with_at_least_2_members <- sum(sub_counts >= 2, na.rm=T)
  split <- n_subgroups_with_at_least_2_members >= 2
  
  #store the split in the dataframe
  speed_df$split[i] <- split
  
  #calculate speed with the xs and ys
  ind_xs_current <- xs[ind, time]
  ind_ys_current <- ys[ind, time]
  
  next_time <- time + 1
  ind_xs_next <- xs[ind, next_time]
  ind_ys_next <- ys[ind, next_time]
  
  dx <- ind_xs_next - ind_xs_current
  dy <- ind_ys_next - ind_ys_current
  
  #getting 
  speed <- (sqrt((dx)^2 + (dy)^2))/dt
  
  speed_df$speed[i] <- speed
  
}







