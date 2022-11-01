#analyse split events


data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

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

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, 50)

# this code is to find the change points when the group splits into smaller units (include all splits)

# #find indexes of subgroup_data where there is 1 group
# one_group_idxs <- which(subgroup_data$n_subgroups==1)
# split_idxs <- c()
# for(i in 1:length(one_group_idxs)){
#   #get index where whole group was together
#   current_idx <- one_group_idxs[i]
#   #get number of subgroups in next index
#   n_subs_next <- subgroup_data$n_subgroups[one_group_idxs[i]+1]
#   
#   #check if there is not an NA, otherwise ignore
#   if(!is.na(n_subs_next)){
#     
#     #check if the subgroups in the next step are > 1
#     if(n_subs_next > 1){
#       
#       #if so, add to the list of split_idxs
#       split_idxs <- c(split_idxs, one_group_idxs[i])
#     }
#   }
# }
# 
# 

#run through time by time and check if each time is a split
splits <- c()
for(t in 1:(n_times-1)){
  
  #get subgroup membership now and later
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(subgroups_now))==0 | sum(!is.na(subgroups_later))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
  n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(subgroups_now)==1)
  singletons_later <- sum(table(subgroups_later)==1)
  
  #determine if this time step is a split
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if(n_subgroups_now==1 
     & n_subgroups_later >1
     & singletons_later==0
     ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & (singletons_now+1) == n_subgroups_now
     & n_subgroups_later > n_subgroups_now
     & singletons_now == singletons_later
  ){
    splits <- c(splits, t)
  }
  
  
}

