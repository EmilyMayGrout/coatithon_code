#this script is to find the call rates for each group size. To do this:
#1. we find consecutive periods when the group doesn't change
#2. count the number of individuals in each of those groupings
#3. get the group size for each of the grouping options

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'
gps_file <- "galaxy_xy_highres_level2.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#list of Rs
R <- 50

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vioplot)
library(plotly)


#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

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

#----------------------------------------------------------------------
#getting subgroup data
subgroup_data <- get_subgroup_data(xs, ys, R)



#----------------------------------------------------------------------
#get the call files 
datadir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
callfile <- 'all_data_hms_all_ml_synched.csv' #made in coati_synch
id_file <- 'galaxy_coati_ids.RData'

#LOAD DATA
setwd(datadir)
#calls <- read.csv(callfile, header=T, sep = ',') 
#calls <- calls[,-1]
load(id_file)


#go through each of the synched times in call dataframe and get the number of individuals in the subgroup of that group member

i = "G9463"
j = 1

# #initialise column to add the count of individuals in each subgroup
# calls$count <- NA
# 
# #this for loop takes 5 hours to run
# # Loop through each unique individual ID
# for (i in unique(calls$id)) {
# 
#   # Subset the data for the current individual
#   ind_i <- calls[calls$id == i, ]
# 
#   # Remove the 'G' to get the collar ID - to find the coati_ids index
#   ind_i$id <- gsub('G', '', ind_i$id)
#   id_indx <- which(coati_ids$tag_id == unique(ind_i$id))
# 
#   # Loop through each of the times in ind_i$datetime_synch to find the number of individuals in the subgroup of that individual
#   for (j in 1:nrow(ind_i)) {
# 
#     print(j)
# 
#     time <- as.POSIXct(ind_i$datetime_synch[j], format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
#     time <- round_date(time, "second")
# 
#     # Find the ts index for this time
#     t <- which(ts == time)
# 
#     # Skip to the next iteration if t is empty
#     if (length(t) == 0) {
#       next
#     }
# 
#     # Get the information of the subgroup membership for each individual for each time
#     subgroup_info_j <- subgroup_data$ind_subgroup_membership[, t]
#     subgroup_id_for_j <- subgroup_data$ind_subgroup_membership[id_indx, t]
# 
#     # If there is any NA in subgroup_info_j, skip to the next row
#     if (any(is.na(subgroup_info_j))) {
#       next
#     }
#     # Count the number of values equal to the subgroup_id_for_j, excluding NAs
#     count <- sum(subgroup_info_j == subgroup_id_for_j, na.rm = TRUE)
# 
#     # Assign the count to the appropriate row in the calls dataframe
#     calls$count[calls$id == i & calls$datetime_synch == ind_i$datetime_synch[j]] <- count
# 
#   }
# }

#save(calls, file = paste0(datadir, "/all_data_hms_all_ml_synched_subgrp_size_noNas.RData"))
#write.csv(calls, file = paste0(datadir, "/all_data_hms_all_ml_synched_subgrp_size_noNas.csv"))

calls <- read.csv(paste0(datadir, "/all_data_hms_all_ml_synched_subgrp_size_noNas.csv"))
#write.csv(calls, file = paste0(datadir, "/all_data_hms_all_ml_synched_subgrp_size_noNas.csv"))
calls <- calls[,-c(1:3)]


#combining the contact calls 
calls$calltype <- NA
calls$calltype[calls$label == "chirpgr"] <- "contact call"
calls$calltype[calls$label == "chirp grunt"] <- "contact call"
calls$calltype[calls$label == "chirp click"] <- "contact call"
calls$calltype[calls$label == "click grunt"] <- "contact call"
calls$calltype[calls$label == "click"] <- "contact call"
calls$calltype[calls$label == "chirp"] <- "contact call"

#combine aggressive calls
calls$calltype[calls$label == "chitter"] <- "aggression call"
#calls$calltype[calls$label == "squeal"] <- "aggression call"
#calls$calltype[calls$label == "squeal chitter"] <- "aggression call"
#calls$calltype[calls$label == "squeal chitter x"] <- "aggression call"
#calls$calltype[calls$label == "squeal chitters"] <- "aggression call"
#calls$calltype[calls$label == "low squeal"] <- "aggression call"
calls$calltype[calls$label == "chitter x"] <- "aggression call"
#calls$calltype[calls$label == "squeal chittering"] <- "aggression call"

# Ensure datetime_synch is in POSIXct format
calls$datetime_synch <- as.POSIXct(calls$datetime_synch, format="%Y-%m-%d %H:%M:%S", tz = "UTC")

# Create 2-minute bins
calls$time_bin <- cut(calls$datetime_synch, breaks="2 mins")

# Convert time_bin to POSIXct for start and end time extraction
calls$time_bin <- as.POSIXct(calls$time_bin, tz = "UTC")

# Define a function to calculate the mode
calculate_mode <- function(x) {
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Generate the complete grid of all combinations of id, date, and time_bin
complete_grid <- expand_grid(
  id = unique(calls$id),
  date = unique(calls$date),
  time_bin = unique(calls$time_bin)
)


# Summarize the data to get mean and mode subgroup sizes
subgroup_size_summary <- calls %>%
  group_by(id, date, time_bin) %>%
  summarise(
    mean_subgroup_size = mean(count, na.rm = TRUE),
    mode_subgroup_size = calculate_mode(count[!is.na(count)]),
    .groups = 'drop'
  )

# Merge the complete grid with the subgroup size data
results_complete <- complete_grid %>%
  left_join(subgroup_size_summary, by = c("id", "date", "time_bin"))

# Summarize call counts
call_count_summary <- calls %>%
  group_by(id, calltype, date, time_bin) %>%
  summarise(
    call_count = n(),
    .groups = 'drop'
  )

# Complete the call count data to ensure all combinations of id, date, time_bin, and calltype
complete_call_count <- complete_grid %>%
  expand_grid(calltype = unique(calls$calltype)) %>%
  left_join(call_count_summary, by = c("id", "date", "time_bin", "calltype")) %>%
  mutate(call_count = replace_na(call_count, 0))


# Merge the call rate data with the subgroup size data
final_results <- results_complete %>%
  left_join(complete_call_count, by = c("id", "date", "time_bin"))


#remove rows in the dataframe where count is NA
callcount_subsize <- final_results[!is.na(final_results$calltype),]
callcount_subsize <- callcount_subsize[!is.na(callcount_subsize$mean_subgroup_size),]

callcount_subsize$call_rate <- callcount_subsize$call_count / 120 #2 minute time bins

#renaming dataframe for Odd
calling_subsize <- callcount_subsize


write.csv(calling_subsize, file = paste0(datadir, "/calling_grpsize_2mins.csv"))



ggplot(calling_subsize, aes(x = mode_subgroup_size, y = call_rate, group = mode_subgroup_size)) +
  geom_boxplot(fill = "steelblue1", color = "black") +
  labs(
    title = "Call Count vs. Mean Subgroup Size",
    y = "Call Count"
  ) +
  scale_x_discrete(name ="Mean Subgroup Size", limits = factor(1:11))+
  facet_wrap(~ calltype) +
  theme_classic()

ggsave(paste0(plot_dir, "callrate_grpsize.png"), width = 15, height = 5)


# Create a ggplot histogram
ggplot(calling_subsize, aes(x = mode_subgroup_size)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of subgroup sizes in 1Hz period",
       x = "Mode Subgroup Size",
       y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(paste0(plot_dir, "hist_subgroupsizes_1Hz.png"), width = 8, height = 5)


#adding speed of travel for each x minute bin
calling_subsize$time_bin_end <- as.POSIXct(calling_subsize$time_bin) + lubridate::seconds(120)
calling_subsize$id <- gsub("G", "", calling_subsize$id)
#add the individual index
calling_subsize$ind_idx <- match(calling_subsize$id, coati_ids$tag_id)

calling_subsize$mean_speed <- NA
calling_subsize$distance <- NA
calling_subsize$duration <- NA


i = 1

for (i in 1:nrow(calling_subsize)){
  
  #get index of the individual in the row
  ind_idx <- calling_subsize$ind_idx[i]
  #get the time index for the duration to extract the speed value
  start_idx <- which(ts == calling_subsize$time_bin[i])
  end_idx <- which(ts == calling_subsize$time_bin_end[i])
  
  # Skip the iteration if xs_sub or ys_sub are all NA's (resulting in empty vectors)
  if (length(start_idx) == 0 || length(end_idx) == 0) {
    next
  }
  
  xs_sub <- na.omit(xs[ind_idx, c(start_idx:end_idx)])
  ys_sub <- na.omit(ys[ind_idx, c(start_idx:end_idx)])
  
  #if (length(xs_sub) == 0 || length(ys_sub) == 0) {
  #  next
  #}
  
  time_interval <- 30  # Interval to filter every 10 seconds
  
  # Filter xs and ys vectors to every 10 seconds
  filtered_xs <- na.omit(xs_sub[seq(1, length(xs_sub), by = time_interval)])
  filtered_ys <- na.omit(ys_sub[seq(1, length(ys_sub), by = time_interval)])
  
  
  # Calculate the distance between consecutive points
  distances <- sqrt(diff(filtered_xs)^2 + diff(filtered_ys)^2)
  
  duration <- end_idx - start_idx
  
  # Calculate speed in meters per second (m/s)
  calling_subsize$mean_speed[i] <- sum(distances)/duration
  calling_subsize$distance[i] <- sum(distances)
  calling_subsize$duration[i] <- duration
  
}

calling_subsize_speed <- calling_subsize



#why is the old and new dataframe different??

save(calling_subsize_speed, file = paste0(datadir, "/calling_grpsize_speed_2mins_chitters.RData"))


#----------------------------------------------------------------------
#now looking at the group spread for each individual 

#get the time between the start and end of each two minute bin
calling_subsize_speed <- calling_subsize
calling_subsize_speed$time_bin_middle <- as.POSIXct(calling_subsize_speed$time_bin) + lubridate::seconds(60)
calling_subsize_speed$mean_diadic_dist_to_others <- NA
calling_subsize_speed$max_diadic_dist <- NA
calling_subsize_speed$group_ids <- NA
calling_subsize_speed$grp_size_mid <- NA
calling_subsize_speed$row_num <- 1:nrow(calling_subsize_speed)

i = 17381
for (i in 1:nrow(calling_subsize_speed)){
  
  #get index of the individual in the row
  ind_idx <- calling_subsize_speed$ind_idx[i]
  
  #get the index of the middle time 
  mid_time <- which(ts == calling_subsize_speed$time_bin_middle[i])
  
  #get all group memberships for the middle time for this row
  subgroup_comp_for_i <- subgroup_data$ind_subgroup_membership[,mid_time]
  
  #get the id of the individual on row i
  group_id_for_i <- subgroup_comp_for_i[ind_idx]
  
  #get the indexes of the individuals also in the same group as i
  group_members_for_i <- which(subgroup_comp_for_i == group_id_for_i)
  
  calling_subsize_speed$grp_size_mid[i] <- length(unique(group_members_for_i))
  calling_subsize_speed$group_ids[i] <- paste(group_members_for_i, collapse = " ")
  
  #get the xs and ys for their locations at mid_idx
  xs_all <- xs[group_members_for_i, mid_time]
  ys_all <- ys[group_members_for_i, mid_time]
  
  # Coordinates of the individual of interest
  xs_i <- xs_all[ind_idx]
  ys_i <- ys_all[ind_idx]
  
 
  # Calculate dyadic distances
  distances <- sqrt((xs_all - xs_i)^2 + (ys_all - ys_i)^2)
  
  #remove own diadic distance to get the mean diadic distance
  distances <- distances[-ind_idx]
  calling_subsize_speed$mean_diadic_dist_to_others[i] <- mean(distances)
  
  # Compute pairwise distances
  all_dists <- dist(cbind(xs_all, ys_all))
  
  # Find the maximum distance
  calling_subsize_speed$max_diadic_dist[i] <- max(all_dists)
  

}


table(is.na(calling_subsize_speed$mean_diadic_dist_to_others))

#what's the relationship between max diadic distance and mean diadic distance
ggplot(calling_subsize_speed, aes(x= max_diadic_dist, y = mean_diadic_dist_to_others))+
  geom_point()+facet_wrap(~grp_size_mid)

calling_subsize_speed_diaddist_allgrps <- calling_subsize_speed

save(calling_subsize_speed_diaddist_allgrps, file = paste0(datadir, "/calling_grpsize_speed_2mins_chitters_allgrps.RData"))


#------------------------------------------------------------------------------------------
#so now want to filter the dataframe to times when the subgroup size doesn't change in each two minute bin

calling_subsize_speed_cut <- calling_subsize_speed_diaddist_allgrps

#if there is no value for the mean diadic distance, then cut
calling_subsize_speed_cut <- calling_subsize_speed_cut[is.finite(calling_subsize_speed_cut[["mean_diadic_dist_to_others"]]), ]

calling_subsize_speed_cut$row_num <- 1:nrow(calling_subsize_speed_cut)

#this for loop is going through each two minute bin to determine whether the membership of the group has changed, and how many individuals have changed
calling_subsize_speed_cut$grp_change <- NA
i = 9401
for (i in 1:nrow(calling_subsize_speed_cut)){
  
  #get index of the individual in the row
  ind_idx <- calling_subsize_speed_cut$ind_idx[i]
  
  #get the index of the first time and last time
  start_time <- which(ts == calling_subsize_speed_cut$time_bin[i])
  end_time <- which(ts == calling_subsize_speed_cut$time_bin_end[i])
  
  if(length(end_time) == 0){
    next
  }
  
  mat <- subgroup_data$ind_subgroup_membership[,start_time:end_time]
  
  #determine whether the subgrouping changed in the two minute bin for individual in row i
   
  # Initialize a vector to store the number of changes for each timepoint
  changes <- numeric(ncol(mat) - 1)
  baseline <- mat[, 1]
    
  # Loop through each column (timepoint) starting from the second column
  for (j in 2:ncol(mat)) {
    # Compare the current column to the baseline column
    changes[j - 1] <- sum(baseline != mat[, j])
    
    # Update the baseline to the current column
    baseline <- mat[, j]
  }
  
  #if there were more than two individuals changing subgroup, then it is classed as a change
  if (max(changes, na.rm = TRUE) > 2) {
    calling_subsize_speed_cut$grp_change[i] <- "change"
  } else if (max(changes, na.rm = TRUE) > 0 & max(changes, na.rm = TRUE) <= 2) {
    calling_subsize_speed_cut$grp_change[i] <- "little change"
  } else {
    calling_subsize_speed_cut$grp_change[i] <- "no change"
  }

}

table(calling_subsize_speed_cut$grp_change)

#cut the times when there is a change in group membership
calling_subsize_speed_cut <- calling_subsize_speed_cut[!calling_subsize_speed_cut$grp_change == "change",]



save(calling_subsize_speed_cut, file = paste0(datadir, "/calling_grpsize_speed_2mins_chitters_cut_nochanges.RData"))


ggplot(calling_subsize_speed_cut, aes(x= max_diadic_dist, y = mean_diadic_dist_to_others))+
  geom_point()+facet_wrap(~grp_size_mid)











