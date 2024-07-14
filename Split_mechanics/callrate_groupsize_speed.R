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
datadir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed"
callfile <- 'all_data_hms_all_ml_synched.csv' #made in coati_synch
id_file <- 'galaxy_coati_ids.RData'

#LOAD DATA
setwd(datadir)
#calls <- read.csv(callfile, header=T, sep = ',')
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
calls$calltype[calls$label == "squeal"] <- "aggression call"
calls$calltype[calls$label == "squeal chitter"] <- "aggression call"
calls$calltype[calls$label == "squeal chitter x"] <- "aggression call"
calls$calltype[calls$label == "squeal chitters"] <- "aggression call"
calls$calltype[calls$label == "low squeal"] <- "aggression call"
calls$calltype[calls$label == "chitter x"] <- "aggression call"
calls$calltype[calls$label == "squeal chittering"] <- "aggression call"

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

# Summarize call rates
call_rate_summary <- calls %>%
  group_by(id, calltype, date, time_bin) %>%
  summarise(
    call_rate = n(),
    .groups = 'drop'
  )

# Complete the call rate data to ensure all combinations of id, date, time_bin, and calltype
complete_call_rate <- complete_grid %>%
  expand_grid(calltype = unique(calls$calltype)) %>%
  left_join(call_rate_summary, by = c("id", "date", "time_bin", "calltype")) %>%
  mutate(call_rate = replace_na(call_rate, 0))


# Merge the call rate data with the subgroup size data
final_results <- results_complete %>%
  left_join(complete_call_rate, by = c("id", "date", "time_bin"))


#remove rows in the dataframe where count is NA
callrate_subsize <- final_results[!is.na(final_results$calltype),]
callrate_subsize <- callrate_subsize[!is.na(callrate_subsize$mean_subgroup_size),]

write.csv(callrate_subsize, file = paste0(datadir, "/callrate_grpsize_noNas_2mins.csv"))



ggplot(callrate_subsize, aes(x = mode_subgroup_size, y = call_rate, group = mode_subgroup_size)) +
  geom_boxplot(fill = "steelblue1", color = "black") +
  labs(
    title = "Call Rate vs. Mean Subgroup Size",
    y = "Call Rate"
  ) +
  scale_x_discrete(name ="Mean Subgroup Size", limits = factor(1:11))+
  facet_wrap(~ calltype) +
  theme_classic()

ggsave(paste0(plot_dir, "callrate_grpsize.png"), width = 15, height = 5)


# Create a ggplot histogram
ggplot(callrate_subsize, aes(x = mode_subgroup_size)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of subgroup sizes in 1Hz period",
       x = "Mean Subgroup Size",
       y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(paste0(plot_dir, "hist_subgroupsizes_1Hz.png"), width = 8, height = 5)


#adding speed of travel for each x minute bin
callrate_subsize$time_bin_end <- as.POSIXct(callrate_subsize$time_bin) + lubridate::seconds(120)
callrate_subsize$id <- gsub("G", "", callrate_subsize$id)
#add the individual index
callrate_subsize$ind_idx <- match(callrate_subsize$id, coati_ids$tag_id)

callrate_subsize$mean_speed <- NA
callrate_subsize$distance <- NA
callrate_subsize$duration <- NA


i = 1

for (i in 1:nrow(callrate_subsize)){
  
  #get index of the individual in the row
  ind_idx <- callrate_subsize$ind_idx[i]
  #get the time index for the duration to extract the speed value
  start_idx <- which(ts == callrate_subsize$time_bin[i])
  end_idx <- which(ts == callrate_subsize$time_bin_end[i])
  
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
  callrate_subsize$mean_speed[i] <- sum(distances)/duration
  callrate_subsize$distance[i] <- sum(distances)
  callrate_subsize$duration[i] <- duration
  
}

write.csv(callrate_subsize, file = paste0(datadir, "/callrate_grpsize_noNas_2mins_speed.csv"))

