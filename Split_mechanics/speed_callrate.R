# This script is looking at travel speed and call rate when the group is travelling

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'
gps_file <- "galaxy_xy_highres_level2.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/call_rates_together_5mincut_gal.RData")

calls <- read.csv("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/all_data_hms_synched.csv", header=T, sep = ',')


#finding the travel speed of each individual for each bout duration

call_rates_together_5mins$mean_speed <- NA
call_rates_together_5mins$max_speed <- NA
call_rates_together_5mins$distance <- NA
call_rates_together_5mins$duration <- NA

#i = 1

for (i in 1:nrow(call_rates_together_5mins)){
  
  #get index of the individual in the row
  ind_idx <- call_rates_together_5mins$ind_idx[i]
  #get the time index for the duration to extract the speed value
  start_idx <- which(ts == call_rates_together_5mins$start_bout[i])
  end_idx <- which(ts == call_rates_together_5mins$end_bout[i])
  
  xs_sub <- na.omit(xs[ind_idx, c(start_idx:end_idx)])
  ys_sub <-na.omit(ys[ind_idx, c(start_idx:end_idx)])
  
  time_interval <- 30  # Interval to filter every 10 seconds
  
  # Filter xs and ys vectors to every 10 seconds
  filtered_xs <- na.omit(xs_sub[seq(1, length(xs_sub), by = time_interval)])
  filtered_ys <- na.omit(ys_sub[seq(1, length(ys_sub), by = time_interval)])
  
  
  # Calculate the distance between consecutive points
  distances <- sqrt(diff(filtered_xs)^2 + diff(filtered_ys)^2)
  
  duration <- end_idx - start_idx
  
  # Calculate speed in meters per second (m/s)
  call_rates_together_5mins$mean_speed[i] <- mean(distances/duration)
  call_rates_together_5mins$max_speed[i]<- max(distances/duration)
  call_rates_together_5mins$distance[i] <- sum(distances)
  call_rates_together_5mins$duration[i] <- duration
  
}


save(call_rates_together_5mins, file = paste0(data_dir, "/calling_5mincut_baseline.RData"))

#just interested in the contact calls
contactcall_speed <- call_rates_together_5mins[call_rates_together_5mins$call_type == "contact call",]

custom_colors <- c("#A6CEE3", "#1F78B4", "#A2CD5A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CDB5CD", "#8B4789", "#00868B")

ggplot(contactcall_speed, aes(mean_speed, call_rate))+
  geom_point(aes(colour = factor(name))) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") + # Line of best fit
  labs(x ="Mean Speed (m/s)", y = "Contact call rate (per minute)", colour = "ID")+
  scale_colour_manual(values = custom_colors) +
  theme_classic()

ggplot(contactcall_speed, aes(mean_speed))+
  geom_histogram() +
  theme_classic()



