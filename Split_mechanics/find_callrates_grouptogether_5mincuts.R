#getting speed and call rate values for 5 minute intervals 


#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'
gps_file <- "galaxy_xy_highres_level1.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#load in events
load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/galaxy_auto_ff_events_characterized.RData") #level 1

#list of Rs
R <- 50

use_machine_labels <- T

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


#finding times when the number of subgroups is 1
together_times <- as.data.frame(ts[which(subgroup_data$n_subgroups == 1)]); names(together_times) = "time"
together_times$tdiff <- c(NA,difftime(together_times$time[2:nrow(together_times)], 
                                      together_times$time[1:(nrow(together_times)-1)], units = "secs"))

# find "bouts" of group being together
together_times$bout <- numeric(nrow(together_times))
together_times$bout[1] <- 1  # Initialize the first row

n = 1
for(i in 2:nrow(together_times)){
  if(together_times$tdiff[i] > 1){ n = n + 1}
  together_times$bout[i] <- n
}

# now get the start and stop time of each together bout (the output is a list with the start and end time for each together bout)
# Extract start and stop times per bout
split_bouts <- lapply(split(together_times, together_times$bout), function(split_bout) {
  split_bout[c(1, nrow(split_bout)), "time"]
})

# Filter events with valid before_time and after_time
et <- events[!is.na(events$before_time) & !is.na(events$after_time), c("before_time", "after_time")]

# Adjust start and stop times based on event overlaps
adjusted_bouts <- lapply(split_bouts, function(together_bout) {
  # Initialize remove flag
  remove <- FALSE
  
  # Loop through events and adjust bout times if they overlap with events
  for (i in seq_len(nrow(et))) {
    bt <- et$before_time[i]
    at <- et$after_time[i]
    
    # Check if start overlaps with end of event
    if (together_bout[1] < at & together_bout[2] > at) {
      together_bout[1] <- at
    }
    
    # Check if end overlaps with before of event
    if (together_bout[2] > bt & together_bout[1] < bt) {
      together_bout[2] <- bt
    }
    
    # Check if the together bout is completely within the events
    if (together_bout[1] > bt & together_bout[2] < at) {
      remove <- TRUE
      break
    }
  }
  
  # Return the adjusted bout or NULL if it should be removed
  if (!remove) {
    return(together_bout)
  } else {
    return(NULL)
  }
})

# Remove NULL entries from the adjusted bouts list
adjusted_bouts <- Filter(Negate(is.null), adjusted_bouts)



#for looking at speed, want the bouts to be the same length, so making a different list with equal length bouts

# Function to split bouts into 5-minute intervals
split_into_intervals <- function(bout, interval_minutes = 3) {
  start_time <- bout[1]
  end_time <- bout[2]
  
  interval_seconds <- interval_minutes * 60
  intervals <- seq(from = start_time, to = end_time, by = interval_seconds)
  
  if (length(intervals) == 1) {
    # The bout is shorter than the interval duration
    return(NULL)
  }
  
  if (tail(intervals, 1) != end_time) {
    intervals <- c(intervals, end_time)
  }
  
  bout_intervals <- list()
  for (i in 1:(length(intervals) - 1)) {
    duration <- as.numeric(difftime(intervals[i + 1], intervals[i], units = "secs"))
    if (duration >= interval_seconds) {
      bout_intervals[[length(bout_intervals) + 1]] <- c(intervals[i], intervals[i + 1])
    }
  }
  
  return(bout_intervals)
}

# Create the 5-minute adjusted bouts list
adjusted_bouts_5min <- lapply(adjusted_bouts, function(bout) {
  split_into_intervals(bout, interval_minutes = 3)
})

# Flatten the list and remove NULL entries
adjusted_bouts_5min <- unlist(adjusted_bouts_5min, recursive = FALSE)
adjusted_bouts_5min <- Filter(Negate(is.null), adjusted_bouts_5min)




#get the call files 
datadir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed"

if(use_machine_labels){
  callfile <- "all_data_hms_all_ml_synched.csv"
} else {
  callfile <- 'all_data_hms_synched.csv'
  }

id_file <- 'galaxy_coati_ids.RData'

#LOAD DATA
setwd(datadir)
calls <- read.csv(callfile, header=T, sep = ',')
load(id_file)


if(use_machine_labels){
  
  bound_info <- unique(paste(calls$file, calls$filestart_UTC_soroka))
  
  # Split bound_info into file and starttime
  split_info <- strsplit(bound_info, " ")
  file <- sapply(split_info, '[', 1)
  date <- sapply(split_info, '[', 2)
  time <- sapply(split_info, '[', 3)
  
  # Combine date and time into starttime
  starttime <- paste(date, time)
  # Extract id from file
  id <- substring(file, 1, 5)
  
  # Create the new dataframe
  labeled_periods <- data.frame(file = file, id = id, starttime = starttime, stringsAsFactors = FALSE)
  
  # Initialize an empty vector to store the latest times
  # I would have got the duration of the file by loading in the wave files but I didn't have permissions to read in these files so this is the second best option - finding the latest label to infer the file duration
  latest_times <- vector("character", length = nrow(labeled_periods))
  labeled_periods$stoptime <- NA
  # Iterate over each unique file in new_data
  for (i in seq_along(labeled_periods$file)) {
    file <- labeled_periods$file[i]
    date <- as.Date(labeled_periods$starttime[i])
    # Subset calls dataframe for the current file
    calls_onefile <- calls[calls$file == file, ]
    # Find the latest time in the Start column
    latest_time <- max(calls_onefile$Start)  # Assuming Start is in POSIXct format
    #put it in the correct time (add 11:00:000)
    #latest_time <- as.character(hms(latest_time) + hms("11:00:00.000"))
    datetime <- paste(date, latest_time)
    # Store the latest time in the corresponding index of latest_times
    labeled_periods$stoptime[i] <- format(as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")+ as.difftime("11:00:00"))
    
  } 
  
} else {  
  #identify labeled periods for each individual
  startstop <- calls[which(calls$label %in% c('start','stop')),]
  startstop <- startstop[order(startstop$file, startstop$datetime_synch),]
  starts <- startstop[which(startstop$label == 'start'),]
  stops <- startstop[which(startstop$label == 'stop'),]
  labeled_periods <- data.frame(file = starts$file, id = starts$id, starttime = starts$datetime_synch, stoptime = NA)
  #for each start marker, search for the next stop marker
  for(i in 1:nrow(labeled_periods)){
    starttime <- labeled_periods$starttime[i]
    stops_for_file <- stops[which(stops$file == labeled_periods$file[i] & stops$datetime_synch > starttime),]
    next_stop <- min(stops_for_file$datetime_synch, na.rm=T)
    labeled_periods$stoptime[i] <- next_stop
  }
  
}



#remove the leading G from the tag ids for matching
labeled_periods$id <- gsub('G', '', labeled_periods$id)
#add the individual index
labeled_periods$ind_idx <- match(labeled_periods$id, coati_ids$tag_id)
rm("startstop","i","stops_for_file","next_stop","starts","stops","starttime")

#parse calls into types
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

#add coati names to column based on IDs
calls$name[calls$id == "G9463"] <- "Estrella"
calls$name[calls$id == "G9476"] <- "Luna"
calls$name[calls$id == "G9474"] <- "Saturno"
calls$name[calls$id == "G9464"] <- "Venus"
calls$name[calls$id == "G9470"] <- "Orbita"
calls$name[calls$id == "G9475"] <- "Pluto"
calls$name[calls$id == "G9480"] <- "Cometa"
calls$name[calls$id == "G9466"] <- "Lucero"
calls$name[calls$id == "G9471"] <- "Planeta"
calls$name[calls$id == "G9460"] <- "Quasar"
calls$name[calls$id == "G9467"] <- "Gus"

calls$datetime_synch_pos<-as.POSIXct(calls$datetime_synch, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 
# labeled_periods$tstart<-as.POSIXct(labeled_periods$starttime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 
# labeled_periods$tstop<-as.POSIXct(labeled_periods$stoptime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 
calls$datetime_synch_pos<-as.POSIXct(calls$datetime_synch, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
together_times$time_pos<-as.POSIXct(together_times$time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

calling_together <- calls[which(calls$datetime_synch_pos %in% together_times$time),]


# calculate call rate for each together bout
# get the duration of the bout

crt_ <- as.data.frame(coati_ids[,1]); names(crt_) <- "name"
crt_$ind_idx <- 1:nrow(crt_)
call_rates_together <- data.frame()
together_bout <- adjusted_bouts_5min[120][[1]]

for(together_bout in adjusted_bouts_5min){
  ct_ <- calling_together[which(calling_together$datetime_synch_pos >= together_bout[1] & calling_together$datetime_synch_pos <= together_bout[2]), ]
  ct_ <- ct_[!is.na(ct_$calltype),]
  
  if(nrow(ct_) > 0){
    calls_ind <- as.data.frame(table(ct_$name, ct_$calltype)); names(calls_ind)<-c("name","call_type","count")
    calls_all <- merge(calls_ind,crt_,by = "name")
    calls_all$bout_dur <- as.numeric(difftime(together_bout[2], together_bout[1], units = "secs"))
    calls_all$call_rate <- calls_all$count/calls_all$bout_dur
    calls_all$n_ind_in_bout <- length(unique(calls_all$name))
    calls_all$bout <- n
    calls_all$start_bout <- together_bout[1]
    calls_all$end_bout <- together_bout[2]
    call_rates_together <- rbind(call_rates_together, calls_all)
    
    n = n+1
  }
}

#kick out bouts that are shorter than 2min
call_rates_together_3mins <- call_rates_together[which(call_rates_together$bout_dur > 120), ]


#find total duration of times when group together to get call rates

sum(unique(call_rates_together_3mins$bout_dur))/60/60 #in hours

save(call_rates_together_3mins, file = "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/call_rates_together_3mincut_gal_ml.RData")
#load("C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/call_rates_together_gal.RData")

#just making the naming the same so I don't have to rewrite all the code
call_rates_together_5mins <- call_rates_together_3mins


#check number of events left - 65
length(unique(call_rates_together_5mins$bout))

par(mar = c(14,4,1,1))
boxplot(call_rate ~ name*call_type, data = call_rates_together_5mins[call_rates_together_5mins$n_ind_in_bout > 5,], las = 2, xlab = "", ylim = c(0, 0.6))

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
  ys_sub <- na.omit(ys[ind_idx, c(start_idx:end_idx)])
  
  # Skip the iteration if xs_sub or ys_sub are all NA's (resulting in empty vectors)
  if (length(xs_sub) == 0 || length(ys_sub) == 0) {
    next
  }
  
  
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

hist(call_rates_together_5mins$mean_speed, breaks = 100)


#no bimodal relationship in the speeds, so to decide if the group was moving or not, we could filter by distance and duration

#0.001 as cut off for not moving

#save this data frame for Odd
save(call_rates_together_5mins, file = paste0(datadir, "/calling_5mincut_baseline.RData"))

call_rates_together_5mins$starttime <- as_hms(call_rates_together_5mins$start_bout) - as_hms("5:00:00")
call_rates_together_5mins$starttime <- as_hms(call_rates_together_5mins$starttime)

#plotting distance travelled over 3 hour period
ggplot(call_rates_together_5mins[call_rates_together_5mins$call_type == "contact call",], aes(x = starttime, y = distance))+
  geom_jitter(width = 150, height = 0.001)+ #jitter by 2.5 minutes
  labs(x = "Time", y = "Distance in 3 minute bins (m)")+
  theme_classic()
ggsave(paste0(plot_dir, "distancetravelled_time_3minbins.png"), width = 10, height = 5)


