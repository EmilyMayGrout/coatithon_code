#this script is for comparing the fission fusion events in the automated detection to the manual labels

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

group <- 'presedente' #subdirectory where the group data is stored

codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
#groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(codedir)


#read in events
events <- read.csv(paste0('Split_mechanics/',group,'_manual_split_merge_clean.csv'), sep=';')

#read in coati ids
setwd(groupdir)
load(file=paste0(group,'_coati_ids.RData'))

#read in timestamp data
load(file=paste0(group,'_xy_highres_level1.RData'))

#PROCESS
#preprocess events to...
events <- events[which(events$fission_time!='before start'),] #remove events where we missed the start
events <- events[which(events$event_type %in% c('fission','fusion')),] #only include fission and fusion events (remove 'almost fusion')

#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})

#create columns for subgroup idxs (initialize with zeros to convince R to let oyu do this)
events$group_A_idxs <- list(c(0,0,0))
events$group_B_idxs <- list(c(0,0,0))

for (i in 1:nrow(events)){
  group_A_names <- events$group_A[i]
  group_B_names <- events$group_B[i]
  
  group_A_idxs <- match_coati_names(group_A_names, coati_ids)
  group_B_idxs <- match_coati_names(group_B_names, coati_ids)
  events$group_A_idxs[i] <- list(group_A_idxs)
  events$group_B_idxs[i] <- list(group_B_idxs)
}

#merge fission_time and fusion_time columns into one
events$time_min <- paste0(events$fission_time,events$fusion_time)

#convert to POSIX
events$datetime <- as.POSIXct(paste(events$date, events$time_min), format = "%Y-%m-%d %H:%M",tz = "UTC")

#match times to get indexes into matrices
events$tidx <- match(events$datetime, ts)

#count up how many individuals are in each group
events$n_A <- unlist(lapply(events$group_A_idxs,length))
events$n_B <- unlist(lapply(events$group_B_idxs,length))

events$before_time <- events$start_time <- events$end_time <- events$after_time <- NA
events$AB_before_disp <- events$A_during_disp <- events$B_during_disp <- NA
events$split_angle <- events$turn_angle_A <- events$turn_angle_B <- NA
for(i in c(1:nrow(events))){
  print(i)
  ff_data <- analyse_ff_event(i, events, xs, ys, ts, plot=F, max_time = 600)
  if(!is.null(ff_data$disps)){
    events$AB_before_disp[i] <- ff_data$disps['AB','before']
    events$A_during_disp[i] <- ff_data$disps['A','during']
    events$B_during_disp[i] <- ff_data$disps['B','during']
  }
  events$split_angle[i] <- ff_data$split_angle
  events$turn_angle_A[i] <- ff_data$turn_angle_A
  events$turn_angle_B[i] <- ff_data$turn_angle_B
  events$before_time[i] <- ff_data$before_time
  events$start_time[i] <- ff_data$start_time
  events$end_time[i] <- ff_data$end_time
  events$after_time[i] <- ff_data$after_time
}

#now have events as the manual labels
#reduce number of columns for easier comparison
events_man <- events[, c("tidx", "datetime", "event_type", "group_A_idxs", "group_B_idxs","group_A", "group_B", "n_A", "n_B" )]
#manual labels start 2 days after gps on - but can't find where in code we excluded the first 2 days?

#make automated labels using Ari's function 
R_inner <- 15
R_outer <- 50

events_aut <- detect_fissions_and_fusions(R_inner = R_inner, R_outer = R_outer,  xs = xs, ys = ys, ts = ts, coati_ids = coati_ids, verbose = T )
events_aut_detected <- events_aut$events_detected

events_man_txt <- events[, c("tidx", "datetime", "event_type", "n_A", "n_B" )]
events_aut_txt <- events_aut_detected[, c("tidx", "datetime", "event_type", "n_A", "n_B" )]


#save the events_man and events_aut_detected to txt files 
setwd('C:/Users/egrout/Dropbox/coatithon/coatithon_code/Split_mechanics/labels/')
write.table(events_man_txt, paste0(group, "_events_man.txt"))
write.table(events_aut_txt, paste0(group, "_events_aut.txt"))


