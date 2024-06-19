#Generating preprocessed data for split and merge events
#Generates a data frame giving properties of the splits and merges and their involved subgroups for manual and automated labels

#should plot the duration of events

#LIBRARY
library(lubridate)
library(scales)

#----------PARAMETERS - MODIFY HERE--------------

#which group (galaxy or presedente)
group <- 'galaxy'

#who is using (ari or emily)
user <- 'emily'

#whether to identify splits and merges automatically (if F) or use manually identified events (if T)
use_manual_events <- F

#---PARAMETERS (probably don't modify)---

#events filename - where to get the split/merge events for manually labeled events
events.filename <- paste0('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/',group,'_manual_split_merge_clean.csv') 

#radii to use
R_inner <- 15
R_outer <- 50

#------TIME ZONE------
#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS

if(user %in% c('Ari','ari')){
  codedir <- '~/Dropbox/code_ari/coatithon_code/'
  dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
  if(group == 'galaxy'){
    groupdir <- '~/Dropbox/coati/processed/galaxy/'
  } else if(group=='presedente'){
    groupdir <- '~/Dropbox/coati/processed/presedente/'
  }
} else{
  codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
  if(group == 'galaxy'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
  }
}

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(groupdir)

#read in coati ids
load(file=paste0(group,'_coati_ids.RData'))

#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})

#read in timestamp data
load(file=paste0(group,'_xy_highres_level2.RData'))

#PROCESS
setwd(codedir)
#read in events if using manual events
if(use_manual_events){
  events <- read.csv(events.filename, sep=';')
} else{ #otherwise do it automatically 
  ff_data <- detect_fissions_and_fusions(R_inner = R_inner, R_outer = R_outer, xs, ys, ts, coati_ids)
  events <- ff_data$events_detected
}



#-------CLEANING EVENTS------------------------

#preprocess events to...
if(use_manual_events){
  events <- events[which(events$fission_time!='before start'),] #remove events where we missed the start
  events <- events[which(events$event_type %in% c('fission','fusion')),] #only include fission and fusion events (remove 'almost fusion')
}

#using google sheets (Finding correct ff from level 1) to remove the events which are not true events (either caused by one individuals moving groups and temporarily joining the groups or due to GPS error when collars started recording)
#also removing events which are due to males ff-ing
#these indexes are for the time index so should work for events in level 1 and 2 data
if(group == 'galaxy'){
events <- events[!events$tidx %in% c(46805,122721,130464,134273,134288,141605,143333),]

} else{(group=='presedente')
  
  #removing false events
  events <- events[!events$tidx %in% c(16099,16123,106812,106820,112549,120551,120691,122811,126046,130953,131232,131657,131736,131745,131746,131795,133444,133926,134385,135523,136192,136279,137678,138010,140811,141521,149395,149442,153490,153653,153948,154409,154493,154661),]
  #these are the event idxs which match the tidx: c(32,33,199,200,211,224,225,235,255,267,269,272:275,277,281,286,289,296,300,301,305,306,315,318,347,348,366,369,370,372,373,374)
  
  #removing the MALE events - got numbers from google sheets (manually found these as male events by going through each event with analyse_ff_events function in the identify_splits_and_merges code)
 #events <- events[!events$tidx %in% c(330,822,2159,2438,2505,6364,6992,9446,9566,9877,10010,10131,10313,12985,16747,17107,19229,19357,20503,21200,23996,24443,26658,28965,30814,31967,33159,33592,42029,44951,45094,45560,49299,49469,54764,55607,55819,57222,57275,57545,59151,64502,64521,64720,64841,65443,65954,66028,74171,74173,74466,76501,76815,76990,77475,77523,77622,81869,83212,83438,86147,86192,86334,86485,88445,88625,91699,92030,92400,92551 ,92804,93041,96484,100519,105525,105798,105823,105971,106059,107389,112609,113319,120961,121008,121359,122213,122773,123393,123704,123854,125431,125626,126237,126514,126592,126693,126784,126849,127456,128094,128694,129221,131458,132232,133963,134222,134626,135406,136771,137198,138010,138477,139508,143063,143532,143846,144407,146494,146900,147804,148288,150425,150527,151380,151833,152384,152766,152780,153279,153411,157533,157535,157976,159280,159767),]
  
   #these are the event ids which match these tidx: c(4,5,6,7,16:23, 27,35,36,44,45,49,50,55,57,60,66:70,77,79,82,84,88,89,102:104,109:111,114,117:123,134:137,139,140,142:144,156:158,160:165,168:174,178,193:197,202,212,214,227:229,231,234,239,241,242,251,252,257:266,270,278,287, 288,291,295,303,304,306,307,310,324,325,327,329,333,334,336,340,350,351,354,356,358:360,363,365,379,380,382,385,387),]
  
  #removing events which are too early to be real
  events <- events[!events$tidx %in% c(14,117,10823,10939,21672,53990,54019,64841,86485,97329,108015,118821,140414),]
  
  }

#-----------------------------------------------

#CHARACTERIZE EVENTS

#create columns for subgroup idxs (initialize with zeros to convince R to let you do this)
events$group_A_idxs <- list(c(0,0,0))
events$group_B_idxs <- list(c(0,0,0))


for (i in 1:nrow(events)){
  
  group_A_names <- events$group_A[i][[1]]
  group_B_names <- events$group_B[i][[1]]
  
  #if auto identified events, put in the right form (string) for matching
  if(!use_manual_events){
    group_A_names <- paste(group_A_names, collapse= ', ')
    group_B_names <- paste(group_B_names, collapse= ', ')
  }
  
  group_A_idxs <- match_coati_names(group_A_names, coati_ids)
  group_B_idxs <- match_coati_names(group_B_names, coati_ids)
  events$group_A_idxs[i] <- list(group_A_idxs)
  events$group_B_idxs[i] <- list(group_B_idxs)
}

#merge fission_time and fusion_time columns into one
if(use_manual_events){
  events$time_min <- paste0(events$fission_time,events$fusion_time)
  #convert to POSIX
  events$datetime <- as.POSIXct(paste(events$date, events$time_min), format = "%Y-%m-%d %H:%M",tz = "UTC")
}

#match times to get indexes into matrices
events$tidx <- match(events$datetime, ts)

#count up how many individuals are in each group
events$n_A <- unlist(lapply(events$group_A_idxs,length))
events$n_B <- unlist(lapply(events$group_B_idxs,length))

events$before_time <- events$start_time <- events$end_time <- events$after_time <- NA
events$AB_before_disp <- events$AB_after_disp <- events$A_during_disp <- events$B_during_disp <- NA
events$split_angle <- events$turn_angle_A <- events$turn_angle_B <- NA

i = 3
for(i in c(1:nrow(events))){
  print(i)
  ff_data <- analyse_ff_event(i, events, xs, ys, ts, plot=F, max_time = 700) #with 700s, catches more start time of events which are still accurate to the event (not picking up a time from a different event)
  if(!is.null(ff_data$disps)){
    events$AB_before_disp[i] <- ff_data$disps['AB','before']
    events$AB_after_disp[i] <- ff_data$disps['AB','after']
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

sum(is.na(events$after_time))


#adding age class as a number in coati_ids
coati_ids$age_class <- NA
coati_ids$age_class[coati_ids$age == "Juvenile"] <- 1
coati_ids$age_class[coati_ids$age == "Sub-adult"] <- 2
coati_ids$age_class[coati_ids$age == "Adult"] <- 3


#this for loop is running through each groups IDs and assigning their age class and getting the groups average age
# events$A_average_grp_age <- NA
# events$A_age_each_ind <- events$group_A_idxs
# events$B_average_grp_age <- NA
# events$B_age_each_ind <- events$group_B_idxs

# for (i in 1:nrow(events)){
#   
#   if (group == "galaxy"){
#     v1 <- unlist(events$A_age_each_ind[i])
#     v1[v1 == 1] <- 1
#     v1[v1 == 2] <- 3
#     v1[v1 == 3] <- 3
#     v1[v1 == 4] <- 3
#     v1[v1 == 5] <- 3
#     v1[v1 == 6] <- 3
#     v1[v1 == 7] <- 2
#     v1[v1 == 8] <- 2
#     v1[v1 == 9] <- 2
#     v1[v1 == 10] <- 3
#     v1[v1 == 11] <- 3
#     events$A_average_grp_age[i] <- mean(v1)
#     events$A_n_adults[i] <- length(v1[v1=="3"])
#     events$A_n_subadults[i] <- length(v1[v1=="2"])
#     events$A_n_juveniles[i] <- length(v1[v1=="1"])
#     events$A_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
#     events$A_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
#     events$A_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
#     events$A_subgroup_size[i] <- length(v1)
#     events$A_age_each_ind[i] <- relist((v1), skeleton=events$A_age_each_ind[i])
#     
#     #for  group
#     v1 <- unlist(events$B_age_each_ind[i])
#     v1[v1 == 1] <- 1
#     v1[v1 == 2] <- 3
#     v1[v1 == 3] <- 3
#     v1[v1 == 4] <- 3
#     v1[v1 == 5] <- 3
#     v1[v1 == 6] <- 3
#     v1[v1 == 7] <- 2
#     v1[v1 == 8] <- 2
#     v1[v1 == 9] <- 2
#     v1[v1 == 10] <- 3
#     v1[v1 == 11] <- 3
#     events$B_average_grp_age[i] <- mean(v1)
#     events$B_n_adults[i] <- length(v1[v1=="3"])
#     events$B_n_subadults[i] <- length(v1[v1=="2"])
#     events$B_n_juveniles[i] <- length(v1[v1=="1"])
#     events$B_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
#     events$B_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
#     events$B_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
#     events$B_subgroup_size[i] <- length(v1)
#     events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])
#     
#   } else if (group == "presedente"){
#     
#     v1 <- unlist(events$A_age_each_ind[i])
#     v1[v1 == 1] <- 3
#     v1[v1 == 2] <- 1
#     v1[v1 == 3] <- 1
#     v1[v1 == 4] <- 1
#     v1[v1 == 5] <- 3
#     v1[v1 == 6] <- 2
#     v1[v1 == 7] <- 2
#     v1[v1 == 8] <- 3
#     v1[v1 == 9] <- 3
#     v1[v1 == 10] <- 2
#     v1[v1 == 11] <- 1
#     v1[v1 == 12] <- 1
#     v1[v1 == 13] <- 3
#     v1[v1 == 14] <- 2
#     v1[v1 == 15] <- 1
#     v1[v1 == 16] <- 1
#     v1[v1 == 17] <- 2
#     v1[v1 == 18] <- 3
#     v1[v1 == 19] <- 2
#     v1[v1 == 20] <- 3
#     v1[v1 == 21] <- 3
#     v1[v1 == 22] <- 1
#     
#     events$A_average_grp_age[i] <- mean(v1)
#     events$A_n_adults[i] <- length(v1[v1=="3"])
#     events$A_n_subadults[i] <- length(v1[v1=="2"])
#     events$A_n_juveniles[i] <- length(v1[v1=="1"])
#     events$A_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
#     events$A_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
#     events$A_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
#     
#     events$A_subgroup_size[i] <- length(v1)
#     events$A_age_each_ind[i] <- relist((v1), skeleton=events$A_age_each_ind[i])
#     
#     #for  group
#     v1 <- unlist(events$B_age_each_ind[i])
#     v1[v1 == 1] <- 3
#     v1[v1 == 2] <- 1
#     v1[v1 == 3] <- 1
#     v1[v1 == 4] <- 1
#     v1[v1 == 5] <- 3
#     v1[v1 == 6] <- 2
#     v1[v1 == 7] <- 2
#     v1[v1 == 8] <- 3
#     v1[v1 == 9] <- 3
#     v1[v1 == 10] <- 2
#     v1[v1 == 11] <- 1
#     v1[v1 == 12] <- 1
#     v1[v1 == 13] <- 3
#     v1[v1 == 14] <- 2
#     v1[v1 == 15] <- 1
#     v1[v1 == 16] <- 1
#     v1[v1 == 17] <- 2
#     v1[v1 == 18] <- 3
#     v1[v1 == 19] <- 2
#     v1[v1 == 20] <- 3
#     v1[v1 == 21] <- 3
#     v1[v1 == 22] <- 1
#     events$B_average_grp_age[i] <- mean(v1)
#     events$B_n_adults[i] <- length(v1[v1=="3"])
#     events$B_n_subadults[i] <- length(v1[v1=="2"])
#     events$B_n_juveniles[i] <- length(v1[v1=="1"])
#     events$B_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
#     events$B_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
#     events$B_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
#     
#     events$B_subgroup_size[i] <- length(v1)
#     events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])
#   }else {print("fail")}
# }

if(use_manual_events){
  save(list = c('events'), file = paste0(groupdir, group,'_manual_ff_events_characterized.RData'))
} else{
  save(list = c('events'), file = paste0(groupdir, group,'_auto_ff_events_characterized.RData'))
}



#Check whether the before period overlaps with a preceding event involving some of the same individuals
events$ovlp_before <- events$ovlp_during <- events$ovlp_after <- list(c(0))
events$n_ovlp_before <- events$n_ovlp_during <- events$n_ovlp_after <- NA
for(i in 1:nrow(events)){
  
  #before time and after time of the current event
  before_time_curr <- events$before_time[i]
  start_time_curr <- events$start_time[i]
  end_time_curr <- events$end_time[i]
  after_time_curr <- events$after_time[i]
  
  #individuals
  inds_in_event <- c(events$group_A_idxs[i][[1]], events$group_B_idxs[i][[1]])
  
  #lists of overlapping events - set to empty lists initially
  ovlp_before <- ovlp_during <- ovlp_after <- list()
  #loop through other events to see if they overlap
  for(j in 1:nrow(events)){
    
    #skip the current event
    if(i==j){
      next
    }
    
    #if the events don't involve any of the same individuals, skip
    inds_in_other_event <- c(events$group_A_idxs[j][[1]], events$group_B_idxs[j][[1]])
    if(length(intersect(inds_in_event, inds_in_other_event))==0){
      next
    }
    
    #get start and end time of the other event
    start_time_other <- events$start_time[j]
    end_time_other <- events$end_time[j]
    
    #check if another event (start - end period) falls within before_time and start_time of the current event
    if(!is.na(before_time_curr) & !is.na(start_time_curr) & !is.na(end_time_other) & !is.na(start_time_other)){
      if(before_time_curr <= end_time_other & start_time_curr >= start_time_other){
        ovlp_before <- c(ovlp_before, j)
      }
    }
    #check if another event (start - end period) falls within start_time and end_time of the current event
    if(!is.na(start_time_curr) & !is.na(end_time_curr) & !is.na(end_time_other) & !is.na(start_time_other)){
      if(start_time_curr <= end_time_other & end_time_curr >= start_time_other){
        ovlp_during <- c(ovlp_during, j)
      }
    }
    #check if another event (start - end period) falls within end_time and after_time of the current event
    if(!is.na(end_time_curr) & !is.na(after_time_curr) & !is.na(end_time_other) & !is.na(start_time_other)){
      if(end_time_curr <= end_time_other & after_time_curr >= start_time_other){
        ovlp_after <- c(ovlp_after, j)
      }
    }
  }
  events$ovlp_before[i] <- list(ovlp_before)
  events$ovlp_during[i] <- list(ovlp_during)
  events$ovlp_after[i] <- list(ovlp_after)
  
}

events$n_ovlp_before <- sapply(events$ovlp_before, length)
events$n_ovlp_during <- sapply(events$ovlp_during, length)
events$n_ovlp_after <- sapply(events$ovlp_after, length)

#look at the speed distributions to decide where to do the cut-off for stationary/slow and moving
dev.off()
hist(rbind(events$A_during_disp, events$B_during_disp), breaks = 40)


#visualising the events with Eli's code

#new function currently get's error after "Identifying changes in group membership"
ff_output <- identify_splits_and_merges(R_inner = R_inner, R_outer = R_outer, xs, ys, ts, breaks = c(1, length(ts)+1), names = coati_ids)
                                        

events <- ff_data$events_detected



#-------------------------------------------------------------------------

# plot_animated_events(1, ff_output, xs, ys, ts, surrounding_mins = 5, 60, color_only_main_groups = F, coati_ids$name, plot_dir, plot_interpolated = F, interpolated_xs = NULL, interpolated_ys = NULL, interpolated_timestamps = NULL,interpolated_fixes_per_minute = NULL)
# 
# 
# fixes_per_minute <- 60
# surrounding_mins = 5
# plot_interpolated <- F
# interpolated_xs <- NULL
# interpolated_ys <- NULL
# interpolated_timestamps <- NULL
# interpolated_fixes_per_minute <- NULL
# timestamps <- ts
# r <- 1
# ids <- coati_ids$name
# 
# plot_animated_events <- function(r, ff_output, xs, ys, timestamps, surrounding_mins = 5, fixes_per_minute, color_only_main_groups = F, ids, plot_dir,
#                                  plot_interpolated = F, interpolated_xs = NULL, interpolated_ys = NULL, interpolated_timestamps = NULL,
#                                  interpolated_fixes_per_minute = NULL){
#   
#   events <- ff_output$events_detected
#   
#   if(plot_interpolated){
#     if(is.null(interpolated_xs) | is.null(interpolated_ys | is.null(interpolated_timestamps) | is.null(interpolated_fixes_per_minute)))
#       stop('if plotting interpolated data, please provide interpolated xs, ys, and timestamps')
#     
#     ## Make the event timestamps and the plotting timestamps and the links between them
#     ## Create mapping between timestamps and interpolated_timestamps
#     timestamps_mapping <- as.vector(tidyr::fill(data.frame(x = match(interpolated_timestamps, timestamps)), x)$x)
#     
#     time_of_event_interpolated <- which(interpolated_timestamps == timestamps[events$tidx[r]])
#     time_of_event <- which(timestamps == timestamps[events$tidx[r]])
#     time_idxs_plotting <- (time_of_event_interpolated - (surrounding_mins*interpolated_fixes_per_minute)) : (time_of_event_interpolated + (surrounding_mins*interpolated_fixes_per_minute))
#     time_idxs_event <- (time_of_event - (surrounding_mins*fixes_per_minute)) : (time_of_event + (surrounding_mins*fixes_per_minute))
#   }else{
#     time_of_event <- which(timestamps == timestamps[events$tidx[r]])
#     time_idxs_plotting <- time_idxs_event <- (time_of_event - (surrounding_mins*fixes_per_minute)) : (time_of_event + (surrounding_mins*fixes_per_minute))
#   }
#   
#   ids_in_event <- ids[events$big_group_idxs[r][[1]]]
#   ids_in_A <-ids[events$group_A_idxs[r][[1]]]
#   ids_in_B <-ids[events$group_B_idxs[r][[1]]]
#   ids_in_C <-ids[events$group_C_idxs[r][[1]]]
#   ids_in_D <-ids[events$group_D_idxs[r][[1]]]
#   
#   ## Data frame for drawing lines between individuals who are within inner threshold
#   togethers <- ff_output$together[,,time_idxs_event]
#   together_df_list <- list()
#   for(tid in time_idxs_event){
#     together_idxs <- which(ff_output$together[,,tid], arr.ind = T)
#     
#     ## get the 1hz time indexes that are greater than tid and less then the next time
#     if(plot_interpolated){
#       matching_1hz_tids <- which(interpolated_timestamps >= timestamps[tid] & interpolated_timestamps < timestamps[tid+1])
#       together_df_list[[length(together_df_list)+1]] <- data.frame(x = as.vector(interpolated_xs[together_idxs[,1], matching_1hz_tids]),xend = as.vector(interpolated_xs[together_idxs[,2], matching_1hz_tids]), y = as.vector(interpolated_ys[together_idxs[,1], matching_1hz_tids]),yend = as.vector(interpolated_ys[together_idxs[,2], matching_1hz_tids]),time = rep(interpolated_timestamps[matching_1hz_tids], each = nrow(together_idxs)))
#     }else{
#       together_df_list[[length(together_df_list)+1]] <- data.frame(x = as.vector(xs[together_idxs[,1], timestamps[tid]]), xend = as.vector(xs[together_idxs[,2], timestamps[tid]]), y = as.vector(ys[together_idxs[,1], timestamps[tid]]), yend = as.vector(ys[together_idxs[,2], timestamps[tid]]), time = rep(timestamps[tid], nrow(together_idxs)))
#     }
#   }
#   together_df <- do.call(rbind, together_df_list)
#   
#   ## Combine into data frame for plotting
#   if(plot_interpolated){
#     plot_df <- ariformat_to_df(time_indices_selected = time_idxs_plotting, xs = interpolated_xs, ys = interpolated_ys, ids = ids, timestamps = interpolated_timestamps)
#   }else{
#     plot_df <- ariformat_to_df(time_indices_selected = time_idxs_plotting, xs = xs, ys = ys, ids = ids, timestamps = timestamps)
#   }
#   ## Assign individuals to subgroups for coloring points
#   plot_df$group <- 'uninvolved'
#   if(events$event_type[r] == 'fission'){
#     ## Before event time, all involved individuals are in the big group
#     plot_df[plot_df$id %in% ids_in_event & plot_df$time < timestamps[events$tidx[r]],'group'] <- 'big_group'
#     
#     ## After (or equal to) assign to the subgroups
#     for(letter in c('A', 'B', 'C', 'D')){
#       subgroup_ids <- get(paste0('ids_in_', letter))
#       
#       ## If subgroup is empty, skip (common for C and D, should never happen for A and B)
#       if(all(is.na(subgroup_ids))){
#         next
#       }
#       plot_df[plot_df$id %in% subgroup_ids & plot_df$time >= timestamps[events$tidx[r]],'group'] <- letter
#     }
#   }else if(events$event_type[r] == 'fusion'){
#     ## Before event time, all involved individuals are in the big group
#     plot_df[plot_df$id %in% ids_in_event & plot_df$time >= timestamps[events$tidx[r]],'group'] <- 'big_group'
#     
#     ## After (or equal to) assign to the subgroups
#     for(letter in c('A', 'B', 'C', 'D')){
#       subgroup_ids <- get(paste0('ids_in_', letter))
#       
#       ## If subgroup is empty, skip (common for C and D, should never happen for A and B)
#       if(all(is.na(subgroup_ids))){
#         next
#       }
#       plot_df[plot_df$id %in% subgroup_ids & plot_df$time < timestamps[events$tidx[r]],'group'] <- letter
#     }
#   }
#   
#   if(color_only_main_groups){
#     color_lookup <- data.frame(col =  c('firebrick2', 'blue1', 'purple4', 'gray', 'orange', 'dodgerblue'),
#                                group = c('A', 'B', 'big_group', 'uninvolved', 'C', 'D'))
#     plot_df$col <- left_join(plot_df, color_lookup)$col
#   }else{
#     color_lookup <- data.frame(id = ids,
#                                r = runif(length(ids), 0,1),
#                                g = runif(length(ids), 0,1),
#                                b = runif(length(ids), 0,1))
#     
#     if(plot_interpolated){
#       ### For every 1 second, get matching 30s timestamp index using the 1hz to 30s mapping
#       ts <- timestamps_mapping[match(plot_df$time, interpolated_timestamps)]
#     }else{
#       ts <- match(plot_df$time, timestamps)
#     }
#     
#     
#     plot_df$col <- NA
#     for(i in 1:nrow(plot_df)){
#       
#       ## Skip for individuals who have no gps data (e.g., LEXI)
#       if(is.na(plot_df[i,'x']))
#         next
#       
#       id <- match(plot_df$id[i], ids)
#       ## Get group number for the row
#       group_num <- ffs$groups[id,ts[i]]
#       
#       ## Find matching group in groups_list
#       group <- ffs$groups_list[[ts[i]]][[group_num]]
#       
#       ## Assign color based on mean of individual values
#       rgbs <- apply(X = color_lookup[group,c('r', 'g', 'b')], FUN = mean, MARGIN = 2)
#       plot_df$col[i] <- rgb(rgbs[1], rgbs[2], rgbs[3])
#       
#     }
#   }
#   
#   ## Center the plot and standardize axes to meters from center
#   centroid_group <- c(mean(plot_df[plot_df$id %in% ids_in_event,'x'], na.rm = T),
#                       mean(plot_df[plot_df$id %in% ids_in_event,'y'], na.rm = T))
#   plot_df$x <- plot_df$x - centroid_group[1]
#   plot_df$y <- plot_df$y - centroid_group[2]
#   
#   together_df$x <- together_df$x - centroid_group[1]
#   together_df$xend <- together_df$xend - centroid_group[1]
#   together_df$y <- together_df$y - centroid_group[2]
#   together_df$yend <- together_df$yend - centroid_group[2]
#   together_df <- together_df[together_df$time <= max(plot_df$time),]
#   
#   p <- ggplot(data = plot_df, aes(x = x, y = y, group = id, color = col))+
#     geom_segment(data = together_df, aes(xend = xend, yend = yend, x = x, y = y), color = 'gray', alpha = 0.1, inherit.aes = F)+
#     geom_point(alpha = 0.3) + 
#     xlim(c(-500, 500))+
#     ylim(c(-500, 500))+
#     scale_color_identity()+
#     theme_classic()+
#     coord_equal()+
#     theme(legend.position = 'none')+
#     transition_time(time)
#   #shadow_wake(wake = 0.3, exclude_layer = 1)
#   
#   
#   
#   pa <- animate(p, height = 4, width = 4, units = 'in', res = 300, type = 'cairo', renderer = gifski_renderer())
#   anim_save(paste0('animated_event_', r,'.gif'), pa, path = plot_dir)
# }
# 
# 
