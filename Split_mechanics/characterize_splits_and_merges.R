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
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
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
  events <- events[!events$tidx %in% c(330,822,2159,2438,2505,6364,6992,9446,9566,9877,10010,10131,10313,12985,16747,17107,19229,19357,20503,21200,23996,24443,26658,28965,30814,31967,33159,33592,42029,44951,45094,45560,49299,49469,54764,55607,55819,57222,57275,57545,59151,64502,64521,64720,64841,65443,65954,66028,74171,74173,74466,76501,76815,76990,77475,77523,77622,81869,83212,83438,86147,86192,86334,86485,88445,88625,91699,92030,92400,92551 ,92804,93041,96484,100519,105525,105798,105823,105971,106059,107389,112609,113319,120961,121008,121359,122213,122773,123393,123704,123854,125431,125626,126237,126514,126592,126693,126784,126849,127456,128094,128694,129221,131458,132232,133963,134222,134626,135406,136771,137198,138010,138477,139508,143063,143532,143846,144407,146494,146900,147804,148288,150425,150527,151380,151833,152384,152766,152780,153279,153411,157533,157535,157976,159280,159767),]
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
  ff_data <- analyse_ff_event(i, events, xs, ys, ts, plot=T, max_time = 700) #with 700s, catches more start time of events which are still accurate to the event (not picking up a time from a different event)
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
events$A_average_grp_age <- NA
events$A_age_each_ind <- events$group_A_idxs
events$B_average_grp_age <- NA
events$B_age_each_ind <- events$group_B_idxs

for (i in 1:nrow(events)){
  
  if (group == "galaxy"){
    v1 <- unlist(events$A_age_each_ind[i])
    v1[v1 == 1] <- 1
    v1[v1 == 2] <- 3
    v1[v1 == 3] <- 3
    v1[v1 == 4] <- 3
    v1[v1 == 5] <- 3
    v1[v1 == 6] <- 3
    v1[v1 == 7] <- 2
    v1[v1 == 8] <- 2
    v1[v1 == 9] <- 2
    v1[v1 == 10] <- 3
    v1[v1 == 11] <- 3
    events$A_average_grp_age[i] <- mean(v1)
    events$A_n_adults[i] <- length(v1[v1=="3"])
    events$A_n_subadults[i] <- length(v1[v1=="2"])
    events$A_n_juveniles[i] <- length(v1[v1=="1"])
    events$A_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
    events$A_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
    events$A_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
    events$A_subgroup_size[i] <- length(v1)
    events$A_age_each_ind[i] <- relist((v1), skeleton=events$A_age_each_ind[i])
    
    #for  group
    v1 <- unlist(events$B_age_each_ind[i])
    v1[v1 == 1] <- 1
    v1[v1 == 2] <- 3
    v1[v1 == 3] <- 3
    v1[v1 == 4] <- 3
    v1[v1 == 5] <- 3
    v1[v1 == 6] <- 3
    v1[v1 == 7] <- 2
    v1[v1 == 8] <- 2
    v1[v1 == 9] <- 2
    v1[v1 == 10] <- 3
    v1[v1 == 11] <- 3
    events$B_average_grp_age[i] <- mean(v1)
    events$B_n_adults[i] <- length(v1[v1=="3"])
    events$B_n_subadults[i] <- length(v1[v1=="2"])
    events$B_n_juveniles[i] <- length(v1[v1=="1"])
    events$B_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
    events$B_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
    events$B_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
    events$B_subgroup_size[i] <- length(v1)
    events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])
    
  } else if (group == "presedente"){
    
    v1 <- unlist(events$A_age_each_ind[i])
    v1[v1 == 1] <- 3
    v1[v1 == 2] <- 1
    v1[v1 == 3] <- 1
    v1[v1 == 4] <- 1
    v1[v1 == 5] <- 3
    v1[v1 == 6] <- 2
    v1[v1 == 7] <- 2
    v1[v1 == 8] <- 3
    v1[v1 == 9] <- 3
    v1[v1 == 10] <- 2
    v1[v1 == 11] <- 1
    v1[v1 == 12] <- 1
    v1[v1 == 13] <- 3
    v1[v1 == 14] <- 2
    v1[v1 == 15] <- 1
    v1[v1 == 16] <- 1
    v1[v1 == 17] <- 2
    v1[v1 == 18] <- 3
    v1[v1 == 19] <- 2
    v1[v1 == 20] <- 3
    v1[v1 == 21] <- 3
    v1[v1 == 22] <- 1
    
    events$A_average_grp_age[i] <- mean(v1)
    events$A_n_adults[i] <- length(v1[v1=="3"])
    events$A_n_subadults[i] <- length(v1[v1=="2"])
    events$A_n_juveniles[i] <- length(v1[v1=="1"])
    events$A_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
    events$A_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
    events$A_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
    
    events$A_subgroup_size[i] <- length(v1)
    events$A_age_each_ind[i] <- relist((v1), skeleton=events$A_age_each_ind[i])
    
    #for  group
    v1 <- unlist(events$B_age_each_ind[i])
    v1[v1 == 1] <- 3
    v1[v1 == 2] <- 1
    v1[v1 == 3] <- 1
    v1[v1 == 4] <- 1
    v1[v1 == 5] <- 3
    v1[v1 == 6] <- 2
    v1[v1 == 7] <- 2
    v1[v1 == 8] <- 3
    v1[v1 == 9] <- 3
    v1[v1 == 10] <- 2
    v1[v1 == 11] <- 1
    v1[v1 == 12] <- 1
    v1[v1 == 13] <- 3
    v1[v1 == 14] <- 2
    v1[v1 == 15] <- 1
    v1[v1 == 16] <- 1
    v1[v1 == 17] <- 2
    v1[v1 == 18] <- 3
    v1[v1 == 19] <- 2
    v1[v1 == 20] <- 3
    v1[v1 == 21] <- 3
    v1[v1 == 22] <- 1
    events$B_average_grp_age[i] <- mean(v1)
    events$B_n_adults[i] <- length(v1[v1=="3"])
    events$B_n_subadults[i] <- length(v1[v1=="2"])
    events$B_n_juveniles[i] <- length(v1[v1=="1"])
    events$B_proportion_adults[i] <- length(v1[v1=="3"])/length(v1)
    events$B_proportion_subadults[i] <- length(v1[v1=="2"])/length(v1)
    events$B_proportion_juveniles[i] <- length(v1[v1=="1"])/length(v1)
    
    events$B_subgroup_size[i] <- length(v1)
    events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])
  }else {print("fail")}
}

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




