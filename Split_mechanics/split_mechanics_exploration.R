#Exploring fission-fusion mechanics

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS
codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'galaxy' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

#for Emily:
codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
#groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'


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
for(i in c(1:nrow(events))[-58]){
  print(i)
  ff_data <- analyse_ff_event(i, events, xs, ys, ts, plot=F, max_time = 600)
  print('')
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


png(height = 1700, width = 2000, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfission.png'))
par(mfrow=c(2,2), mar = c(10,10,10,10),(mgp=c(3,5,1))) #bottom, left, top, right)
plot(events$A_during_disp, events$B_during_disp, pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Distance travelled during fission event (m)', cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "aquamarine3",mgp=c(5,2,.5))
plot(events$turn_angle_A, events$turn_angle_B, pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Turn angle after fission event (degrees)', cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral3",mgp=c(5,2,.5))
hist(events$AB_before_disp, main = "Distance travelled 10 minutes before displacement (m)", xlab = '', col = "aquamarine3",cex.axis = 3, cex.lab = 3, cex.main = 3,mgp=c(5,2,.5))
hist(events$split_angle, main = "Split angle between sub-groups (degrees)", xlab = '', cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral3",mgp=c(5,2,.5))

dev.off()

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

i = 2
for (i in 1:nrow(events)){
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
  events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])

  
}



