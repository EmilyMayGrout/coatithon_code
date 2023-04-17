#Exploring fission-fusion mechanics

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS
codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'presedente' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

#for Emily:
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
  events$B_subgroup_size[i] <- length(v1)
  events$B_age_each_ind[i] <- relist((v1), skeleton=events$B_age_each_ind[i])
  
  
}


### FISSION PLOTTING ###

#combine the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp to make abline
combined_1 <- data.frame(events$A_subgroup_size[fis])
colnames(combined_1) <- "sub_size"
disp_1 <- data.frame(events$A_during_disp[fis])
colnames(disp_1) <- "disp"
combined_2 <- data.frame(events$B_subgroup_size[fis])
colnames(combined_2) <- "sub_size"
disp_2 <- data.frame(events$B_during_disp[fis])
colnames(disp_2) <- "disp"

combined <- rbind(combined_1,combined_2)
disp <- rbind(disp_1, disp_2)
comb <- cbind(combined, disp)

fis <- events$event_type == "fission"

png(height = 1800, width = 2000, units = 'px', filename = paste0(plot_dir,'subgroup_dist_travelled_duringfission.png'))
par(mfrow=c(2,2), mar = c(10,10,10,10),(mgp=c(3,5,1))) #bottom, left, top, right)
#plot 1:
plot(events$A_during_disp[fis], events$B_during_disp[fis], pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Distance traveled during fission event (m)', cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "aquamarine3",mgp=c(5,2,.5))
abline(lm(events$B_during_disp[fis] ~ events$A_during_disp[fis]))

#plot 2:
#combined the A_subgroup_size and B_sub_group_size with A_during_disp and B_during_disp so I could add the abline
plot(comb$disp, comb$sub_size, pch = 20, ylim = c(0, 10), xlab = "Distance traveled during fission event (m)", ylab = "Sub-group size", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2",mgp=c(5,2,.5))
abline(lm(comb$sub_size ~ comb$disp))

hist(events$AB_before_disp[fis], main = "Distance traveled 10 minutes before displacement (m)", xlab = '', col = "aquamarine3",cex.axis = 3, cex.lab = 3, cex.main = 3,mgp=c(5,2,.5))

hist(events$split_angle[fis], main = "Split angle between sub-groups after fission (degrees)", xlab = '', cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral3",mgp=c(5,2,.5))

dev.off()



#removed plot as the turn_angle gives less info than the split angle (because not positive or negative)
plot(events$turn_angle_A[fis], events$turn_angle_B[fis], pch = 20, xlab = "sub-group A", ylab = "sub-group B", main = 'Turn angle after fission event (degrees)',col = "coral3",mgp=c(5,2,.5))

#same as plot 2 from above for loop, but the data here is directly from the events dataframe (rather than my rbind data frame to make the abline)
plot(events$A_during_disp[fis],events$A_subgroup_size[fis], pch = 20, ylim = c(0, 10), xlab = "Distance traveled during fission event (m)", ylab = "Sub-group size", cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2",mgp=c(5,2,.5))
points(events$B_during_disp[fis], events$B_subgroup_size[fis], pch = 20, cex = 4, cex.axis = 3, cex.lab = 3, cex.main = 3,col = "coral2")
abline(lm(events$B_subgroup_size[fis] ~ events$B_during_disp[fis]))


#plotting sub-group age and distance traveled during event
#combining the A group and the B group as it's not about the group id, but the average age (or number of adults) in those groups
plot(events$A_during_disp[fis], events$A_average_grp_age[fis], xlim = c(0,51), ylim = c(2,3), pch = 20, ylab = "Group average age class (1 = juvenile, 2 = subadult, 3 = adult)", xlab = "Distance traveled during fission event (m)")
points(events$B_during_disp[fis], events$B_average_grp_age[fis], pch = 20)


#plot distance traveled depending on number of adults in subgroup
plot(events$A_during_disp[fis], events$A_n_adults[fis], pch = 20, xlab = "Distance traveled during fission event (m)", ylab = "Number of adults in subgroup")
points(events$B_during_disp[fis], events$B_n_adults[fis], pch = 20)

#plot distance traveled depending on proportion of adults in subgroup
plot(events$A_during_disp[fis], events$A_proportion_adults[fis], pch = 20)
points(events$B_during_disp[fis], events$B_proportion_adults[fis], pch = 20)

#plot distance traveled depending on number of subadults in subgroup
plot(events$A_during_disp[fis], events$A_n_subadults[fis], pch = 20)
points(events$B_during_disp[fis], events$B_n_subadults[fis], pch = 20)

#plot distance traveled depending on number of juveniles in subgroup
plot(events$A_during_disp[fis], events$A_n_juveniles[fis], pch = 20)
points(events$B_during_disp[fis], events$B_n_juveniles[fis], pch = 20)

#plot distance traveled depending on number of individuals in subgroup
plot(events$A_during_disp[fis],events$A_subgroup_size[fis], pch = 20, ylim = c(0, 10))
points(events$B_during_disp[fis], events$B_subgroup_size[fis], pch = 20)





### FUSION PLOTTING ###

fus <- events$event_type == "fusion"

#need to think about what sort of plots to do...






