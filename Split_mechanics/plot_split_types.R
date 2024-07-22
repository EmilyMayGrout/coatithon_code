#this script is making plots of the different split types for the paper


#LIBRARY
library(lubridate)
library(scales)
library(ggplot2)
library(patchwork)
library(dplyr)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#DIRECTORIES AND PARAMETERS

#who is using (ari or emily)
user <- 'emily'

#which group - galaxy or presedente
group <- 'galaxy' #subdirectory where the group data is stored

#whether to identify splits and merges automatically (if F) or use manually identified events (if T)
use_manual_events <- F

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
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
    plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
  }
}

#FUNCTIONS
#read in functions
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(codedir)


if(use_manual_events){
  events <- read.csv(paste0('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/',group,'_manual_split_merge_clean2.csv'),  sep=",", header=TRUE)
  #preprocess events to...
  events <- events[which(events$fission_time!='before start'),] #remove events where we missed the start
  events <- events[which(events$event_type %in% c('fission','fusion')),] #only include fission and fusion events (remove 'almost fusion')
  
  
} else{ #otherwise load in automated events
  #read in automated events - df made in characterize_splits_amd_merges 
  load(paste0(groupdir, group,'_auto_ff_events_characterized.RData'))
  
}

#read in coati ids
setwd(groupdir)

#read in timestamp data
load(file=paste0(group,'_xy_highres_level2.RData'))

load(file=paste0(group,'_coati_ids.RData'))
#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})


#plot different event types - code copied from analyse_ff_event script


fusions <- events[events$event_type == "fusion",]

for (i in 1:nrow(fusions)){
analyse_ff_event(i, events = fusions, xs, ys, ts, max_time = 900, plot = T)
}

dev.off()

#add index to the df to find events easier
index <- 1:nrow(events)
events <- cbind(index, events)

#Presidente fissions
#02.02 11:58:14
#Pres fusions
#02.01 12:12:00
#02-01 13:07:49 
#01-28 12:54:59
# i = 300, 178, 235, 139

#for Galaxy fissions
#01-09 11:33:28
#01-07 11:27:13
#01-06 11:55:29
#01-05 12:09:12
#01-03 12:48:47
#01-05 11:32:35
#01-07 13:24:12
#12-28 11:57:12 - nice fission

# i = 130, 108, 94, 76, 61, 73,112, 27

#Galaxy fusions
#01-06 11:58:34
#01-03 13:49:40

#i = 64, 59

i = 139
#events[i,]
max_time = 850
thresh_h = 50
thresh_l = 15
time_window = 120

t_event <- events$tidx[i] #time of the event
group_A <- events$group_A_idxs[i][[1]] #group A individual idxs
group_B <- events$group_B_idxs[i][[1]] #group B individual idxs
group_A_names <- events$group_A[i]
group_B_names <- events$group_B[i]
ti <- t_event - max_time #initial time to plot
tf <- t_event + max_time #final time to plot

if(ti <1 | tf > length(ts)){
  out <- list()
  out$start_time <- out$end_time <- out$before_time <- out$after_time <- NA
  out$disps <- out$speeds <- NULL
  out$turn_angle_A <- out$turn_angle_B <- out$split_angle <- NA
  return(out)
}

event_type <- events$event_type[i]
datetime <- events$datetime[i]
nA <- length(group_A)
nB <- length(group_B)
nT <- ncol(xs)

#get x and y coordinates of the relevant individuals in each subgroup
xA <- matrix(xs[group_A,],nrow=nA,ncol=nT)
xB <- matrix(xs[group_B,],nrow=nB,ncol=nT)
yA <- matrix(ys[group_A,],nrow=nA,ncol=nT)
yB <- matrix(ys[group_B,],nrow=nB,ncol=nT)

#get centroids of the two subgroups
xcA <- colMeans(xA, na.rm=T)
ycA <- colMeans(yA, na.rm=T)
xcB <- colMeans(xB, na.rm=T)
ycB <- colMeans(yB, na.rm=T)

#get distance between centroids
dyad_dist <- sqrt((xcA - xcB)^2 + (ycA - ycB)^2)

#classify the dyadic distance into categories:
#0 = below lower thershold
#1 = between thresholds
#2 = above higher threshold
dyad_dist_event <- dyad_dist[ti:tf]

#first consider modifying thresholds according to subtlety 1 above
upper <- thresh_h
lower <- thresh_l
after_idxs <- (max_time+1):(2*max_time+1) #indexes after the marked event
middle_idxs <- (max_time / 2):(max_time*3/2)
before_idxs <- 1:max_time #indexes before the marked event
if(event_type == 'fission'){
  if(sum(!is.na(dyad_dist_event[after_idxs]))>1){
    if(max(dyad_dist_event[after_idxs],na.rm=T) < thresh_h){
      upper <- max(dyad_dist_event[after_idxs],na.rm=T) - .001
    } else{
      upper <- thresh_h
    }
  }
  if(sum(!is.na(dyad_dist_event[middle_idxs]))){
    if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
      lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
    } else{
      lower <- thresh_l
    }
  }
}
if(event_type == 'fusion'){
  if(sum(!is.na(dyad_dist_event[before_idxs]))>1){
    if(max(dyad_dist_event[before_idxs],na.rm=T) < thresh_h){
      upper <- max(dyad_dist_event[before_idxs],na.rm=T) - .001
    } else{
      upper <- thresh_h
    }
  }
  if(sum(!is.na(dyad_dist_event[middle_idxs]))>1){
    if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
      lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
    } else{
      lower <- thresh_l
    }
  }
}

#if the upper bound was changed to something < thresh_l, move threshold back to thresh_h
if(upper <= thresh_l){
  upper <- thresh_h
}
#likewise for lower bound
if(lower >= thresh_h){
  lower <- thresh_l
}

#get category of each moment in time
#0 = below lower, 1 = middle, 2 = above upper
category <- rep(NA, length(dyad_dist_event))
category[which(dyad_dist_event < lower)] <- 0
category[which(dyad_dist_event >= lower & dyad_dist_event < upper)] <- 1
category[which(dyad_dist_event >= upper)] <- 2
category[which(is.na(dyad_dist_event))] <- 3 #NAs are denoted with 3

#run length encoding to get sequences of each category
seqs <- rle(category)

#find sequences of high-middle-low (2,1,0) or low-mid-high (0,1,2)
seqs_str <- paste0(as.character(seqs$values), collapse = '') #convert to string
if(event_type=='fission'){
  event_loc <- as.data.frame(str_locate_all(seqs_str,'012')[[1]])
}
if(event_type == 'fusion'){
  event_loc <- as.data.frame(str_locate_all(seqs_str,'210')[[1]])
}

#for seuqneces of hml or lmh (for fission and fusion respectively), get the time index when they start and end
if(nrow(event_loc)>0){
  for(r in 1:nrow(event_loc)){
    event_loc$start_time[r] <- ti + sum(seqs$lengths[1:event_loc$start[r]])
    event_loc$end_time[r] <- ti + sum(seqs$lengths[1:event_loc$end[r]-1])
  }
}

# par(mfrow=c(2,1))
# plot(ti:tf, dyad_dist[ti:tf],type='l', main = paste(event_type, datetime),xlab='Time (min)',ylab = 'Distance apart (m)')
# abline(v=t_event,col='black', lty = 2)
# abline(h = thresh_h, col = 'darkorange1')
# abline(h = thresh_l, col = 'magenta')
# for(r in 1:nrow(event_loc)){
#   abline(v=event_loc$start_time[r], col = 'green')
#   abline(v=event_loc$end_time[r], col = 'green')
# }

xmin <- min(min(xA[,ti:tf],na.rm=T),min(xB[,ti:tf],na.rm=T))
xmax <- max(max(xA[,ti:tf],na.rm=T),max(xB[,ti:tf],na.rm=T))
ymin <- min(min(yA[,ti:tf],na.rm=T),min(yB[,ti:tf],na.rm=T))
ymax <- max(max(yA[,ti:tf],na.rm=T),max(yB[,ti:tf],na.rm=T))

xdistline_start <- (ceiling(xmin / 10) * 10)
xdistline_end <- xdistline_start + 50
yline <- (floor(ymin / 10) * 10)+10


#xA colors are: 
#'#CAFF7080' 'olivedrab3'
#xB colors are:
#'lightblue1' 'dodgerblue3'


#png(height = 800, width = 1000, units = 'px', res = 90, filename = paste0(plot_dir, group, "_event_",t_event, ".png"))
plot(NULL, xlim=c(xmin,xmax),ylim=c(ymin,ymax), asp=1, xlab='Easting', ylab = 'Northing',main = paste('(Green =', group_A_names, '), (Blue =', group_B_names,')'))
for(j in 1:nrow(xA)){
  lines(xA[j,ti:tf],yA[j,ti:tf],type='l',col='#CAFF7080')
}
for(j in 1:nrow(xB)){
  lines(xB[j,ti:tf],yB[j,ti:tf],type='l',col='thistle1')
}
lines(xcA[ti:tf],ycA[ti:tf], col = 'yellow3', lwd = 3)
lines(xcB[ti:tf],ycB[ti:tf], col = 'mediumpurple2', lwd = 3)

cex = 3

points(xcA[t_event], ycA[t_event], pch = 8, col = 'black', cex = cex) # star is the event
points(xcB[t_event], ycB[t_event], pch = 8, col = 'black', cex = cex)
points(xcA[ti], ycA[ti], pch = 1, col = 'black', cex = cex) #o is the start
points(xcB[ti], ycB[ti], pch = 1, col = 'black', cex = cex)
points(xcA[tf], ycA[tf], pch = 4, col = 'black', cex = cex) # x is the end
points(xcB[tf], ycB[tf], pch = 4, col = 'black', cex = cex)

lines(c(xdistline_start, xdistline_end), c(yline, yline), type='l', lwd = 3)
text((xdistline_start + 25), (yline + 6), "50 m", cex = 2)

dev.off()

#algorithm-identified start and end
# points(xcA[event_loc$start_time],ycA[event_loc$start_time],pch = 1, col = 'green')
# points(xcB[event_loc$start_time],ycB[event_loc$start_time],pch = 1, col = 'green')
# points(xcA[event_loc$end_time],ycA[event_loc$end_time],pch = 4, col = 'green')
# points(xcB[event_loc$end_time],ycB[event_loc$end_time],pch = 4, col = 'green')

