#What are the relevant spatial scales associated with coati groups?

#LIBRARY
library(lubridate)

#set time zone to UTC to avoid confusing time zone issues
Sys.setenv(TZ='UTC')

#FUNCTIONS
source('~/Dropbox/code_ari/move_comm_analysis/audio_gps_processing/spatially_discretized_headings.R')

#DIRECTORIES AND PARAMETERS
codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'galaxy' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

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

#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})


#HEADING CORRELATION VS DYADIC DISTANCE

n_inds <- nrow(xs)
n_times <- ncol(xs)
heads <- speeds <- matrix(NA, nrow = n_inds, ncol = n_times)
#day start indexes
days <- date(ts)
day_start_idxs <- c(1, which(diff(days)==1)+1)

#spatially discretized headings of individuals and speeds
speed_dt <- 60
for(i in 1:n_inds){
  for(d in 1:(length(day_start_idxs)-1)){
    tday <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
    headsday <- spatial.headings(xs[i,tday],ys[i,tday],10)
    heads[i,tday] <- headsday
    idxs_fut <- tday[seq(speed_dt+1,length(tday),1)]
    idxs_now <- tday[seq(1,length(tday)-speed_dt,1)]
    speeds[i,idxs_now] <- sqrt((xs[i, idxs_fut] - xs[i, idxs_now])^2 + (ys[i, idxs_fut] - ys[i, idxs_now])^2)
  }
}

#dyadic distances and heading correlations
dyad_dists <- dyad_dist_changes <- head_corrs <- speed_diffs <- log_speed_diffs <- array(NA, dim = c(n_inds,n_inds,n_times))
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    dyad_dists[i,j,] <- sqrt((xs[i,]-xs[j,])^2 + (ys[i,]-ys[j,])^2)
    head_corrs[i,j,] <- cos(heads[i,])*cos(heads[j,]) + sin(heads[i,])*sin(heads[j,])
    log_speed_diffs[i,j,] <- abs(log(speeds[i,]) - log(speeds[j,]))
    speed_diffs[i,j,] <- abs(speeds[i,] - speeds[j,])
    for(d in 1:(length(day_start_idxs)-1)){
      tday <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
      idxs_fut <- tday[seq(speed_dt+1,length(tday),1)]
      idxs_now <- tday[seq(1,length(tday)-speed_dt,1)]
      dyad_dist_changes[i,j,idxs_now] <- dyad_dists[i,j,idxs_fut] - dyad_dists[i,j,idxs_now]
    }
  }
}

dist_bins <- c(0, 10^seq(0,2.6,.1))
heads_mean <- speeds_mean <- log_speeds_mean <- change_dyad_dist_mean <- rep(NA, length(dist_bins)-1)
for(i in 1:(length(dist_bins)-1)){
  idxs <- which(dyad_dists >= dist_bins[i] & dyad_dists < dist_bins[i+1])
  heads_mean[i] <- mean(head_corrs[idxs], na.rm=T)
  speeds_mean[i] <- mean(speed_diffs[idxs], na.rm=T)
  log_speeds_mean[i] <- mean(log_speed_diffs[idxs], na.rm=T)
  change_dyad_dist_mean[i] <- mean(dyad_dist_changes[idxs], na.rm=T)
}

#PLOTS

#Plot 1: Distribution of log(dyadic distances) between coatis
quartz()
bins <- seq(-2,3,length=30)
histo <- hist(log(dyad_dists, 10), plot=T, breaks=30, xlab = 'Log dyadic distance (m)',main='')
abline(v=log(50,10), col ='red', lwd = 3, lty=2)
text(y = max(histo$counts)*.95, x = log(50,10)+.4, labels = '50 m', col = 'red',cex=1.5)

#Plot 2: Mean heading correlation as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],heads_mean, xlab = 'Distance apart (m)', ylab = 'Mean heading correlation', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,1))
abline(h=0, lty = 2)
abline(v=50, col = 'red', lwd = 3, lty = 2)

#Plot 3: Mean difference in speed as a function of distance 
quartz()
plot(dist_bins[2:length(dist_bins)],speeds_mean, xlab = 'Distance apart (m)', ylab = 'Mean absolute speed difference (m/min)', pch = 19, col = '#00000066', cex = 1.5, ylim = c(0,max(speeds_mean)))
abline(h=0, lty = 2)
abline(v=50, col = 'red', lwd = 3, lty = 2)

#Plot 4: Mean change in dyadic distance
quartz()
plot(dist_bins[2:length(dist_bins)],change_dyad_dist_mean, xlab = 'Distance apart (m)', ylab = 'Mean change in dyadic distance per minute (m)', pch = 19, col = '#00000066', cex = 1.5, ylim = c(min(change_dyad_dist_mean),max(change_dyad_dist_mean)))
abline(h=0, lty = 2)
abline(v=50, col = 'red', lwd = 3, lty = 2)
