#Identify splits and merges using algorith from Della Libera et al. submitted

library(lubridate)
library(dbscan)

#DIRECTORIES AND PARAMETERS
R_inner <- 15
R_outer <- 50

codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'galaxy' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

#SOURCE FUNCTIONS
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

ff_data <- detect_fissions_and_fusions(R_inner, R_outer, xs, ys, ts, coati_ids)

analyse_ff_event(11, events = ff_data$events_detected, xs, ys, max_time = 600)

#animate events
i <- 10
max_time <- 600
step <- 5

ti <- events_detected$tidx[i] - max_time
tf <- events_detected$tidx[i] + max_time
group_A_idxs <- events_detected$group_A_idxs[i][[1]]
group_B_idxs <- events_detected$group_B_idxs[i][[1]]
all_idxs <- 1:nrow(xs)
not_involved_idxs <- setdiff(all_idxs, c(group_A_idxs, group_B_idxs))
event_type <- events_detected$event_type[i]
xmin <- min(xs[,ti:tf],na.rm=T)
xmax <- max(xs[,ti:tf],na.rm=T)
ymin <- min(ys[,ti:tf],na.rm=T)
ymax <- max(ys[,ti:tf],na.rm=T)
quartz()
tseq <- seq(ti,tf,step)
for(t in tseq){
  #clear()
  if(t >= events_detected$tidx[i]){
    points((xmax+xmin)/2,(ymin+ymax)/2,col='orange',cex = 5)
  }
  plot(NULL, xlim=c(xmin,xmax),ylim=c(ymin,ymax),asp=1,xlab='Easting',ylab='Northing',main = event_type)
  if(length(not_involved_idxs)>0){
    points(xs[not_involved_idxs,t],ys[not_involved_idxs,t],col='#00000033',pch=19)
  }
  points(xs[group_A_idxs,t],ys[group_A_idxs,t],pch=19, col = 'red')
  points(xs[group_B_idxs,t],ys[group_B_idxs,t],pch=19, col = 'blue')
  Sys.sleep(.01)
}
