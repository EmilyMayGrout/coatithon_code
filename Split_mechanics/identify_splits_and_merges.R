#Identify splits and merges using algorith from Della Libera et al. submitted

library(lubridate)
library(dbscan)

#DIRECTORIES AND PARAMETERS
R_inner <- 15
R_outer <- 50

codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'presedente' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

#for Emily:
codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
#groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level2/'


#SOURCE FUNCTIONS
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#navigate into directory
setwd(codedir)

#read in events
events <- read.csv(paste0('C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/',group,'_manual_split_merge_clean.csv'), sep=';')

#read in coati ids
setwd(groupdir)
load(file=paste0(group,'_coati_ids.RData'))

#read in timestamp data - notice the level number, if 2, then its not got GPS speed errors
load(file=paste0(group,'_xy_highres_level2.RData'))

ff_data_50 <- detect_fissions_and_fusions(R_inner = 15, R_outer = 50, xs, ys, ts, coati_ids)

#saving the events_detected dataframe of the ff_data_50 list
if(group == 'galaxy'){
gal_events_detected <- ff_data_50$events_detected
save(gal_events_detected, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/gal_events_detected.Rda") 
}else if(group == 'presedente'){
  pres_events_detected <- ff_data_50$events_detected
      save(pres_events_detected, file = "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/pres_events_detected.Rda") }


#looking at one event
analyse_ff_event(12, events = ff_data_50$events_detected, xs, ys, ts, max_time = 600)

#animate events
i <- 10
max_time <- 600
step_t <- 5

#change gal to pres if presidente group
ti <- gal_events_detected$tidx[i] - max_time
tf <- gal_events_detected$tidx[i] + max_time
group_A_idxs <- gal_events_detected$group_A_idxs[i][[1]]
group_B_idxs <- gal_events_detected$group_B_idxs[i][[1]]
all_idxs <- 1:nrow(xs)
not_involved_idxs <- setdiff(all_idxs, c(group_A_idxs, group_B_idxs))
event_type <- gal_events_detected$event_type[i]
xmin <- min(xs[,ti:tf],na.rm=T)
xmax <- max(xs[,ti:tf],na.rm=T)
ymin <- min(ys[,ti:tf],na.rm=T)
ymax <- max(ys[,ti:tf],na.rm=T)
#quartz()
# tseq <- seq(ti,tf,step_t)
# 
# for(t in tseq){
#   #clear()
#   if(t >= gal_events_detected$tidx[i]){
#     points((xmax+xmin)/2,(ymin+ymax)/2,col='orange',cex = 5)
#   }
#   plot(NULL, xlim=c(xmin,xmax),ylim=c(ymin,ymax),asp=1,xlab='Easting',ylab='Northing',main = event_type)
#   if(length(not_involved_idxs)>0){
#     points(xs[not_involved_idxs,t],ys[not_involved_idxs,t],col='#00000033',pch=19)
#   }
#   points(xs[group_A_idxs,t],ys[group_A_idxs,t],pch=19, col = 'red')
#   points(xs[group_B_idxs,t],ys[group_B_idxs,t],pch=19, col = 'blue')
#   Sys.sleep(.01)
# }







