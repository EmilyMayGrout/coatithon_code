#fission fusion analysis script for Galaxy group 
# We try to find times when the group is all together and coati calls were labeled and extract 5min call periods for each group member as a baseline call rate. Best across multiple days to correct for differences in "mood"? of each coati

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
gps_file <- "galaxy_xy_highres_level1.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#list of Rs
R <- 50

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
#getting subgorup data
subgroup_data <- get_subgroup_data(xs, ys, R)


#finding times when the number of subgroups is 1
together_times<-as.data.frame(ts[which(subgroup_data$n_subgroups ==1)]); names(together_times) = "time"
together_times$tdiff<-c(NA,difftime(together_times$time[2:nrow(together_times)], together_times$time[1:(nrow(together_times)-1)], units = "secs"))

# find "bouts" of group being together
n = 1
for(i in 2:nrow(together_times)){
  if(together_times$tdiff[i]>1){ n = n + 1}
  together_times$bout[i]<-n
}

# now get the start and stop time of each together bout
tts<-split(together_times,together_times$bout)
ttsl<-lapply(tts, function(tts_){
  tts_ss<-tts_[c(1, nrow(tts_)),"time"]
  
  return(tts_ss)
})



#get the call files 
datadir <- "C:/Users/egrout/Dropbox/coatithon/calling_during_fissions_and_fusions/data"
callfile <- 'all_data_hms_synched.csv'
id_file <- 'galaxy_coati_ids.RData'

#LOAD DATA
setwd(datadir)
calls <- read.csv(callfile, header=T, sep = ',')
load(id_file)

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

calling_together<-calls[which(calls$datetime_synch_pos %in% together_times$time),]





# calculate call rate for each together bout
# get the duration of the bout

# together_bout<-ttsl[[34]]
crt_<-as.data.frame(coati_ids[,1]); names(crt_)<-"name"
crt_$ind_idx <- 1:nrow(crt_)
call_rates_together<-data.frame()
n = 1
for(together_bout in ttsl){
  ct_<-calling_together[which(calling_together$datetime_synch_pos >= together_bout[1] & calling_together$datetime_synch_pos <= together_bout[2]), ]
  ct_<-ct_[!is.na(ct_$calltype),]
  
  if(nrow(ct_)> 0){
    calls_ind<-as.data.frame(table(ct_$name, ct_$calltype)); names(calls_ind)<-c("name","call_type","count")
    calls_all<-merge(calls_ind,crt_,by = "name")
    calls_all$bout_dur<-as.numeric(difftime(together_bout[2], together_bout[1], units = "secs"))
    calls_all$call_rate<-calls_all$count / calls_all$bout_dur
    calls_all$n_ind_in_bout<-length(unique(calls_all$name))
    calls_all$bout <-n
    call_rates_together<-rbind(call_rates_together, calls_all)
    
    n = n+1
  }
}

par(mar = c(14,4,1,1))
boxplot(call_rate ~ name*call_type, data = call_rates_together[call_rates_together$n_ind_in_bout > 5,], las = 2, xlab = "")

#could make aggression calls binary




