#This script is for looking at how call rates relate to fission-fusion dynamics in coatis

library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)


#directory holding all the data
#datadir <- '~/Dropbox/coatithon/calling_during_fissions_and_fusions/data'
datadir <- "C:/Users/egrout/Dropbox/coatithon/calling_during_fissions_and_fusions/data"
callfile <- 'all_data_hms_synched.csv'
ff_file <- 'galaxy_auto_ff_events_characterized.RData'
gps_file <- 'galaxy_xy_highres_level1.RData'
id_file <- 'galaxy_coati_ids.RData'

#LOAD DATA
setwd(datadir)
calls <- read.csv(callfile, header=T, sep = ',')
load(id_file)
load(gps_file)
rm("xs","ys")  #only need the time indices

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

#Load fission-fusion events
load(ff_file)
group_events_data <- events
rm('events')

#Ari's script for creating the df with call rates before, during, and after event. Cini's version starts below commented code

# #create a new data frame with data at the individual level
# ind_events_data <- data.frame()
# for(i in 1:nrow(group_events_data)){
#   event_idx <- group_events_data$event_idx[i]
#   ind_idxs <- 1:nrow(coati_ids)
#   rows <- data.frame(event_idx = rep(event_idx, length(ind_idxs)),
#                      ind_idx = ind_idxs)
#   rows_all <- rbind(rows, rows, rows)
#   rows_all$period <- c(rep('before',length(ind_idxs)), 
#                        rep('during',length(ind_idxs)),
#                        rep('after',length(ind_idxs)))
#   ind_events_data <- rbind(ind_events_data, rows_all)
# }
# 
#
# #compute the call rate during an event period for a given individual
# get_call_rates_for_event <- function(ind_idx, event_idx, period, group_events_data_tmp = group_events_data, labeled_periods_tmp = labeled_periods){
#   event_data_curr <- group_events_data_tmp[which(group_events_data$event_idx == event_idx),]
#   before_time <- ts[event_data_curr$before_time]
#   start_time <- ts[event_data_curr$start_time]
#   end_time <- ts[event_data_curr$end_time]
#   after_time <- ts[event_data_curr$after_time]
#   
#   if(period == 'before'){
#     t0 <- before_time
#     tf <- start_time
#   } 
#   if(period == 'during'){
#     t0 <- start_time
#     tf <- end_time
#   }
#   if(period == 'after'){
#     t0 <- end_time
#     tf <- after_time
#   }
#   
#   #duration of event
#   dur <- as.numeric(difftime(tf, t0, units = 'secs'))
#   
#   #check that t0 and tf exist, if not return NA
#   if(is.na(t0) | is.na(tf)){
#     out <- list()
#     out$calls_contact <- NA
#     out$calls_agg <- NA
#     out$duration <- dur
#     return(out)
#   }
#   #check whether that individual has calls labeled within that time window, if not, return NA
#   labeled_periods_ind <- labeled_periods_tmp[which(labeled_periods_tmp$ind_idx == ind_idx),]
#   if(sum(labeled_periods_ind$starttime <= t0 & labeled_periods_ind$stoptime >= tf) == 0){
#     out <- list()
#     out$calls_contact <- NA
#     out$calls_agg <- NA
#     out$duration <- dur
#     return(out)
#   } 
#   
#   #get call rates
#   calls_agg <- sum(calls$calltype == 'aggression call' & 
#                      calls$name == coati_ids$name[ind_idx] &
#                      calls$datetime_synch >= t0 &
#                      calls$datetime_synch < tf, na.rm=T)
#   calls_contact <- sum(calls$calltype == 'contact call' & 
#                          calls$name == coati_ids$name[ind_idx] &
#                          calls$datetime_synch >= t0 &
#                          calls$datetime_synch < tf, na.rm=T)
#   
#   out <- list()
#   out$calls_contact <- calls_contact
#   out$calls_agg <- calls_agg
#   out$duration <- dur
#   
#   return(out)
# }
# 
# ind_events_data$agg_calls <- ind_events_data$contact_calls <- ind_events_data$duration <- ind_events_data$event_type <- ind_events_data$subgroup <-  NA
# for(i in 1:nrow(ind_events_data)){
#   calls_out <- get_call_rates_for_event(ind_idx = ind_events_data$ind_idx[i], 
#                                         event_idx = ind_events_data$event_idx[i],
#                                         period = ind_events_data$period[i],
#                                         group_events_data_tmp = group_events_data,
#                                         labeled_periods_tmp = labeled_periods)
#   ind_events_data$contact_calls[i] <- calls_out$calls_contact
#   ind_events_data$agg_calls[i] <- calls_out$calls_agg
#   ind_events_data$duration[i] <- calls_out$duration
#   
#   #info from the group table to transfer to the individual table
#   group_events_data_tmp <- group_events_data[which(group_events_data$event_idx==ind_events_data$event_idx[i]),]
#   
#   #event type (fission or fusion)
#   ind_events_data$event_type[i] <- group_events_data_tmp$event_type
#   
#   #subgroup 
#   subA <- group_events_data_tmp$group_A_idxs[[1]]
#   subB <- group_events_data_tmp$group_B_idxs[[1]]
#   ind_events_data$subgroup[i] <- NA
#   if(ind_events_data$ind_idx[i] %in% subA){
#     ind_events_data$subgroup[i] <- 'A'
#   } 
#   if(ind_events_data$ind_idx[i] %in% subB){
#     ind_events_data$subgroup[i] <- 'B'
#   }
#   
# }
# 
# #get call rates by dividing by duration
# ind_events_data$agg_call_rate <- ind_events_data$agg_calls / ind_events_data$duration
# ind_events_data$contact_call_rate <- ind_events_data$contact_calls / ind_events_data$duration


#Cini's script for creating the df with call rates before, during, and after event.

# get times to posix format and to UTC time zone
calls$datetime_synch_pos<-as.POSIXct(calls$datetime_synch, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 
labeled_periods$tstart<-as.POSIXct(labeled_periods$starttime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 
labeled_periods$tstop<-as.POSIXct(labeled_periods$stoptime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") 

# get the call rate for each individual for each event before, during and after the event (also some additional info, like the subgroup id and the distance moved)
ind_events_data<-data.frame()
head(group_events_data)
for(i in 1:nrow(group_events_data)){
  event.times<-group_events_data[i,c(10:14)]
  times<-ts[as.numeric(event.times[,c(2:5)])]  # start - end backwards
  date<-substring(times[1], 1,10) # date of the event for a quick check of label presence
  
  if(all(!is.na(times))){

    # find individuals where labels are available for duration of event.
    l.date<-labeled_periods[which(substring(labeled_periods$starttime,1,10) == date),] 
    
    # there are quite a few dates with no labeled sound data, so ignore those
    if(nrow(l.date) >0){
      l.inds<-l.date[which(l.date$tstart <= times[4] & l.date$tstop >= times[1]),"ind_idx"]
      
      # only use events were at least one individual has labeled sound data
      if(length(l.inds) >0){
        # get all calls during event periods
        before<-calls[which(calls$datetime_synch_pos >= times[4] & calls$datetime_synch_pos < times[3]),]
        during<-calls[which(calls$datetime_synch_pos >= times[3] & calls$datetime_synch_pos < times[2]),]
        after<-calls[which(calls$datetime_synch_pos >= times[2] & calls$datetime_synch_pos < times[1]),]
        
        for(j in l.inds){
          event<-group_events_data[i,c("event_idx","event_type")] 
          event<-rbind(event,event,event) # as there are three time periods I will rbind the events 3 times
          event$period<-c("before","during","after")
          event$duration<-c(as.numeric(difftime(times[3],times[4], units = "secs")),
                            as.numeric(difftime(times[2],times[3], units = "secs")),
                            as.numeric(difftime(times[1],times[2], units = "secs")))
          event$ind_idx <- j
          
          # get the call number for each call type (contact/ agression) for each period
          # before
          event[event$period == "before","contact_calls"] <- nrow(before[which(before$name == coati_ids[j,"name"] & before$calltype == "contact call"),])
          event[event$period == "before","agg_calls"] <-nrow(before[which(before$name == coati_ids[j,"name"] & before$calltype == "aggression call"),])
          
          # during
          event[event$period == "during","contact_calls"] <-nrow(during[which(during$name == coati_ids[j,"name"] & during$calltype == "contact call"),])
          event[event$period == "during","agg_calls"] <-nrow(during[which(during$name == coati_ids[j,"name"] & during$calltype == "aggression call"),])
          
          #after
          event[event$period == "after","contact_calls"] <-nrow(after[which(after$name == coati_ids[j,"name"] & after$calltype == "contact call"),])
          event[event$period == "after","agg_calls"] <-nrow(after[which(after$name == coati_ids[j,"name"] & after$calltype == "aggression call"),])
          
          # get the subgroup for each ID
          a_inds<-unlist(group_events_data[group_events_data$event_idx == unique(event$event_idx),"group_A_idxs"])
          b_inds<-unlist(group_events_data[group_events_data$event_idx == unique(event$event_idx),"group_B_idxs"])
          event$subgroup<-ifelse(j %in% a_inds, "A",ifelse(j %in% b_inds,"B","NA"))
          
          # # Optional: get distance the subgroup the individual is in moved during the event
          # a_move<-group_events_data[group_events_data$event_idx == unique(event$event_idx),"A_during_disp"]
          # b_move<-group_events_data[group_events_data$event_idx == unique(event$event_idx),"B_during_disp"]
          # 
          # event$sub_move_dist_during<-ifelse(j %in% a_inds,a_move,b_move)
          
          event$n_ind_labeled<-length(l.inds) # total number of individuals that have labeled data for this event
          ind_events_data<-rbind( ind_events_data,event)
        }
      }
    }
  }
}
rm(list = setdiff(ls(),c("ind_events_data","calls","group_events_data")))
ind_events_data$agg_call_rate<-ind_events_data$agg_calls / ind_events_data$duration
ind_events_data$contact_call_rate<-ind_events_data$contact_calls / ind_events_data$duration
#TODO check events where the before and start time are the same

leadership_metric<-c("position","crosstime","crosstime_ownfinishline")
n = 2
filename<-paste0("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_LeaderRank_",leadership_metric[n],".RData")

load(filename) 


fission_leaders_rank <- as.data.frame(out$fission_leaders) 
fission_leaders_rank<-as.data.frame(t(fission_leaders_rank)) 
names(fission_leaders_rank)<-1:length(fission_leaders_rank) 
fission_leaders_rank$event_idx<-1:nrow(fission_leaders_rank)
fission_leaders_rank<-melt(fission_leaders_rank, id.vars = "event_idx")
colnames(fission_leaders_rank) <- c("event_idx", "ind_idx", "leader_rank") 
fission_leaders_rank$event_type<-"fission"
fission_leaders_rank<-fission_leaders_rank[!is.na(fission_leaders_rank$leader_rank),]

fusion_leaders_rank <- as.data.frame(out$fusion_leaders) 
fusion_leaders_rank<-as.data.frame(t(fusion_leaders_rank)) 
names(fusion_leaders_rank)<-1:length(fusion_leaders_rank) 
fusion_leaders_rank$event_idx<-1:nrow(fusion_leaders_rank)
fusion_leaders_rank<-melt(fusion_leaders_rank, id.vars = "event_idx")
colnames(fusion_leaders_rank) <- c("event_idx", "ind_idx", "leader_rank") 
fusion_leaders_rank$event_type<-"fusion"
fusion_leaders_rank<-fusion_leaders_rank[!is.na(fusion_leaders_rank$leader_rank),]


leader_ranks<-rbind(fission_leaders_rank,fusion_leaders_rank)

ind_events_data_merge <- merge(ind_events_data, leader_ranks, by = c("event_idx","event_type","ind_idx"), all.x = T) 


#change factor levels so before is shown before after in plot 
ind_events_data_merge$period <- factor(ind_events_data_merge$period, levels = c("before","during" ,"after")) 

#removing events where there more than 2 individuals are in the event 
ind_events_data_merge <- ind_events_data_merge[ind_events_data_merge$n_ind_labeled >2, ] 

#rounding the leader rank for plotting 
ind_events_data_merge$rounded_leader_rank <- round(ind_events_data_merge$leader_rank, 1) 

event_type<-"fission"
ggplot(data = ind_events_data_merge[ind_events_data_merge$event_type == event_type,],  
       aes(x = period, y = agg_call_rate, col = as.factor(rounded_leader_rank), group = as.factor(ind_idx)))+ 
  geom_line(size = 1.5)+
  labs(color = "Leader rank")+
  ggtitle(paste(event_type, "aggression",leadership_metric[n]))+
  facet_wrap(~event_idx) 

ggplot(data = ind_events_data_merge[ind_events_data_merge$event_type == event_type,],  
       aes(x = period, y = contact_call_rate, col = as.factor(rounded_leader_rank), group = as.factor(ind_idx)))+ 
  geom_line(size = 1.5)+
  labs(color = "Leader rank")+
  ggtitle(paste(event_type, "contact",leadership_metric[n]))+
  facet_wrap(~event_idx) 



#TODO get the lm line to work!! Add color for moving group
event_type<-"fission"
ggplot(data = ind_events_data_merge[ind_events_data_merge$event_type == event_type,],  
       aes(x = leader_rank, y = agg_call_rate, group = as.factor(ind_idx)))+ 
  geom_point(size = 1.5)+
  # geom_line(aes(contact_call_rate~leader_rank))+
  labs(color = "Leader rank")+
  ggtitle(paste(event_type, "aggression",leadership_metric[n]))+
  facet_wrap(~period) 



## figure out a way to define fission types
#look at the difference in travel distance 

# event<-data.frame()
# for(i in 1:nrow(group_events_data)){
#   event.times<-group_events_data[i,c(10:14)]
#   times<-ts[as.numeric(event.times[,c(2:5)])]  # start - end backwards
#   event_<-group_events_data[i,c("event_idx","event_type")]
#   event_$b.duration<-as.numeric(difftime(times[3],times[4], units = "secs"))
#   event_$d.duration<-as.numeric(difftime(times[2],times[3], units = "secs"))
#   event<-rbind(event,event_)
# }
# 
# dist_travel_df <- group_events_data[, c("event_type", "event_idx", "B_during_disp", "A_during_disp", "AB_before_disp" )]
# 
# dist_travel_df<-merge(dist_travel_df, event, by = c("event_idx","event_type"), all = T)
# 
# fission_dist <- dist_travel_df[dist_travel_df$event_type == "fission",]
# 
# fission_dist <- fission_dist[!is.na(fission_dist$B_during_disp),]
# 
# fission_dist$move_diff <- abs(fission_dist$B_during_disp - fission_dist$A_during_disp)
# fission_dist$B_move_diff <- fission_dist$B_during_disp - fission_dist$AB_before_disp
# fission_dist$A_move_diff <- fission_dist$A_during_disp - fission_dist$AB_before_disp
# 
# fission_dist$B_speed <- fission_dist$B_during_disp / fission_dist$d.duration
# fission_dist$A_speed <- fission_dist$A_during_disp / fission_dist$d.duration
# fission_dist$AB_speed <- fission_dist$AB_before_disp / fission_dist$b.duration
# 
# hist(fission_dist$A_speed, breaks = 20)
# 
# #for error checking, should remove the fissions where speed is incredibly high (e.g. more than 2.5m/second) in a short duration (6 seconds)
# #for now, we remove these events here
# fission_dist <- fission_dist[fission_dist$B_speed < 2.5,]
# fission_dist <- fission_dist[fission_dist$A_speed < 2.5,]
# 
# #removing rows where the before speed is 0
# fission_dist <- fission_dist[which(fission_dist$b.duration > 0),]
# 
# #calculate speed diff
# fission_dist$speed_diff <- abs(fission_dist$B_speed - fission_dist$A_speed)
# fission_dist$B_speed_diff <- fission_dist$B_speed - fission_dist$AB_speed
# fission_dist$A_speed_diff <- fission_dist$A_speed - fission_dist$AB_speed
# 
# plot(fission_dist$AB_speed*60, fission_dist$speed_diff)
# 
# 
# 
# plot(fission_dist$AB_before_disp, fission_dist$move_diff)
# points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff >20 & fission_dist$AB_before_disp > 25,], col = "steelblue2", pch = 19)
# points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff >20 & fission_dist$AB_before_disp <= 25,], col = "indianred", pch = 19)
# points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff <20 & fission_dist$AB_before_disp <= 25,], col = "gold2", pch = 19)
# 
# 
# 
# ind_events_data_long <- ind_events_data %>%
#   pivot_longer(c(agg_call_rate, contact_call_rate), names_to = "call", values_to = "rate")
# 
# #change factor levels so before is shown before after in plot
# ind_events_data_long$period <- factor(ind_events_data_long$period, levels = c("before","during" ,"after"))
# 
# #get the points to correspond to the period
# 
# ggplot(data = ind_events_data_long[ind_events_data_long$event_type == "fission",], 
#             aes(x = call, y = rate, fill = period))+
#   geom_boxplot(outlier.shape = NA)+ 
#   geom_jitter(position = position_dodge(width = 0.75), size = 0.5, color = "gray3", aes(group = interaction(call, period))) +
#   ylim(c(0,0.75))+
#   
#   #scale_fill_manual(values=c("indianred1", "indianred3", "indianred4"))+
#   #theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor =     element_blank(), panel.background = element_blank(), axis.text=element_text(size=14), axis.title = element_text(size = 14), legend.text = element_text(size = 14))+ 
#   xlab(" ") +
#   ylab("Call rate (per minute)")+
#   scale_alpha(guide = 'none')+
#   #facet_wrap(~event_idx)
# 
#   NULL








