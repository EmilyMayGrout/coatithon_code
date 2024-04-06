#This script is for looking at how call rates relate to fission-fusion dynamics in coatis
#need to do this with level2 data

library(dplyr)
library(ggplot2)
library(tidyverse)
library(reshape2)

#directory holding all the data
#datadir <- '~/Dropbox/coatithon/calling_during_fissions_and_fusions/data'
datadir <- "C:/Users/egrout/Dropbox/coatithon/calling_during_fissions_and_fusions/data"
callfile <- 'all_data_hms_synched.csv'
ff_file <- 'galaxy_auto_ff_events_characterized.RData' #this has been rerun with level2 data
gps_file <- 'galaxy_xy_highres_level2.RData'
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

rm(list = setdiff(ls(),c("ind_events_data","calls","group_events_data", "ts")))

ind_events_data$agg_call_rate <- ind_events_data$agg_calls/ind_events_data$duration
ind_events_data$contact_call_rate <- ind_events_data$contact_calls/ind_events_data$duration
#TODO check events where the before and start time are the same - because overlap with previous event


leadership_metric <- c("position","crosstime","crosstime_ownfinishline")

n = 2
filename <- paste0("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_LeaderRank_",leadership_metric[n],".RData")
load(filename) 


fission_leaders_rank <- as.data.frame(out$fission_leaders) 
fission_leaders_rank <- as.data.frame(t(fission_leaders_rank)) 
names(fission_leaders_rank) <- 1:length(fission_leaders_rank) 
fission_leaders_rank$event_idx <- 1:nrow(fission_leaders_rank)
fission_leaders_rank <- melt(fission_leaders_rank, id.vars = "event_idx")
colnames(fission_leaders_rank) <- c("event_idx", "ind_idx", "leader_rank") 
fission_leaders_rank$event_type <- "fission"
fission_leaders_rank <- fission_leaders_rank[!is.na(fission_leaders_rank$leader_rank),]

fusion_leaders_rank <- as.data.frame(out$fusion_leaders) 
fusion_leaders_rank <- as.data.frame(t(fusion_leaders_rank)) 
names(fusion_leaders_rank) <- 1:length(fusion_leaders_rank) 
fusion_leaders_rank$event_idx <- 1:nrow(fusion_leaders_rank)
fusion_leaders_rank <- melt(fusion_leaders_rank, id.vars = "event_idx")
colnames(fusion_leaders_rank) <- c("event_idx", "ind_idx", "leader_rank") 
fusion_leaders_rank$event_type <- "fusion"
fusion_leaders_rank <- fusion_leaders_rank[!is.na(fusion_leaders_rank$leader_rank),]


leader_ranks <- rbind(fission_leaders_rank, fusion_leaders_rank)
ind_events_data_merge <- merge(ind_events_data, leader_ranks, by = c("event_idx","event_type","ind_idx"), all.x = T) 


#change factor levels so before is shown before after in plot 
ind_events_data_merge$period <- factor(ind_events_data_merge$period, levels = c("before","during" ,"after")) 

#removing events where there more than 2 individuals are in the event 
ind_events_data_merge <- ind_events_data_merge[ind_events_data_merge$n_ind_labeled > 2, ] 

#rounding the leader rank for plotting 
ind_events_data_merge$rounded_leader_rank <- round(ind_events_data_merge$leader_rank, 1) 

event_type <- "fission"
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
event_type <- "fission"
ggplot(data = ind_events_data_merge[ind_events_data_merge$event_type == event_type,],  
       aes(x = leader_rank, y = agg_call_rate, group = as.factor(ind_idx)))+ 
  geom_point(size = 1.5)+
  # geom_line(aes(contact_call_rate~leader_rank))+
  labs(color = "Leader rank")+
  ggtitle(paste(event_type, "aggression",leadership_metric[n]))+
  facet_wrap(~period) 

## figure out a way to define fission types
#look at the difference in travel distance 

event <- data.frame()
for(i in 1:nrow(group_events_data)){
  event.times <- group_events_data[i,c(10:14)]
  times <- ts[as.numeric(event.times[,c(2:5)])]  # start - end backwards
  event_ <- group_events_data[i,c("event_idx","event_type")]
  event_$b.duration <- as.numeric(difftime(times[3],times[4], units = "secs"))
  event_$d.duration <- as.numeric(difftime(times[2],times[3], units = "secs"))
  event <- rbind(event,event_)
}


dist_travel_df <- group_events_data[, c("event_type", "event_idx", "B_during_disp", "A_during_disp", "AB_before_disp" )]

dist_travel_df<-merge(dist_travel_df, event, by = c("event_idx","event_type"), all = T)

#---------FISSIONS------------------------------------------------------------------------------------------------------------

fission_dist <- dist_travel_df[dist_travel_df$event_type == "fission",]
fission_dist <- fission_dist[!is.na(fission_dist$B_during_disp),]

#get distance travelled
fission_dist$move_diff <- abs(fission_dist$B_during_disp - fission_dist$A_during_disp)
fission_dist$B_move_diff <- fission_dist$B_during_disp - fission_dist$AB_before_disp
fission_dist$A_move_diff <- fission_dist$A_during_disp - fission_dist$AB_before_disp
#get speed of group before and speeds of subgroups after
fission_dist$B_speed <- fission_dist$B_during_disp / fission_dist$d.duration
fission_dist$A_speed <- fission_dist$A_during_disp / fission_dist$d.duration
fission_dist$AB_speed <- fission_dist$AB_before_disp / fission_dist$b.duration

plotdir <- "C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/"
png(height = 400, width = 580, units = 'px', filename = paste0(plotdir,'speed_before_split.png'))
hist(fission_dist$AB_speed, breaks = 30, main = "", xlab = "Speed before split (m/s)", col = "aquamarine4")
dev.off()

#removing rows where the before speed is 0
fission_dist <- fission_dist[which(fission_dist$b.duration > 0),]

#calculate speed diff
fission_dist$speed_diff <- abs(fission_dist$B_speed - fission_dist$A_speed)
fission_dist$B_speed_diff <- fission_dist$B_speed - fission_dist$AB_speed
fission_dist$A_speed_diff <- fission_dist$A_speed - fission_dist$AB_speed

hist(fission_dist$speed_diff)
plot(fission_dist$AB_speed*60, fission_dist$speed_diff)
plot(fission_dist$B_during_disp, fission_dist$A_during_disp)

plot(fission_dist$AB_before_disp, fission_dist$move_diff)
points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff >20 & fission_dist$AB_before_disp > 25,], col = "steelblue2", pch = 19)
points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff >20 & fission_dist$AB_before_disp <= 25,], col = "indianred", pch = 19)
points(move_diff~AB_before_disp, fission_dist[fission_dist$move_diff <20 & fission_dist$AB_before_disp <= 25,], col = "gold2", pch = 19)


#looks like we can categorize based on some criteria:
#for fissions:
#if the group were moving before the split, then the stopping group and continue moving group - so could look at which subgroup changed their speed of travel most in response to the full group before

#make a column to get the subgroup that changed distance more - either by stopping when the group were moving or from moving when the group were relatively stationary
fission_dist$change_dist_subgroup <- NA
fission_dist$change_dist_subgroup[which(abs(fission_dist$B_move_diff) > abs(fission_dist$A_move_diff))] <- "b_dist_changed"
fission_dist$change_dist_subgroup[which(abs(fission_dist$B_move_diff) < abs(fission_dist$A_move_diff))] <- "a_dist_changed"
fission_dist$change_dist_subgroup[which(abs(fission_dist$move_diff) < 10)] <- "ab_dist_same"

#get the subgroup that changed their speed from the full groups speed
fission_dist$change_speed_subgroup <- NA
fission_dist$change_speed_subgroup[which(abs(fission_dist$B_speed_diff) > abs(fission_dist$A_speed_diff))] <- "b_speed_changed"
fission_dist$change_speed_subgroup[which(abs(fission_dist$B_speed_diff) < abs(fission_dist$A_speed_diff))] <- "a_speed_changed"
fission_dist$change_speed_subgroup[which(abs(fission_dist$speed_diff) < 0.01)] <- "ab_speed_same"

#want to bind the info on which subgroups changed speed to the call rates df
speed_filt <- fission_dist[,c("event_idx", "change_speed_subgroup", "A_speed_diff", "B_speed_diff")]

ind_fission_data <- merge(ind_events_data, speed_filt, by = "event_idx")

ind_fission_data_long <- ind_fission_data %>%
  pivot_longer(c(agg_call_rate, contact_call_rate), names_to = "call", values_to = "rate")

#change factor levels so before is shown before after in plot
ind_fission_data_long$period <- factor(ind_fission_data_long$period, levels = c("before","during" ,"after"))


#get the points to correspond to the period
ggplot(data = ind_fission_data_long[ind_fission_data_long$event_type == "fission",],
            aes(x = call, y = rate, fill = period))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_dodge(width = 0.75), size = 0.5, color = "gray3", aes(group = interaction(call, period))) +
  ylim(c(0,0.75))+
  facet_wrap(~subgroup+change_speed_subgroup)

ind_fission_data_long_bef <- ind_fission_data_long[ind_fission_data_long$period == "before",]
#get the A group and B group speeds in the same column
ind_fission_data_long_bef <- ind_fission_data_long_bef %>%
  pivot_longer(c("A_speed_diff", "B_speed_diff"), names_to = "group", values_to = "speed_diff")
#just looking at contact calls for plotting
ind_fission_data_long_bef <- ind_fission_data_long_bef[ind_fission_data_long_bef$call == "contact_call_rate",]

png(height = 400, width = 580, units = 'px', filename = paste0(plotdir,'callrate_speed.png'))
plot(ind_fission_data_long_bef$speed_diff, ind_fission_data_long_bef$rate, xlab = "speed difference", ylab = "call rate", pch = 16, cex = 0.5)
dev.off()

#next thing to do is remove the NA subgroups as this doesn't help us, and perhaps look at filtering for fission when the distance the subgroups move from one another significantly changes (as the group was likely not fully together before the fission event so their calling behaviour may not show any change)

#removing rows where the subgroup ID is NA
ind_fission_data_long_filt <- ind_fission_data_long[!(ind_fission_data_long$subgroup== "NA"),]
ind_fission_data_long_filt$period <- factor(ind_fission_data_long_filt$period, levels = c("before","during" ,"after"))

#want to combine the group A with b_speed_changed and group B with a_speed_changed 
#want to combine the group A with a_speed_changed and group B with b_speed_changed
#this is so there's just two plots where we have the calling rate of the changing group and the non-changing group

ind_fission_data_long_filt$change_group <- NA

ind_fission_data_long_filt <- within(ind_fission_data_long_filt,{
  change_group = NA
  change_group[subgroup == "A" & change_speed_subgroup == "a_speed_changed"] = "change"
  change_group[subgroup == "B" & change_speed_subgroup == "b_speed_changed"] = "change"
  change_group[subgroup == "A"& change_speed_subgroup == "b_speed_changed"] = "not_change"
  change_group[subgroup == "B"& change_speed_subgroup == "a_speed_changed"] = "not_change"
})

#for now: remove cases where the ab_speed stayed the same
ind_fission_data_long_filt <- ind_fission_data_long_filt[!(ind_fission_data_long_filt$change_speed_subgroup== "ab_speed_same"),]

g <- ggplot(data = ind_fission_data_long_filt[ind_fission_data_long_filt$event_type == "fission",],
       aes(x = call, y = rate, fill = period))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_dodge(width = 0.75), size = 0.5, color = "gray3", aes(group = interaction(call, period))) +
  ylim(c(0,0.75))+
  theme_classic()+
  facet_wrap(~change_group)

g
#here we can see an increase in contact call rate in the changing speed group 
ggsave(paste0(plotdir, "call_change.png"), width = 10, height = 5)



#---------FUSIONS------------------------------------------------------------------------------------------------------------
#need to get the distance travelled of the full group after the fusion and update the scripts here for this to work...

fusion_dist <- dist_travel_df[dist_travel_df$event_type == "fusion",]
fusion_dist <- fusion_dist[!is.na(fusion_dist$B_during_disp),]

#get distance travelled
fusion_dist$move_diff <- abs(fusion_dist$B_during_disp - fusion_dist$A_during_disp)
fusion_dist$B_move_diff <- fusion_dist$B_during_disp - fusion_dist$AB_before_disp
fusion_dist$A_move_diff <- fusion_dist$A_during_disp - fusion_dist$AB_before_disp
#get speed of group before and speeds of subgroups after
fusion_dist$B_speed <- fusion_dist$B_during_disp / fusion_dist$d.duration
fusion_dist$A_speed <- fusion_dist$A_during_disp / fusion_dist$d.duration
fusion_dist$AB_speed <- fusion_dist$AB_before_disp / fusion_dist$b.duration

plotdir <- "C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/"
png(height = 400, width = 580, units = 'px', filename = paste0(plotdir,'speed_before_split.png'))
hist(fusion_dist$AB_speed, breaks = 30, main = "", xlab = "Speed before split (m/s)", col = "aquamarine4")
dev.off()

#removing rows where the before speed is 0
fusion_dist <- fusion_dist[which(fusion_dist$b.duration > 0),]

#calculate speed diff
fusion_dist$speed_diff <- abs(fusion_dist$B_speed - fusion_dist$A_speed)
fusion_dist$B_speed_diff <- fusion_dist$B_speed - fusion_dist$AB_speed
fusion_dist$A_speed_diff <- fusion_dist$A_speed - fusion_dist$AB_speed

hist(fusion_dist$speed_diff)
plot(fusion_dist$AB_speed*60, fusion_dist$speed_diff)
plot(fusion_dist$B_during_disp, fusion_dist$A_during_disp)

plot(fusion_dist$AB_before_disp, fusion_dist$move_diff)
points(move_diff~AB_before_disp, fusion_dist[fusion_dist$move_diff >20 & fusion_dist$AB_before_disp > 25,], col = "steelblue2", pch = 19)
points(move_diff~AB_before_disp, fusion_dist[fusion_dist$move_diff >20 & fusion_dist$AB_before_disp <= 25,], col = "indianred", pch = 19)
points(move_diff~AB_before_disp, fusion_dist[fusion_dist$move_diff <20 & fusion_dist$AB_before_disp <= 25,], col = "gold2", pch = 19)





