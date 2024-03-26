#This script is for looking at how call rates relate to fission-fusion dynamics in coatis

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
calls$datetime_synch_pos<-as.POSIXct(calls$datetime_synch, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") + 5*60*60
labeled_periods$tstart<-as.POSIXct(labeled_periods$starttime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") + 5*60*60
labeled_periods$tstop<-as.POSIXct(labeled_periods$stoptime, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC") + 5*60*60

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

