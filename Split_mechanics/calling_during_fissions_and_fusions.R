#This script is for looking at how call rates relate to fission-fusion dynamics in coatis

#directory holding all the data
datadir <- '~/Dropbox/coatithon/calling_during_fissions_and_fusions/data'

callfile <- 'all_data_hms_synched.csv'
ff_file <- 'galaxy_auto_ff_events_characterized.RData'
gps_file <- 'galaxy_xy_highres_level1.RData'
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

#create a new data frame with data at the individaul level
ind_events_data <- data.frame()
for(i in 1:nrow(group_events_data)){
  event_idx <- group_events_data$event_idx[i]
  ind_idxs <- 1:nrow(coati_ids)
  rows <- data.frame(event_idx = rep(event_idx, length(ind_idxs)),
                     ind_idx = ind_idxs)
  rows_all <- rbind(rows, rows, rows)
  rows_all$period <- c(rep('before',length(ind_idxs)), 
                       rep('during',length(ind_idxs)),
                       rep('after',length(ind_idxs)))
  ind_events_data <- rbind(ind_events_data, rows_all)
}

#compute the call rate during an event period for a given individual
ind_idx <- 1
event_idx <- 1



