#For Galaxy group, cleaning the events to investigate ff dynamics - as the NA's in the times are edge cases which were not true events. Also removing duplicate events caused by an individual joining the group together


#read in functions
codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"

#navigate into directory
setwd(groupdir)

#read in coati ids
load('galaxy_coati_ids.RData')
#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})
#read in timestamp data
load('galaxy_xy_highres_level2.RData')

setwd(codedir)
source('coati_function_library.R')

load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_auto_ff_events_characterized.RData")

events_clean <- events

#removing duplicate events - events: 29 (2021-12-28 12:00:04) and 68 (2022-01-04 12:05:20), 80 (2022-01-05 12:17:52) and 81 (2022-01-05 12:18:07) are a dup of the previous and after event 
events_clean <- events_clean[!(events_clean$event_idx == 29),]
events_clean <- events_clean[!(events_clean$event_idx == 68),]
events_clean <- events_clean[!(events_clean$event_idx == 80),]
events_clean <- events_clean[!(events_clean$event_idx == 81),]

#t <- events_clean[,c(1,2,3,8,9,10)]
#write.csv(t, paste0(groupdir,"events_cleaned_level2.csv"))

#going through these events on google sheets to check - they all look good!

#event 93 on 6th 11:20:04 is wrong time - centroid doesn't change as subgroup is spread around Gus's location, so will replace the time with NAs

events_clean[events_clean$event_idx == 93, 11:21] <- NA
events_clean[events_clean$event_idx == 99, 11:21] <- NA #the time is wrong for this event

#remove events where there are NA's in the times
events_clean <- events_clean[!is.na(events_clean$after_time),]

for(i in c(1:nrow(events_clean))){
  print(i)
  analyse_ff_event(i, events_clean, xs, ys, ts, plot=T, max_time = 700)
}






