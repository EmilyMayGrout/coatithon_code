#For Presidente group, cleaning the events to investigate ff dynamics - as the NA's in the times are edge cases which were not true events. Also removing duplicate events caused by an individual joining the group together


#read in functions
codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"

#navigate into directory
setwd(groupdir)

#read in coati ids
load('presedente_coati_ids.RData')
#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})
#read in timestamp data
load('presedente_xy_highres_level2.RData')

setwd(codedir)
source('coati_function_library.R')

load("C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/presedente_auto_ff_events_characterized.RData")

events_clean <- events

#remove events where there are NA's in the times
events_clean <- events_clean[!is.na(events_clean$after_time),]

#although removing the NAs reduced the number of errors from the start of recording, it didn't remove all of them
#need to remove the errors from start of recording 
events_clean <- events_clean[!(events_clean$event_idx == 25),]
events_clean <- events_clean[!(events_clean$event_idx == 217),]


for(i in c(nrow(events_clean):220)){
  print(i)
  analyse_ff_event(i, events_clean, xs, ys, ts, plot=T, max_time = 700)
}


t <- events_clean[,c(1,2,3,8,9,10)]
write.csv(t, paste0(groupdir,"events_cleaned_level2.csv"))

#Where are the events on the 2nd of Feb??

















