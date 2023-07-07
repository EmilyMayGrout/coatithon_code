#GALAXY GROUP - this script is for reading in all the times of fission and fusion events, then doing FLICA analysis on each event

##ISSUE that the mFLICA reverts id's back to 1,2,3 once individuals are filtered to the event inds. How is this fixed?

#read in libraries
library(mFLICA)
library(R.matlab)
library(raster)
library(plotly)

#ID file
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_coati_ids.RData")
#GPS file
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_xy_highres_level1.RData")
#gal_events_detected (made in identify_splits_merges.R automated events)
load('C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/gal_events_detected.Rda')

#making data into an array
dat <- array(NA, dim = c(nrow(xs), ncol(xs), 2))
dat[,,1] <- xs
dat[,,2] <- ys

spd <- sqrt((apply(xs, 1, diff))^2 + (apply(ys, 1, diff))^2)
spd <- t(spd)
spd_a <- array(NA, dim = c(nrow(xs), 1))
spd <- cbind(spd, spd_a) #adding the NA's to the last column so they match

#need to remove events that are less that 10 minutes after 11am for each day
gal_events_detected$time <- format(gal_events_detected$datetime, format="%H:%M:%S")
gal_events_detected <- subset(gal_events_detected, time > "11:10:00") 

#time in seconds before and after event 
max_time <- 600

#looking at 1 event:
i = 5

#get first time indx and last time indx to extract the fission and fusion time
ti <- gal_events_detected$tidx[i] - max_time
tf <- gal_events_detected$tidx[i] + max_time

eventtime <- ts[ti:tf]
event <- dat[,ti:tf,]

event_dat <- event

sample_vals = 1:length(eventtime)

#replace NAs with 0 
event_dat[is.na(event_dat)] <- 0

## function from raster - not working
## sampling individuals that are within dist_thres
dist_thres <- 100 # removing individuals that do not travel with the group based on median over the all period
dst <- pointDistance(cbind(c(apply(event_dat[,,1],c(1),median)), 
                           c(apply(event_dat[,,2],c(1),median))),
                     c(median(event_dat[,,1], na.rm = TRUE), 
                       median(event_dat[,,2], na.rm = TRUE)),lonlat = FALSE)
focal_tags <- dst < dist_thres

# run flica code with tags that are in the group
#better to have a smaller time shift for more accurate results - 10 is good, but the smaller the number the longer it takes to run
#sigma = 0.75 #lower sigma is lower number of factions
#timeWindow = 240 #usually 1/10th of duration
obj1 <- mFLICA(TS=event_dat[focal_tags,sample_vals,],timeWindow=60,timeShift=10,sigma=0.75)



#------------------------------------------------------------------------------------------
#PUTTING ALL FLICA RESULTS INTO A LIST for each event
max_time <- 300

all_events_speed = list()
all_events_loc = list()
for (i in 1:nrow(gal_events_detected)){
  #get first time indx and last time indx to extract the fission and fusion time
  ti <- gal_events_detected$tidx[i] - max_time
  tf <- gal_events_detected$tidx[i] + max_time
  
  eventtime <- ts[ti:tf]
  event <- dat[,ti:tf,]
  event_dat <- event
  sample_vals = 1:length(eventtime)
  
  #replace NAs with 0 
  event_dat[is.na(event_dat)] <- 0
  
  #because the group were in 2 subgroups the dst was too high (despite the distance within the group being smaller than 100m), so need to find new way of extracting the indx of the individuals involved in the event. 
 groupA <-  gal_events_detected$group_A_idxs[i][[1]]
 groupB <-  gal_events_detected$group_B_idxs[i][[1]]
 both_groups <-  sort(c(groupA, groupB))
 

 name <- paste("event_", i, "_indx", gal_events_detected$tidx[i], "_", gal_events_detected$event_type[i], sep="")
 #for location
 all_events_loc[[name]] <- mFLICA(TS=event_dat[both_groups,sample_vals,],timeWindow=60,timeShift=10,sigma=0.75)
 
 #for speed
 all_events_speed[[name]] <- mFLICA(TS=spd[both_groups,sample_vals ],timeWindow=60,timeShift=10,sigma=0.75)
 
 
 for (i in 1:length(all_events_loc[[name]]$leadersTimeSeries)){
   
   all_events_loc[[name]]$leadersTimeSeries[i] = list(both_groups[unlist(all_events_loc[[name]]$leadersTimeSeries[i])])
   
   all_events_loc[[name]]$factionMembersTimeSeries[i] = list(both_groups[unlist(all_events_loc[[name]]$factionMembersTimeSeries[i])])
   
   all_events_speed[[name]]$leadersTimeSeries[i] = list(both_groups[unlist(all_events_speed[[name]]$leadersTimeSeries[i])])
   
   all_events_speed[[name]]$factionMembersTimeSeries[i] = list(both_groups[unlist(all_events_speed[[name]]$factionMembersTimeSeries[i])])
   
 }
  
}

#saved all_events with max_time = 300 and timeShift = 10 #sigma is 0.75 for events_loc and 0.45 for events_speed
#save(all_events_loc, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_loc.Rdata")
#save(all_events_speed, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_speed.Rdata")
#load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_loc.Rdata")
#load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_speed.Rdata")


#plotting one event:
plotMultipleTimeSeries(TS=all_events_loc$event_1_indx2132_fission$dyNetOut$dyNetBinDensityVec, strTitle="Location Network Density")
plotMultipleTimeSeries(TS=all_events_speed$event_1_indx2132_fission$dyNetOut$dyNetBinDensityVec, strTitle="Speed Network Density")

#to look at subgroup (faction) membership - first number is the leader of that subgroup
all_events_loc$event_1_indx2132_fission$factionMembersTimeSeries[36]
all_events_speed$event_1_indx2132_fission$factionMembersTimeSeries[36]

plotMultipleTimeSeries(TS=all_events_loc$event_1_indx2132_fission$factionSizeRatioTimeSeries, strTitle="Faction Size Ratios",TSnames = coati_ids$name[both_groups]) + scale_color_brewer(palette="Paired")

plotMultipleTimeSeries(TS=all_events_speed$event_3_indx2211_fission$factionSizeRatioTimeSeries, strTitle="Faction Size Ratios",TSnames = coati_ids$name[both_groups]) + scale_color_brewer(palette="Paired")

#plotting all inds except individuals not in the group
dat_filt <- event_dat[both_groups,,]

plotMultipleTimeSeries(TS=event_dat[both_groups,,1],strTitle="x axis", TSnames = coati_ids$name[both_groups])
plotMultipleTimeSeries(TS=event_dat[both_groups,,2],strTitle="y axis", TSnames = coati_ids$name[both_groups])



#finding the individuals that had most influence for each event
#faction size ratio is a number of edges that connect between faction-member nodes divided by a number of total nodes within a following network. If a leader has a higher faction-size ratio, then it has more followers than a leader with a lower faction-size ratio. A faction-size ratio has a value between 0 and 1.
l <- all_events_loc$event_1_indx2132_fission$factionSizeRatioTimeSeries[,]
s <- all_events_speed$event_1_indx2132_fission$factionSizeRatioTimeSeries[,]

#get the individual
i = 3
both_groups <- sort(c(gal_events_detected$group_A_idxs[i][[1]], gal_events_detected$group_B_idxs[i][[1]]))


#plotting one individuals influence
plot(l[1, 1:300])

df <- data.frame(matrix(nrow = length(l[1,]), ncol = 2))
#for loop for getting the individual with the greatest influence for that event and putting it in dataframe over time

for (i in 1:ncol(l)){
  col <- l[,i]
  id <- which.max(col)
  df[i,1] <-  id
  df[i,2] <- coati_ids$name[both_groups[id]]
}

tabl <- table(df$X2)
pie(tabl)

#------------------------------------------------------------------------------
#this for loop is getting the overall leader of the event plus the leader before and after the event
df_all_events <- data.frame(event_id = 1:length(all_events_loc), leader_indx = NA, leader_id = NA)
df_single <- data.frame(matrix(nrow = length(l[1,]), ncol = 2))

#for loop through each event to get the individual who lead before the event and after the event - save into df_all_events
#can look at speed influence by changing all_events_loc to all_events_speed

#something is wrong as Quasar is in the groups where she wasn't involved in that event!!
#issue is that when the mFLICA is run on the individuals interested, the index of those individuals is lost (goes back to 1,2,3 ect...) Not sure how to fix this??

for (j in 1:length(all_events_loc)){

  l <- all_events_loc[j][[1]] #assign each event to "l"
  l <- l$factionSizeRatioTimeSeries[,]
  
  for (i in 1:ncol(l)){
    col <- l[,i]
    id <- which.max(col)
    df_single[i,1] <-  id
    df_single[i,2] <- coati_ids$name[both_groups[id]]
  }
  
  before_event <- df_single[1:max_time,]
  after_event <- df_single[(max_time+1):(max_time*2+1),]
  #get the coati_indx
  df_all_events$leader_indx[j] <- tail(names(sort(table(df_single$X1))), 1)
  #get the coati_name
  df_all_events$leader_id[j] <- tail(names(sort(table(df_single$X2))), 1)
  #get the coati_indx before the event
  df_all_events$before_event[j] <- tail(names(sort(table(before_event$X1))), 1)
  #get coati_id before event
  df_all_events$before_event_id[j] <- tail(names(sort(table(before_event$X2))), 1)
  #get the coati_indx after the event
  df_all_events$after_event[j] <- tail(names(sort(table(after_event$X1))), 1)
  #get coati_id after event
  df_all_events$after_event_id[j] <- tail(names(sort(table(after_event$X2))), 1)
}


hist(as.numeric(df_all_events$leader_indx))
hist(as.numeric(df_all_events$before_event))
hist(as.numeric(df_all_events$after_event))

df_all_events$group_A <- gal_events_detected$group_A
df_all_events$group_B <- gal_events_detected$group_B













