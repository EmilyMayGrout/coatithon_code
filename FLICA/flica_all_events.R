#GALAXY GROUP - this script is for reading in all the times of fission and fusion events, then doing FLICA analysis on each event

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
#sigma = 0.75
#timeWindow = 240 #usually 1/10th of duration
obj1 <- mFLICA(TS=event_dat[focal_tags,sample_vals,],timeWindow=60,timeShift=1,sigma=0.75)


#------------------------------------------------------------------------------------------
#PUTTING ALL FLICA RESULTS INTO A LIST for each event
max_time <- 300

all_events =list()
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
 both_groups <-  c(groupA, groupB)

 name <- paste("event_", i, "_indx", gal_events_detected$tidx[i], "_", gal_events_detected$event_type[i], sep="")
 all_events[[name]] <- mFLICA(TS=event_dat[both_groups,sample_vals,],timeWindow=60,timeShift=10,sigma=0.75)
  
}

#saved all_events with max_time = 300 and timeShift = 10
#save(all_events, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_list.Rdata")
#load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/all_events_list.Rdata")


#plotting one event:
plotMultipleTimeSeries(TS=all_events$event_148_indx172154_fusion$dyNetOut$dyNetBinDensityVec, strTitle="Network Density")

plotMultipleTimeSeries(TS=all_events$event_148_indx172154_fusion$factionSizeRatioTimeSeries, strTitle="Faction Size Ratios",TSnames = coati_ids$name[both_groups]) + scale_color_brewer(palette="Paired")

#plotting all inds except individuals not in the group
dat_filt <- event_dat[both_groups,,]

plotMultipleTimeSeries(TS=event_dat[both_groups,,1],strTitle="x axis", TSnames = coati_ids$name[both_groups])
plotMultipleTimeSeries(TS=event_dat[both_groups,,2],strTitle="y axis", TSnames = coati_ids$name[both_groups])



#finding the individuals that had most influence for each event
d <- all_events$event_148_indx172154_fusion$factionSizeRatioTimeSeries
#plotting one individuals influence
plot(d[5, 1:300])

df <- data.frame(matrix(nrow = length(d[1,]), ncol = 2))

#for loop for getting the individual with the greatest influence for that event and putting it in dataframe over time
for (i in 1:ncol(d)){
  col <- d[,i]
  id <- which.max(col)
  df[i,1] <-  id
  df[i,2] <- coati_ids$name[id]
}

tabl <- table(df$X2)
pie(tabl)






























