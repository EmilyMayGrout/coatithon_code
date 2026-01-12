#This script is looking at the coatis Rich collared on BCI to see whether they did fission and fusion
library(lubridate)
library(sf)
library(tidyverse)


indir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2024/gps/"
outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2024/rich_project/"
metadatadir <- "C:/Users/egrout/Dropbox/coatithon/rawdata/2024/metadata/"

setwd(indir)

#downloaded file from movebank
#first read in csv and split into files for each individual and save as txt file into the indir
bci_all <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2024/rich_project_bci.csv")
colnames(bci_all)

unique(bci_all$individual.taxon.canonical.name)

coati_all <- bci_all[bci_all$individual.taxon.canonical.name == "Nasua narica",]

#filter bci_all to "name", "lon", "lat", "date", "time")

coati_all_filt <- coati_all[,c(36,4,5,3)]

colnames(coati_all_filt) <- c("name", "lon", "lat", "datetime")

#split by id
split_bci <- split(coati_all_filt, coati_all_filt$name, drop = F)

#can make each individual as own dataframe in global enviro
#list2env(lapply(split_bci, as.data.frame.list), .GlobalEnv)

#save each individual as own txt file in dropbox rawdata file
allNames <- names(split_bci)
for(i in allNames){
  #get name of individual from i of list
  saveName <- paste0(i ,".txt")
  write.table(split_bci[[i]], file = saveName)

}

test <- split_bci$Ginny

for (i in 1:nrow(test)){
  
  test$time_diff[i] <- difftime(test$datetime[i+1], test$datetime[i])

}

hist(as.numeric(test$time_diff), breaks = 100)

#so the sampling rate is a GPS burst every 4 mins and if they aren't moving, then its 30 min sampling 
#Time of First Deployed Location	2024-03-09 16:46:27.000
#Time of Last Deployed Location	2024-07-11 20:35:28.000


firsttime <- as.POSIXct('2024-05-04 11:00', tz = 'UTC')
lasttime <-  as.POSIXct('2024-05-30 23:00', tz = 'UTC')


#create times vector - only for the times that are in the DAY 11:00 (6am) - 23:00 (6pm)
#first make list of 10 min intervals then remove the times overnight
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '30 min')
minhr <- 11
maxhr <- 23
#minhr <- 00 #to get sleep times for josephine
#maxhr <- 10

ts <- ts[which((hour(ts) <= maxhr) & (hour(ts) >= minhr))]

#making a list of all individuals from Rich's project on BCI
all_files <- sort(list.files())

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))


for(i in 1:length(all_files)){
  
  #filter columns to ones needed
  tagdata <- read.table(all_files[i])
  #getting the correct dat and time column
  tagdata <- tagdata[, c(2, 3, 4)]
  colnames(tagdata) <- c("lon", "lat", "datetime")
  #make datetime format
  tagdata$datetime <- as.POSIXct(tagdata$datetime, format="%Y-%m-%d %H:%M:%OS", tz="UTC")
  #getting the last gps point per burst
  tagdata$test <- !duplicated(tagdata$datetime, fromLast = T )
  
  #remove rows with 0's from the lat and lon
  tagdata <- tagdata[which(tagdata$lon != 0),]
  
  #high res data
  highresdata <- tagdata[which(hour(tagdata$datetime) >= 11 & hour(tagdata$datetime) < 14),]
  
  #getting the difference between each time to get the last gps point from each burst
  diffs <- difftime(tagdata$datetime[2:nrow(tagdata)], tagdata$datetime[1:(nrow(tagdata)-1)])
  diffs <- c(diffs, 600)
  tagdata$diffs <- diffs
  highresdata$diffs <- NA #need to add these columns for the rbind later
  
  #filter the gps points to the last fix
  tagdata <- tagdata[which(tagdata$diffs > 160) ,]
  
  #take the modulus 10 (divide by 10 and get remainder) to find values which are within 2 min of a 10 min interval
  tagdata <- tagdata[which(minute(tagdata$datetime) %% 10 <= 2),]

  #round down the times to the nearest 4 mins
  tagdata$roundtime <- floor_date(tagdata$datetime,unit="30 minutes")
  highresdata$roundtime <- highresdata$datetime #need this for the matching later
  
  #add in the high res data
  tagdata <- rbind(tagdata, highresdata)
  
  #match times to get lons and lats at each time for that individual
  lon <- tagdata$lon[match(ts, tagdata$roundtime)]
  lat <- tagdata$lat[match(ts, tagdata$roundtime)]
  
  lats[i,] <- lat
  lons[i,] <- lon
  
  #to convert to UTM, need to use sf function and for this we need to combine the matrices rows to a dataframe
  combined_df <- data.frame(lat = lats[i,], lon = lons[i,])
  
  #convert NA's to zero for sf functions to work
  combined_df$lon[is.na(combined_df$lon)] <- 0
  combined_df$lat[is.na(combined_df$lat)] <- 0
  
  #convert to UTM
  #first need to give the latlon data the correct CRS so it converts to UTM correctly
  latlon_i <- st_as_sf(x=combined_df, coords=c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
  #convert to UTM
  utm_i <- st_transform(latlon_i, crs="+proj=utm +zone=17 +north +datum=WGS84 +units=m") #convert UTM to lat/long - coordinates are stored in a geometry column of class 'sfc_POINT'
  
  #store eastings and northings in xs and ys matrices
  xs[i,] <- unlist(map(utm_i$geometry,1))
  ys[i,] <- unlist(map(utm_i$geometry,2))
  
}


setwd(metadatadir)
coati_ids <- read.csv("coati_ids.csv", header = F)
colnames(coati_ids) <- c("name", "tag_id", "age", "sex")

coati_ids <- read.csv("coati_ids.csv", header = F)
colnames(coati_ids) <- c("name", "tag_id", "age", "sex")
coati_ids$color <- '#0000FF'
coati_ids$color[which(coati_ids$age == 'Adult' & coati_ids$sex == 'Female')] <- '#FF0000'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Female')] <- '#FFAA66'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Male')] <- '#66AAFF'
coati_ids$color[which(coati_ids$age == 'Juvenile')] <- '#666666'
save(coati_ids, file = paste0(outdir, 'richcoati_ids.RData'))


save(list=c('xs','ys','ts'), file = paste0(outdir,'richproject_xy_30min_level0.RData'))
save(list=c('lats','lons','ts'), file = paste0(outdir,'richproject_latlon_30min_level0.RData'))  



