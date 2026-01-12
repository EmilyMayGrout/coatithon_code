#Preprocessing the coati GPS data from Presidente group whom were collared again in Dec 2024 to March 2025 on BCI

library(lubridate)
library(sf)
library(tidyverse)

#load useful functions
source('C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/coati_function_library_V1.R')

firsttime <- as.POSIXct('2024-12-16 11:00', tz = 'UTC')
lasttime <-  as.POSIXct('2025-02-15 23:00', tz = 'UTC')

indir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/gps/"
outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2025/presidente/"
metadatadir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/metadata/"

#setting the working directory to the raw data file in dropbox
setwd(indir)

#------------------------------------------------------------------------------------------
#BELOW CODE NOW COMPLETE SO DOES NOT NEED RERUNNING TO GET TXT FILES

#NEED TO DO THIS IF DOWNLOADING DATA FROM MOVEBANK RATHER THAN USING DECODER
#first read in csv and split into files for each individual and save as txt file into the indir
# presidente_all <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2025/presidente/movebank_nutritional_landscapes.csv")
# 
# #filter to just coatis, as this is also spider monkeys
# presidente_all <- presidente_all[presidente_all$individual.taxon.canonical.name == "Nasua narica",]

#filter presidente_all to "name", "lon", "lat", "date", "time")
# 
# presidente_all_filter <- presidente_all[,c(36,4,5,3)]
# 
# colnames(presidente_all_filter) <- c("name", "lon", "lat", "datetime")
# 
# #split by id
# split_presidente <- split(presidente_all_filter, presidente_all_filter$name, drop = F)
# 
# #save each individual as own txt file in dropbox rawdata file
# 
# allNames <- names(split_presidente)
# for(i in allNames){
#   #get name of individual from i of list
#   saveName <- paste0(i ,".txt")
#   write.table(split_presidente[[i]], file = saveName)
# 
# }


#once an hour from 23-11, every 5 mins on 1C collars from 11-23, and every 10 mins on 1A collars from 11-23

#create times vector - only for the times that are in the DAY 11:00 (6am) - 23:00 (6pm)
#first make list of 10 min intervals then remove the times overnight
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '5 min')
minhr <- 11
maxhr <- 23
#minhr <- 00 #to get sleep times for josephine
#maxhr <- 10

ts <- ts[which((hour(ts) <= maxhr) & (hour(ts) >= minhr))]

#making a list of all individuals in presidente
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
  
  #take the modulus 10 (divide by 10 and get remainder) to find values which are within 2 min of a 10 min interval
  tagdata <- tagdata[which(minute(tagdata$datetime) %% 5 <= 2),]
  
  #round down the times to the nearest 10 mins
  tagdata$roundtime <- floor_date(tagdata$datetime, unit = "5 minutes")
  
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

#looking for the first and last dates the collars were on for:
# NonNAindex <- which(!is.na(xs[1,]))
# firstNonNA <- min(NonNAindex)
# lastNonNA <- max(NonNAindex)
# 
# ts[firstNonNA]
# ts[lastNonNA]

setwd(metadatadir)
coati_ids <- read.csv("coati_ids.csv", header = F)
colnames(coati_ids) <- c("name", "tag_id", "age", "sex")
coati_ids$color <- '#0000FF'
coati_ids$color[which(coati_ids$age == 'Adult' & coati_ids$sex == 'Female')] <- '#FF0000'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Female')] <- '#FFAA66'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Male')] <- '#66AAFF'
coati_ids$color[which(coati_ids$age == 'Juvenile')] <- '#666666'
save(coati_ids, file = paste0(outdir, 'presidente2025_coati_ids.RData'))

save(list=c('xs','ys','ts'), file = paste0(outdir,'presidente2025_xy_10min_level0.RData'))
save(list=c('lats','lons','ts'), file = paste0(outdir,'presidente2025_latlon_10min_level0.RData'))  



#processing the ripeness csv data sent to me on 03.09.2025

ripe <- read.csv("C:/Post-doc/nutriscapes/dipx_ripeness.csv")
ripe$treeID <- sub("_.*", "",ripe$ID)
ripe$fruitID <- sub("^[^_]+_[^_]+_(\\d+).*", "\\1", ripe$ID)

ripe <- ripe[,c(1,3,6,7)]
#remove invalid rows
ripe <- ripe %>%
  filter(!is.na(dmy(Date))) 

ripe$Ripeness[ripe$Ripeness == "ripper"] <- "riper"
ripe$Ripeness[ripe$Ripeness == "ripe/riper"] <- "riper"
ripe$Ripeness[ripe$Ripeness == "ripe/old"] <- "ripe"
ripe$Ripeness[ripe$Ripeness == "riper/old"] <- "ripe"

ripe <- ripe[!(ripe$Ripeness == "" | ripe$Ripeness == "eaten" | ripe$Ripeness == "empty" | ripe$Ripeness == "Lacmellea panamensis") ,]



#remove the replicates
ripe <- ripe[!duplicated(ripe), ]

#make treeID column numeric and remove the non-dipteryx fruit
ripe$treeID <- as.numeric(ripe$treeID)
#remove rows with NA in treeID (as were the other species names)
ripe <- ripe[!is.na(ripe$treeID),]

ripe$Date <- as.Date(ripe$Date, format = "%d-%m-%Y")

save(ripe, file = paste0(outdir,'dipx_ripeness.RData'))






















