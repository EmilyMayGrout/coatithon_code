#this script will read in the raw txt files form e-obs collars and
#output a preprocess r data file for the low res data (1 gps point/10 mins)
library(lubridate)
library(sf)
library(tidyverse)


#load functions
source('C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/coati_function_library_V1.R')

firsttime <- as.POSIXct('2021-12-24 11:00', tz = 'UTC')
lasttime <-  as.POSIXct('2022-01-13 23:00', tz = 'UTC')

indir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/gps/binded_files/"
outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
metadatadir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/metadata/"

#setting the working directory to the raw data file in dropbox
setwd(indir)

#the file structure for the rbind code below has been changed so the commented code no longer works, but the binded data is in the binded_files folder
#rbind the last day to each gps txt file
#9460
#all_files[2]
#all_files[3]
#tag_9460 <- read.table(all_files[2])
#tag_9460_2 <- read.table(all_files[3])
#tag_9460_all <- rbind(tag_9460, tag_9460_2)
#write.table(tag_9460_all,"binded_files/tag9460_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9463
#all_files[4]
#all_files[5]
#tag_9463 <- read.table(all_files[4])
#tag_9463_2 <- read.table(all_files[5])
#tag_9463_all <- rbind(tag_9463, tag_9463_2)
#write.table(tag_9463_all,"binded_files/tag9463_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9464
#all_files[6]
#all_files[7]
#tag_9464 <- read.table(all_files[6])
#tag_9464_2 <- read.table(all_files[7])
#tag_9464_all <- rbind(tag_9464, tag_9464_2)
#write.table(tag_9464_all,"binded_files/tag9464_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9466
#all_files[8]
#all_files[9]
#tag_9466 <- read.table(all_files[8])
#tag_9466_2 <- read.table(all_files[9])
#tag_9466_all <- rbind(tag_9466, tag_9466_2)
#write.table(tag_9466_all,"binded_files/tag9466_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9467
#all_files[10]
#all_files[11]
#tag_9467 <- read.table(all_files[10])
#tag_9467_2 <- read.table(all_files[11])
#tag_9467_all <- rbind(tag_9467, tag_9467_2)
#write.table(tag_9467_all,"binded_files/tag9467_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9470
#all_files[12]
#all_files[13]
#tag_9470 <- read.table(all_files[12])
#tag_9470_2 <- read.table(all_files[13])
#tag_9470_all <- rbind(tag_9470, tag_9470_2)
#write.table(tag_9470_all,"binded_files/tag9470_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9471
#all_files[14]
#all_files[15]
#tag_9471 <- read.table(all_files[14])
#tag_9471_2 <- read.table(all_files[15])
#tag_9471_all <- rbind(tag_9471, tag_9471_2)
#write.table(tag_9471_all,"binded_files/tag9471_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9474
#all_files[16]
#all_files[17]
#tag_9474 <- read.table(all_files[16])
#tag_9474_2 <- read.table(all_files[17])
#tag_9474_all <- rbind(tag_9474, tag_9474_2)
#write.table(tag_9474_all,"binded_files/tag9474_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9475
#all_files[18]
#all_files[19]
#tag_9475 <- read.table(all_files[18])
#tag_9475_2 <- read.table(all_files[19])
#tag_9475_all <- rbind(tag_9475, tag_9475_2)
#write.table(tag_9475_all,"binded_files/tag9475_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9476
#all_files[20]
#all_files[21]
#tag_9476 <- read.table(all_files[20])
#tag_9476_2 <- read.table(all_files[21])
#tag_9476_all <- rbind(tag_9476, tag_9476_2)
#write.table(tag_9476_all,"binded_files/tag9476_gps.txt", row.names = F, col.names = F,  quote = FALSE)

#9480
#all_files[22]
#all_files[23]
#tag_9480 <- read.table(all_files[22])
#tag_9480_2 <- read.table(all_files[23])
#tag_9480_all <- rbind(tag_9480, tag_9480_2)
#write.table(tag_9480_all,"binded_files/tag9480_gps.txt", row.names = F, col.names = F,  quote = FALSE)



#----------------------------------------------------------------------------------

#create times vector - only for the times that are in the DAY 11:00 (6am) - 23:00 (6pm)
#first make list of 10 min intervals then remove the times overnight
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '10 min')
minhr <- 11
maxhr <- 23
ts <- ts[which((hour(ts) <= maxhr) & (hour(ts) >= minhr))]

#making a list of all individuals in galaxy
all_files <- sort(list.files())

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))

for(i in 1:length(all_files)){

  #filter columns to ones needed
  tagdata <- read.table(all_files[i], sep =  ",")
  #getting the correct dat and time column
  tagdata <- tagdata[, c( 6, 7, 14, 16)]
  colnames(tagdata) <- c("lon", "lat", "date", "time")
  #make datetime format
  tagdata$datetime <- paste(tagdata$date, tagdata$time)
  tagdata$datetime <- as.POSIXct(tagdata$datetime, format="%d.%m.%Y %H:%M:%S", tz="UTC")
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
  
  
  #round down the times to the nearest 10 mins
  tagdata$roundtime <- floor_date(tagdata$datetime,unit="10 minutes")
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
coati_ids <- read.csv("coati_id.csv", header = F)
colnames(coati_ids) <- c("name", "tag_id", "age", "sex")
coati_ids$color <- '#0000FF'
coati_ids$color[which(coati_ids$age == 'Adult' & coati_ids$sex == 'Female')] <- '#FF0000'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Female')] <- '#FFAA66'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Male')] <- '#66AAFF'
coati_ids$color[which(coati_ids$age == 'Juvenile')] <- '#666666'
save(coati_ids, file = paste0(outdir, 'galaxy_coati_ids.RData'))

save(list=c('xs','ys','ts'), file = paste0(outdir,'galaxy_xy_10min_level0.RData'))
save(list=c('lats','lons','ts'), file = paste0(outdir,'galaxy_latlon_10min_level0.RData'))  

#calculate the proportion of missing fixes
each_sum <- data.frame(sum = rowSums(!is.na(xs)))
max(each_sum$sum)
min(each_sum$sum)
each_sum$missing <- max(each_sum$sum)- each_sum$sum
each_sum$prop <- (each_sum$missing/max(each_sum$sum))*100
100 - mean(each_sum$prop) #96.46707
100 - max(each_sum$prop) #93.02885
sd(each_sum$prop)#2.280898



