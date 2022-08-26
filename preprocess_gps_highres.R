#this script will read in the raw txt files form e-obs collars and
#output a preprocess r data file for the HIGH res data (1 gps point/sec)
library(lubridate)

#load useful functions
source('C:/Users/egrout/Dropbox/coatithon/coatithon_code/coati_function_library.R')

firsttime <- as.POSIXct('2021-12-24 11:00', tz = 'UTC')
lasttime <-  as.POSIXct('2022-01-09 23:00', tz = 'UTC')
#to do: get the last day of data and rbind it to the same file - then change the last date

indir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/gps/binded_files"
outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
metadatadir <-  "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/metadata/"

#setting the working directory to the raw data file in dropbox
setwd(indir)

#create times vector - only for the times that are in the DAY 11:00 (6am) - 23:00 (6pm)
#first make list of 10 min intervals then remove the times overnight
ts <- seq.POSIXt(from = firsttime, to = lasttime,  by = '1 sec')
minhr <- 11
maxhr <- 14
ts <- ts[which((hour(ts) <= maxhr) & (hour(ts) >= minhr))]

#making a list of all individuals in galaxy
all_files <- sort(list.files())

lats <- lons <- xs <- ys <- matrix(NA, nrow = length(all_files), ncol = length(ts))

i = 1

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
  
  #to do: check whether the 160s time between fixes is right - to account for time taken to get GPS fix
  #filter the gps points to the last fix
  tagdata <- tagdata[which(tagdata$diffs > 160) ,]
  
  #take the modulus 10 (divide by 10 and get remainder) to find values which are within 2 min of a 10 min interval
  tagdata <- tagdata[which(minute(tagdata$datetime) %% 10 <= 2),]
  
  
  #round down the times to the nearest 10 mins
  tagdata$roundtime <- floor_date(tagdata$datetime,unit="10 minutes")
  highresdata$roundtime <- highresdata$datetime #need this for the matching later
  
  #add in the high res data
  tagdata <- highresdata
  
  
  #match times to get lons and lats at each time for that individual
  lon <- tagdata$lon[match(ts, tagdata$roundtime)]
  lat <- tagdata$lat[match(ts, tagdata$roundtime)]
  
  lats[i,] <- lat
  lons[i,] <- lon
  
  #convert to UTM
  eastsNorths <- latlon.to.utm(LonsLats = cbind(lon, lat), utm.zone = 17, southern_hemisphere = F)
  
  #store eastings and northings in xs and ys matrices
  xs[i,] <- eastsNorths[,1]
  ys[i,] <- eastsNorths[,2]
  
}

setwd(metadatadir)
coati_ids <- read.csv("coati_id.csv", header = F)
colnames(coati_ids) <- c("name", "tag_id", "age", "sex")
coati_ids$color <- '#0000FF'
coati_ids$color[which(coati_ids$age == 'Adult' & coati_ids$sex == 'Female')] <- '#FF0000'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Female')] <- '#FFAA66'
coati_ids$color[which(coati_ids$age == 'Sub-adult' & coati_ids$sex == 'Male')] <- '#66AAFF'
coati_ids$color[which(coati_ids$age == 'Juvenile')] <- '#666666'
save(coati_ids, file = paste0(outdir, 'coati_ids.RData'))

save(list=c('xs','ys','ts'), file = paste0(outdir,'galaxy_xy_highres_level0.RData'))
save(list=c('lats','lons','ts'), file = paste0(outdir,'galaxy_latlon_highres_level0.RData'))  






