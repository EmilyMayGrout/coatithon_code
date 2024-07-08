#Coati time synch
#Script to synchronize sorokas on coatis to UTC time
#This script produces an offset value for each date / soroka ID combo, based on a linear fit to the synch calls
#It then uses these offset values (stored in the table offset_table) to synchronize all the call labels to UTC
#The final synched calls are stored in a csv with the same name as the original input file plus "_synched", in the same folder

#The set of synch calls comes from an input file stored at path_to_synch_file, where the time in file and time in UTC of each synch call are noted.
#The set of calls comes from an input file stored at path_to_call_labels_file

#The first synch call in the morning and afternoon are used, if present. Others are not used, for now
#Things that are NOT taken into account in this script are:
# -- speed of sound / distance to sound source
# -- within-file drift (this is low, approx 200 ms/hr), probably negligible compared to the uncertainty in UTC times from the talking clock
# -- the original synch UTC times are only accurate within ~ 1 sec
#Because of this, the adjusted times will probably be accurate to within a few sec, but this should be measured somehow

#----PARAMETERS---
utc_offset <- -5 #local time offset from UTC, in hours - should be negative 5 for panama
#path_to_synch_file <- '~/Downloads/Synch calls coati collar project - Sheet1.csv' #where the synch table is located
path_to_synch_file <- 'C:/Users/egrout/Dropbox/coatithon/synch_table_all.csv'

use_machine_labels <- T

if(use_machine_labels){
path_to_call_labels_file <- 'C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/all_data_hms_all_ml.csv' #made in calls_around_fission_fusions_galaxy_old with labels which were cleaned in the cleaning_labels 
} else {
  path_to_call_labels_file <- 'C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/all_data_hms.csv'
}

#-----SETUP----
#LIBRARIES
library(lubridate)
library(stringr)
library(tools)
library(hms)

#FUNCTIONS
#---Convert audition label to seconds into the file
#x: the string from audition
#sec: number of seconds into the file
convert_audition_to_sec <- function(x){
  
  #if there is a space between sec and millisec, turn it into a period (some labels contained spaces rather than periods)
  x <- gsub(' ', '.', x)
  
  #get number of colons 
  n_colons <- str_count(x, ':')
  if(length(n_colons)==0 | is.na(n_colons)){
    return(NA)
  }
  
  #if there are too many or too few colons (should be 1 or 2), return NA
  if(n_colons != 1 & n_colons != 2){
    return(NA)
  }
  
  #if 1 colon, assume it is in the format m:s
  if(n_colons==1){
    sec <- ms(x)
  }else{
    
    #if 1 colons, assume it is in the format h:m:s
    sec <- as_hms(x)
  }
  
  return(as.numeric(sec))
}

#------MAIN-----
#load synch calls
synch_dat <- read.csv(path_to_synch_file, stringsAsFactors = F)

#add column for filename_short
pattern <- "G\\d+_\\d+_FL\\d+" #need to do this incase the collar ID is 5 numbers not 4
synch_dat$filename_short <- regmatches(synch_dat$filename, regexpr(pattern, synch_dat$filename))

#convert date to Y-M-D
synch_dat$date <- as.Date(sub("(..)$", "20\\1", synch_dat$date), "%d.%m.%Y")

#remove unknowns (for now) - may want to return to this later
synch_dat$X.1..synch.time.in.file[grepl("unk", synch_dat$X.1..synch.time.in.file)] <- NA
synch_dat$X.4..synch.time.in.file[grepl("unk", synch_dat$X.4..synch.time.in.file)] <- NA

#get soroka number (1 or 2)
synch_dat$soroka_num <- sapply(synch_dat$filename, FUN = function(x){strsplit(x, '_')[[1]][2]})

#---REAL TIMES IN UTC

#Convert start time in file to UTC
synch_dat$filestart_UTC <- as.POSIXct(paste(synch_dat$date, 
                                            synch_dat$start.time.in.file..st.), 
                                      tz = 'UTC',
                                      format = '%Y-%m-%d %H:%M:%S') - utc_offset*60*60

#Convert end time in file to UTC
synch_dat$fileend_UTC <- as.POSIXct(paste(synch_dat$date, 
                                            synch_dat$end.time.in.file..et.), 
                                      tz = 'UTC',
                                      format = '%Y-%m-%d %H:%M:%S') - utc_offset*60*60

#Convert first synch time to UTC
synch_dat$synch1_UTC <- as.POSIXct(paste(synch_dat$date, 
                                          synch_dat$X.1..real.time), 
                                    tz = 'UTC',
                                    format = '%Y-%m-%d %H:%M:%S') - utc_offset*60*60

#Convert second synch time to UTC
synch_dat$synch2_UTC <- as.POSIXct(paste(synch_dat$date, 
                                          synch_dat$X.4..real.time), 
                                    tz = 'UTC',
                                    format = '%Y-%m-%d %H:%M:%S') - utc_offset*60*60

#Convert file start and end times (time in file) to seconds since beginning of file
synch_dat$filestart_sec <- sapply(synch_dat$time.in.file.of.st, convert_audition_to_sec)
synch_dat$fileend_sec <- sapply(as_hms(synch_dat$time.in.file.of.et), convert_audition_to_sec)

#Convert synch times in file to seconds into the file
synch_dat$synch1_sec <- sapply(synch_dat$X.1..synch.time.in.file, convert_audition_to_sec)
synch_dat$synch2_sec <- sapply(synch_dat$X.4..synch.time.in.file, convert_audition_to_sec)

#get the time elapsed between start of file and synch time according to the time in file and the UTC difference
dt_start_to_synch1_sec <- synch_dat$synch1_sec - synch_dat$filestart_sec
dt_start_to_synch1_UTC <- as.numeric(difftime(synch_dat$synch1_UTC, synch_dat$filestart_UTC, tz = 'UTC', units = 'secs'))
offsets_start_to_synch1 <- dt_start_to_synch1_sec - dt_start_to_synch1_UTC
synch_dat$offset_1 <- offsets_start_to_synch1

#get the time elapsed between start of file and second synch time according to the time in file and the UTC difference
dt_start_to_synch2_sec <- synch_dat$synch2_sec - synch_dat$filestart_sec
dt_start_to_synch2_UTC <- as.numeric(difftime(synch_dat$synch2_UTC, synch_dat$filestart_UTC, tz = 'UTC', units = 'secs'))
offsets_start_to_synch2 <- dt_start_to_synch2_sec - dt_start_to_synch2_UTC
synch_dat$offset_2 <- offsets_start_to_synch2

#drift within file - we are going to assume this is negligible which seems reasonable based on the numbers all being within 0-2 sec
#dt_synch1_to_synch2_sec <- synch_dat$synch2_sec - synch_dat$synch1_sec
#dt_synch1_to_synch2_UTC <- as.numeric(difftime(synch_dat$synch2_UTC, synch_dat$synch1_UTC, tz = 'UTC', units = 'secs'))
#offsets_synch1_to_synch2 <- dt_synch1_to_synch2_sec - dt_synch1_to_synch2_UTC

#plot offset vs actual UTC time within a collar and fit linear models to offsets
synch_dat$soroka_id <- paste(synch_dat$id, synch_dat$soroka_num, sep='_')
soroka_ids <- unique(synch_dat$soroka_id)
plot(NULL, xlim = c(min(synch_dat$synch1_UTC,na.rm=T), max(synch_dat$synch1_UTC,na.rm=T)), ylim = c(-120,0))
colors <- rainbow(length(soroka_ids))
for(i in 1:length(soroka_ids)){
  soroka_id <- soroka_ids[i]
  idxs <- which(synch_dat$soroka_id == soroka_id)
  times_UTC <- c(synch_dat$synch1_UTC[idxs], synch_dat$synch2_UTC[idxs])
  offsets <- c(synch_dat$offset_1[idxs], synch_dat$offset_2[idxs])
  points(times_UTC, col = colors[i],pch =i, offsets)
}

#run fits for collars with at least 2 points and get coeffs
intercepts <- slopes <- rep(NA, length(soroka_ids))
names(intercepts) <- names(slopes) <- soroka_ids
for(i in 1:length(soroka_ids)){
  soroka_id <- soroka_ids[i]
  idxs <- which(synch_dat$soroka_id == soroka_id)
  
  #collect up data of offset vs utc time
  times_UTC <- c(synch_dat$synch1_UTC[idxs], synch_dat$synch2_UTC[idxs])
  offsets <- c(synch_dat$offset_1[idxs], synch_dat$offset_2[idxs])
  
  min_time <- min(times_UTC, na.rm=T)
  max_time <- max(times_UTC, na.rm=T)
  time_elapsed <- as.numeric(difftime(max_time, min_time, tz = 'UTC', units = 'days'))
  
  
  if(time_elapsed >= 1){
    fit <- lm(offsets ~ times_UTC)
    intercepts[i] <- fit$coefficients[1]
    slopes[i] <- fit$coefficients[2]
  }
}

#make a table of soroka id and date vs offset and use to correct labels
#pull the dates and soroka ids from the table of labelled calls in the future, but for now use the synch table
dates <- seq.Date(min(synch_dat$date), max(synch_dat$date), by = "1 day")
soroka_ids <- unique(synch_dat$soroka_id)

#create the table
offset_table <- data.frame(soroka_id = rep(soroka_ids, each = length(dates)), date = rep(dates, length(soroka_ids)))
offset_table$offset <- NA

#for each row, find soroka id, get model fit for that soroka id, then interpolate to get predicted offset
for(i in 1:nrow(offset_table)){
  
  #get date and convert to seconds since 1970-01-01
  date <- offset_table$date[i]
  date_time <- as.POSIXct(paste(date, '12:30:00'), tz = 'UTC') #middle of the file
  date_time_numeric <- as.numeric(date_time)
  
  #get soroka id
  soroka_id <- offset_table$soroka_id[i]
  
  #get the correct slope and intercept
  slope <- slopes[which(names(slopes)==soroka_id)]
  intercept <- intercepts[which(names(intercepts)==soroka_id)]
  
  #get the offset
  offset <- intercept + slope*date_time_numeric
  
  #store the offset
  offset_table$offset[i] <- offset
  
}

#if there are more than one offset (based on the two synch calls)

#----PART 2: Use synch table to synchronize call labels
calls <- read.csv(path_to_call_labels_file, sep=',')

#get soroka id of calls in calls table
calls$soroka_id <- sapply(calls$file, function(x){return(paste(strsplit(x, '_')[[1]][1:2],collapse='_'))})
calls$soroka_id_date <- paste(calls$soroka_id, calls$date, sep = '_')

#get file start and end of each corresponding file according to the soroka clock
#use the marker in the file which gives the UTC time according to soroka clock (filestart_UTC) at a given (early) time in the file (filestart_sec)
calls$filestart_UTC_soroka <- synch_dat$filestart_UTC[match(calls$file, synch_dat$filename_short)] - synch_dat$filestart_sec[match(calls$file, synch_dat$filename_short)]

#match to offset table to get the predicted offset for each file
offset_table$soroka_id <- paste0('G', offset_table$soroka_id)
offset_table$soroka_id_date <- paste(offset_table$soroka_id, offset_table$date, sep='_')
calls$offset <- offset_table$offset[match(calls$soroka_id_date, offset_table$soroka_id_date)]

#subtract the offset to get the synched time
calls$datetime_synch <- as.POSIXct(calls$datetime, tz = 'UTC', format = "%Y-%m-%d %H:%M:%OS") - calls$offset
#the milliseconds are not displayed but they are there: format(calls$datetime_synch, "%Y-%m-%d %H:%M:%OS3" )


#save synched calls table
base <- file_path_sans_ext(path_to_call_labels_file)
outfile_name <- paste(base, 'synched.csv', sep = '_')
write.csv(calls, file = outfile_name, quote = F)

#TODO: 
#get more synchs for soroka ids that are missing them
#check and fix the wonky ones
#check whether the offset table is consistent with the original synch table
#CHECK whether the final synching makes sense (maybe with another synch sound, and also looking at aggressive sequences to see if they line up across tags?)
#eventually: label synchs in all files - label earliest and latest synch you can find