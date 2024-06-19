#FOR GALAXY GROUP

#This script:
#	Removes unrealistic speeds and replaces with NAs ( > 99.95 percentile = 14.8 m/s)
#	Fills in missing data gaps less than max.interp = 5 sec with linear interpolation
#	Finds instances where coati did not move during an NA gap of < max.move time = 5 minutes (moved < max.move=5 m) and replaces 
#		them with the mean location of the coati between start and end of seq
# Find extreme distances (0.01% quantile or > 99.99% quantile of xs or ys for that individual), if they are not next to other values within 1000 m of that point, replace with NA
# Find extreme xs and ys above mean + sd * 10 for each ind and remove those

#output of this code will make a new gps_file called "presedente_xy_highres_level2.RData"
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/'
gps_file <- "galaxy_xy_highres_level1.RData" 
id_file <- 'galaxy_coati_ids.RData'

outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#---------------PARAMETERS-------------
max.interp <- 5 #maximum number of seconds of NAs to interpolate (default = 5)
max.speed.quantile <- 0.9995 #maximum speed quantile to allow - speeds above this will be replaced with NAs (default = 0.9995, corresponds to 14.8 m/s)
max.move <- 5 #maximum distance moved (in meters) during an NA gap of duration <  max.move.time seconds, to replace all with the mean of the start and end point (default = 5 m)
max.move.time <- 5*60 #maximum NA gap time (in seconds) to replace with mean value if the individual did not move more than max.move meters (default = 5*60 sec = 5 min)
max.dist.quantile <- .9999 #maximum distance quantile. distances beyond this quantile in both x and y direction for each individual will be replaced with NAs (defaults to 0.9999)
min.dist.quantile <- .0001 #minimum distance quantile. distances below this quantile in both x and y direction for each individual will be replaced with NAs (defaults to 0.0001)

#LIBRARIES
library(lubridate)

#number of indidivudals
n.inds <- nrow(xs)
n.times <- ncol(xs)


#Find extreme speeds and replace with NAs
speeds <- matrix(NA,nrow=n.inds,ncol=n.times)

for(i in 1:n.inds){
  
  speeds[i, seq(1,n.times-1)] <- sqrt(diff(xs[i,])^2 + diff(ys[i,])^2)

}

#finding the maximum speed they are likely to travel
max.speed <- quantile(speeds, max.speed.quantile, na.rm=T)

#need to remove the high speed calculated when the day changes - these values could already have an NA?
#first get the difference in time between ts to find the change in days
diff_time <- vector("numeric", length(ts) - 1)

# Loop over each timestamp (except the last one)
for (i in 1:(length(ts) - 1)){
  
  # Calculate the time difference between consecutive timestamps
  ts_diff <- difftime(ts[i+1], ts[i])
  
  # Store the time difference in the vector
  diff_time[i] <- ts_diff
  
}

unique(diff_time)
#it worked!
#now get the indexes of day changes

speeds[, (which(diff_time > 1))]
#looks like speed isn't calculated when the time changes to the next day...

#to double check this, finding the times when the speed is high
#looking at the first individual

#for (i in 1:n.inds){
#print(ts[which(speeds[i,] > 15)])
#}
#from inspecting these timestamps, the high speed isn't calculated between the change in days so this code shouldn't remove incorrect gps 



#replace unrealistic speeds with NAs
for(i in 1:n.inds){
  unrealistic.speed.idxs <- which(speeds[i,] > max.speed)
  xs[i,unrealistic.speed.idxs] <- NA
  ys[i,unrealistic.speed.idxs] <- NA
  xs[i,unrealistic.speed.idxs + 1] <- NA
  ys[i,unrealistic.speed.idxs + 1] <- NA
}

#Interpolate through seqs of NAs of length < max.interp
for(i in 1:n.inds){
  
  #data for a single individual
  x <- xs[i,]
  y <- ys[i,]
  
  #vectors to hold interpolated data
  x.interp <- x
  y.interp <- y
  
  #get runs of NAs
  runs <- rle(is.na(x)) #get runs of NAs
  vals <- runs$values
  lens <- runs$lengths
  idxs <- c(0, cumsum(runs$lengths))
  
  #for each run of NAs, fill in with linearly interp values if length is less than max.interp
  for(j in 1:length(vals)){
    if(vals[j]==TRUE){
      first.idx <- idxs[j]+1
      last.idx <- idxs[j+1]
      
      #If not too near the beginning or the end of the sequence...
      if((first.idx > 1) & (last.idx < length(x))){
        
        #Get values before and after NA sequence
        prev.val.x <- x[first.idx-1]
        next.val.x <- x[last.idx + 1]
        prev.val.y <- y[first.idx-1]
        next.val.y <- y[last.idx +1]
        
        #Fill in with linear interpolation if the NA sequence is short (< 5)
        if(lens[j] <= max.interp){
          interp.vals.x <- seq(prev.val.x,next.val.x,length.out=lens[j]+2)
          interp.vals.y <- seq(prev.val.y,next.val.y,length.out=lens[j]+2)
          x.interp[seq(first.idx-1,last.idx+1)] <- interp.vals.x
          y.interp[seq(first.idx-1,last.idx+1)] <- interp.vals.y
        }
        
        #Otherwise...if less than 5 minutes gap...
        if((lens[j] > max.interp) & (lens[j] < max.move.time)){
          #Fill in with mean value at start and end if they are close enough ( <= max.move)
          dist.moved <- sqrt((next.val.x - prev.val.x)^2 + (next.val.y - prev.val.y)^2)
          time.elapsed <- last.idx - first.idx
          if(dist.moved < max.move){
            mean.x <- mean(c(next.val.x,prev.val.x))
            mean.y <- mean(c(next.val.y,prev.val.y))
            x.interp[seq(first.idx,last.idx)] <- mean.x
            y.interp[seq(first.idx,last.idx)] <- mean.y
          }
        }
      }
    }
  }
  
  xs[i,] <- x.interp
  ys[i,] <- y.interp
  
}

#Find and remove points that are most likely GPS error
#We do this in two phases, 
#     1. first looking for some extreme x or y values (determined by max/min.dist.quantile; more extreme than 0.9999 or 0.0001 quantile) 
#        and checking if the previous and next point are in a similar area (w/in 1km). If not, this extreme point is replaced with NA
#     2. second we find SUPER extreme x or y values (>mean + 10*sd for each x or y) and remove them no matter what. 

##--------- STEP 1 of removing extreme  points --
for(i in 1:n.inds){
  print('ind:')
  print(i)
  xi <- xi.new <- xs[i,]
  yi <- yi.new <- ys[i,]
  non.nas <- which(!is.na(xi))
  
  #get very large or very small values of x and y (outside their normal range)
  bigs <- unique(c(which(xi>quantile(xi,max.dist.quantile,na.rm=T)),which(yi>quantile(yi,max.dist.quantile,na.rm=T))))
  smalls <- unique(c(which(xi<quantile(xi,min.dist.quantile,na.rm=T)),which(yi<quantile(yi,min.dist.quantile,na.rm=T))))
  extremes <- unique(c(bigs,smalls))
  
  
  #for each extreme value, find the previous and next non-NA data point
  #if this previous point is more than 1 km away, just replace the unrealistic location with NA
  #otherwise, leave it
  for(j in 1:length(extremes)){
    t.idx <- extremes[j]
    prev.t <- max(non.nas[which(non.nas < t.idx)])
    next.t <- min(non.nas[which(non.nas > t.idx)])
    
    dist.prev <- sqrt((xs[i,prev.t]-xs[i,t.idx])^2 + (ys[i,prev.t]-ys[i,t.idx])^2)
    dist.next <- sqrt((xs[i,next.t]-xs[i,t.idx])^2 + (ys[i,next.t]-ys[i,t.idx])^2)
    
    if(dist.prev > 1000 & dist.next > 1000){
      xi.new[t.idx] <- NA
      yi.new[t.idx] <- NA
    }
    
  }
  
  xs[i,] <- xi.new
  ys[i,] <- yi.new
  
}

##--------- STEP 2 of removing extreme  points --
#Find and remove super unrealistic locations (>mean + 10*sd for each ind in x or y)
for(i in 1:n.inds){
  xi <- xs[i,]
  yi <- ys[i,]
  extremes.high.x <- which(xi > mean(xi,na.rm=T)+10*sd(xi,na.rm=T))
  extremes.low.x <- which(xi < mean(xi,na.rm=T)-10*sd(xi,na.rm=T))
  extremes.high.y <- which(yi > mean(yi,na.rm=T)+10*sd(yi,na.rm=T))
  extremes.low.y <- which(yi < mean(yi,na.rm=T)-10*sd(yi,na.rm=T))
  extremes <- c(extremes.high.x,extremes.high.y,extremes.low.x,extremes.low.y)
  
  if(length(extremes)>0){
    xs[i,extremes] <- NA
    ys[i,extremes] <- NA
  }
  
}

#still has the GPS error after the cleaning
#this only goes when we change the max dist quantile to 0.9350 (0.77m/s) and this removes real data so might not be the best option to filter the errors
#could there be a way of filtering the really fast looping movement - which seems to be the way the GPS errors go
plot(xs[, 26000:28000], ys[, 26000:28000])


#Save as RData

save(file=paste0(data_dir, 'galaxy_xy_highres_level2.RData'),list=c('xs','ys','ts'))




