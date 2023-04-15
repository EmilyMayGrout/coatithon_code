#library of general functions

#LIBRARIES
library(dbscan)
library(rgdal)
library(lubridate)
library(stringr)


#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	EastingsCol1: whether eastings should be given in first column of output (default) or not
#Outputs:
#	EastNorths or NorthEasts: [N x 2 matrix] of Eastings and Northings - eastings are first column by default
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='34',southern_hemisphere=TRUE){
  latlons <- data.frame(X=LonsLats[,2],Y=LonsLats[,1])
  non.na.idxs <- which(!is.na(latlons$X) & !is.na(latlons$Y))
  len <- nrow(latlons)
  non.na.latlons <- latlons[non.na.idxs,]
  coordinates(non.na.latlons) <- ~Y + X
  proj4string(non.na.latlons) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  utm <- spTransform(non.na.latlons,CRS(projection.string))
  EastNorths <- matrix(NA,nrow=len,ncol=2)
  EastNorths[non.na.idxs,] <- utm@coords
  if(!EastingsCol1){
    NorthEasts <- EastNorths[,c(2,1)]
    return(NorthEasts)
  } else{
    return(EastNorths)
  }
}

#Converts a matrix of eastings and northings (eastings first column, northings second column) to latlong
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	LonsCol1: whether lons should be given in first column of output (default) or not
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column by default 
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '34',southern_hemisphere=TRUE){
  utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
  non.na.idxs <- which(!is.na(utms$X) & !is.na(utms$Y))
  len <- nrow(utms)
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  non.na.utms <- SpatialPoints(utms[non.na.idxs,],proj4string=CRS(projection.string))
  lonlat <- spTransform(non.na.utms,CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
  LonLats <- matrix(NA,nrow=len,ncol=2)
  LonLats[non.na.idxs,] <- lonlat@coords
  if(!LonsCol1){
    LatLons <- LonLats[,c(2,1)]
    return(LatLons)
  } else{
    return(LonLats)
  }	
}


#This function computes information about subgroup membership over time
#Inputs:
# R [numeric]: radius for DBSCAN
# xs, ys [n_inds x n_times matrix]: x and y locations of individuals 
#Outputs:
# subgroup_data: object that contains:
#   $ind_subgroup_membership [n_inds x n_times matrix]: matrix giving the subgroup id for each individual at each time
#   $n_subgroups [n_times vector]: vector of the number of subgroups over time
#   $subgroup_counts [max_n_subgroups x n_times matrix]: matrix giving the number of individuals in each subgroup (with NAs when there are fewer subgroups than the max)
#   $R [numeric]: radius used in DBSCAN
get_subgroup_data <- function(xs, ys, R){
  
  
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  sub_groups <- matrix(NA, nrow = n_inds, ncol = n_times)
  
  for (t_idx in 1:n_times){
    
    #getting the x and y coordinates for each time point
    x <- xs[ , t_idx]
    y <- ys[, t_idx]
    
    non_nas <- which(!is.na(x))
    
    x_non_na <- x[non_nas]
    y_non_na <- y[non_nas]
    
    n_tracked <- length(non_nas)
    
    if(n_tracked >= 1){
      #scan for each time column to get the cluster number for each individual
      groups <- dbscan(x= cbind(x_non_na, y_non_na), eps = R, minPts = 1)
      #saving the cluster for each time into the matrix
      sub_groups[non_nas, t_idx] <- groups$cluster
    }
    
  }
  
  #get max number of groups per time
  ngroups <- suppressWarnings(apply(sub_groups, 2, max, na.rm = T)) #suppress warnings because if all NAs it gives a warning, this is fine because we deal with it in the next step by removing the infinite maxes
  #remove infinite values
  ngroups[is.infinite(ngroups)] <- NA

  
  
  #now want to look at size of the sub-groups as the histogram shows 2 subgroups being most frequent 
  #but this might be because of one individual not being with the group
  
  group_max <- max(ngroups, na.rm = T)
  
  subgroup_counts <- matrix(NA, nrow = group_max, ncol = n_times)
  
  #for loop for getting the number of individuals in each sub group every 10 mins (when the radius for defining within group is 50m, there are maximum 5 groups in full dataset)
  for(i in 1:n_times){
    
    #going through each time stamp
    current_group <- sub_groups[, i]
    #number of individuals in each subgroup 
    counts <- table(current_group)
    #to remove NA's when there is no data from all individuals
    if(length(counts) >= 1){
      #putting the data into the matrix
      subgroup_counts[1:length(counts),i] <- sort(counts)
      
    }
  }
  
  #store all data in a list
  subgroup_data <- list()
  subgroup_data$ind_subgroup_membership <- sub_groups
  subgroup_data$n_subgroups <- ngroups
  subgroup_data$subgroup_counts <- subgroup_counts
  subgroup_data$R <- R
  
  return(subgroup_data)
  
}

#Make a network visualization with colors
#INPUT:
#net[adjacency matrix of the network]
#coati_ids[dataframe of the coati ids]
#OUTPUT:
#plot
visualize_network_matrix <- function(net, coati_ids){
  
  zmin <- min(net, na.rm=T)
  zmax <- max(net, na.rm=T)
  par(mgp=c(3, 1, 0), mar=c(11,10,6,4)) #bottom, left, top, and right ##CHANGE THESE VALUES IF DOING FORLOOP OF THE MATRIX PLOTS AND THEY DON'T FIT
  image.plot(net, col = viridis(256), zlim=c(zmin,zmax), xaxt= 'n', yaxt = 'n', legend.cex = 5, legend.width = 5,legend.mar = 6, axis.args=list(cex.axis=2))
  axis(1, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2, cex.axis=1.5)
  axis(2, at = seq(0,1,length.out= nrow(net)), labels = coati_ids$name, las = 2,  cex.axis=1.5)
  
  points(rep(-.04, nrow(net)),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19)
  points(seq(0,1,length.out=nrow(net)),rep(-.04,n_inds),col=coati_ids$color, xpd = T, pch = 19)
}

#changed the size of the labels with cex.axis = 1.5, default is 1
#also changed the size of the legend axis values with: axis.args=list(cex.axis=2), remove if want default
#legend width was 1.3, changed to 5 for the multiple matrix plots




#Function for getting distance between individuals over time
#INPUT:
#xs, ys, radius
#OUTPUT
#object that contains:
  #   $distance over time: distance between each individual at each time point
  #   $proximity_network: probability that i and j are within a distance of R
  #   $r_within [numeric]: radius used for proximity network

get_proximity_data <- function(xs, ys, r_within){

  #get n_inds and n_times
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #get distance between individuals at these subset times
  
  #make array to store data into
  dist_over_time <- array(NA, dim = c(n_inds, n_inds, n_times))
  
  #for loop through each individual with each individual for each time point (where all individuals are together and have gps point)
  for(t in 1:ncol(xs)){
    for(i in 1:n_inds){
      for(j in 1:n_inds){
        
        #get x and y coordinates to calculate distance using pythagorus
        xcoord_i <- xs[i,t]
        xcoord_j <- xs[j,t]
        dx <- (xcoord_i - xcoord_j)
        ycoord_i <- ys[i,t]
        ycoord_j <- ys[j,t]
        dy <- (ycoord_i - ycoord_j)
        dist <- sqrt((dx)^2 + (dy)^2)
        
        #put the statements into the array at the correct position in the correct dimension
        dist_over_time[i, j, t] <- dist
      }
    }
  }
  
  #threshold the distances to determine if individuals are in proximity at each time
  net_over_time <- dist_over_time < r_within
  
  #construct proximity network
  proximity_net <- apply(net_over_time, MARGIN = c(1,2), FUN = mean, na.rm=T)
  diag(proximity_net) <- NA
  
  #store output
  proximity_data <- list()
  proximity_data$proximity_net <- proximity_net
  proximity_data$dist_over_time <- dist_over_time
  proximity_data$r_within <- r_within
  
  #return output
  return(proximity_data)

}


#Randomise the splits and construct new randomized data frame
#Keep original group the same, keep subgroup sizes the same, but randomize who goes to which subgroup
#INPUTS:
# splits_df: data frame with information about splits from real data
#OUTPUTS:
# splits_df_random: randomised version of splits data frame
randomise_splits <- function(splits_df){
  
  #start with the splits_df data frame from the data
  splits_df_rand <- splits_df
  
  #loop over each row and randomize who goes where
  for (i in 1:nrow(splits_df_rand)){
    
    #get the original group and sizes of subgroups for one row
    row <- splits_df[i,]
    orig_group <- row$orig_group[[1]]
    n_sub1 <- row$n_sub1
    n_sub2 <- row$n_sub2
    n_sub3 <- row$n_sub3
    
    #shuffle the original group to a random order
    orig_group_shuffled <- sample(orig_group)
    
    #allocate them to groups
    sub1 <- orig_group_shuffled[1:n_sub1]
    sub2 <- orig_group_shuffled[(length(sub1)+1):(length(sub1)+n_sub2)]
    
    splits_df_rand$sub1[i] <- list(sub1)
    splits_df_rand$sub2[i] <- list(sub2)
    
    if(n_sub3 > 0){
      sub3 <- orig_group_shuffled[(length(sub1)+length(sub2)+1):(n_sub1+n_sub2+n_sub3)]
      splits_df_rand$sub3[i] <- list(sub3)
    }
    
  }
  
  return(splits_df_rand)
  
}

#Helper function for get_consistency
#Returns q if q <= 0.5 and 1-q if q > 0.5
dist_to_0_or_1 <- function(q){
  
  if(is.na(q)){
    return(NA)
  }
  if(q > 0.5) {return(1-q)
   }else{return(q)}
  
}

#Compute our metric of consistency, C
#Let q_dyad = p_dyad_together
#m_dyad = 1-q_dyad if q_dyad > 0.5 or = q_dyad if q_dyad < 0.5 (so it's the distance to 0 or 1, whichever is closer)
#C = 1 - 2*mean(m_dyad)
#At the end, C is a measure of consistency where C = 0 if all q_dyad are at 0.5 and C = 1 if all q_dyad are either 0 or 1
#Why do we multiply by 2 and subtract from 1? The multiplying by 2 puts the number between 0 and 1 instead of between 0 and 0.5. The 1 minus makes it be consistency instead of inconsistency
#INPUTS:
# p_dyad_together: [matrix] of probabilities of splitting together given you were both in the original group
#OUTPUTS:
# C: [numeric] metric of consistency
get_consistency <- function(p_dyad_together){
  
  m_dyad <- sapply(p_dyad_together, FUN = dist_to_0_or_1)
  C <- 1 - 2*mean(m_dyad, na.rm=T)
  return(C)
}

#function that takes in a splits dataframe and outputs a p_dyad_together matrix (probability of splitting together for each dyad)
#INPUTS:
# splits_df_local: [data frame] containing information on group splits
# n_inds_local: [numeric] number of individuals in the full group
#OUTPUTS:
# p_dyad_together_local: [matrix n_inds_local x n_inds_local]: probability of being together in a split, given you were both in the original group 
get_p_dyad_together <- function(splits_df_local, n_inds_local){
  
  #COMPUTE METRIC OF P(STAY TOGETHER | originally together) for each dyad
  p_dyad_together_local <- array(NA, dim = c(n_inds_local,n_inds_local))
  #loop over dyads
  for(i in 1:(n_inds_local-1)){
    for(j in (i+1):n_inds_local){
      
      #find rows where they were both in the original group
      originally_together_rows <- which(unlist(lapply(splits_df_local$orig_group, FUN = function(x){return(i %in% x & j %in% x)})))
      
      #how many times do they end up in the same group
      still_together <- 0
      for(r in originally_together_rows){
        if(i %in% splits_df_local$sub1[r][[1]] & j%in% splits_df_local$sub1[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub2[r][[1]] & j %in% splits_df_local$sub2[r][[1]]){
          still_together <- still_together + 1
        }
        if(i %in% splits_df_local$sub3[r][[1]] & j %in% splits_df_local$sub3[r][[1]]){
          still_together <- still_together + 1
        }
        
      }
      
      p_dyad_together_local[i,j] <- still_together / length(originally_together_rows)
      
    }
  }
  return(p_dyad_together_local)
  
}

#Function to match names in fission/fusion manual events table to coati ids idxs
#INPUTS:
# subgroup_names: a character string containing the names of the coatis in a given subgroup
# coati_ids: table of coati ids to match to
#OUTPUTS:
# subgroup_idxs: a vector of the indexes associated with the subgroup members
match_coati_names <- function(subgroup_names, coati_ids){
  names_split <- strsplit(subgroup_names,',')[[1]]
  names_split <- sub(' ','',names_split)
  names_short <- sapply(names_split, function(x){return(substr(x,1,3))})
  
  #match names to get individual indexes
  subgroup_idxs <- match(names_short, coati_ids$name_short)
  return(subgroup_idxs)
  
}

#analyse and make a visualization of the fission-fusion event
#INPUTS:
# events: a data frame of manually-labeled ff events
# i: index to the row to use
# xs, ys: matrices of x and y coordinates
#OUTPUTS:
# some info (TBD) about the ff event
# a plot showing (top) dyadic distance over time and (bottom) a visualization of trajectories
analyse_ff_event <- function(i, events, xs, ys, max_time = 1200, thresh_h = 50, thresh_l = 15, plot = T){
  t_event <- events$tidx[i] #time of the event
  group_A <- events$group_A_idxs[i][[1]] #group A individual idxs
  group_B <- events$group_B_idxs[i][[1]] #group B individual idxs
  group_A_names <- events$group_A[i]
  group_B_names <- events$group_B[i]
  ti <- t_event - max_time #initial time to plot
  tf <- t_event + max_time #final time to plot
  event_type <- events$event_type[i]
  datetime <- events$datetime[i]
  nA <- length(group_A)
  nB <- length(group_B)
  nT <- ncol(xs)
  
  #get x and y coordinates of the relevant individuals in each subgroup
  xA <- matrix(xs[group_A,],nrow=nA,ncol=nT)
  xB <- matrix(xs[group_B,],nrow=nB,ncol=nT)
  yA <- matrix(ys[group_A,],nrow=nA,ncol=nT)
  yB <- matrix(ys[group_B,],nrow=nB,ncol=nT)
  
  #get centroids of the two subgroups
  xcA <- colMeans(xA, na.rm=T)
  ycA <- colMeans(yA, na.rm=T)
  xcB <- colMeans(xB, na.rm=T)
  ycB <- colMeans(yB, na.rm=T)
  
  #get distance between centroids
  dyad_dist <- sqrt((xcA - xcB)^2 + (ycA - ycB)^2)
  
  #classify the dyadic distance into categories:
  #0 = below lower thershold
  #1 = between thresholds
  #2 = above higher threshold
  dyad_dist_event <- dyad_dist[ti:tf]
  after_idxs <- (max_time+1):(2*max_time+1) #indexes after the marked event
  middle_idxs <- (max_time / 2):(max_time*3/2)
  before_idxs <- 1:max_time #indexes before the marked event
  if(event_type == 'fission'){
    if(max(dyad_dist_event[after_idxs],na.rm=T) < thresh_h){
      upper <- max(dyad_dist_event[after_idxs],na.rm=T) - .001
    } else{
      upper <- thresh_h
    }
    if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
      lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
    } else{
      lower <- thresh_l
    }
  }
  if(event_type == 'fusion'){
    if(max(dyad_dist_event[before_idxs],na.rm=T) < thresh_h){
      upper <- max(dyad_dist_event[before_idxs],na.rm=T) - .001
    } else{
      upper <- thresh_h
    }
    if(min(dyad_dist_event[middle_idxs],na.rm=T) > thresh_l){
      lower <- min(dyad_dist_event[middle_idxs],na.rm=T) + .001
    } else{
      lower <- thresh_l
    }
  }
  
  
  category <- rep(NA, length(dyad_dist_event))
  category[which(dyad_dist_event < lower)] <- 0
  category[which(dyad_dist_event >= lower & dyad_dist_event < upper)] <- 1
  category[which(dyad_dist_event >= upper)] <- 2
  
  #run length encoding to get sequences of each category
  seqs <- rle(category)
  
  #find sequences of high-middle-low (2,1,0) or low-mid-high (0,1,2)
  seqs_str <- paste0(as.character(seqs$values), collapse = '') #convert to string
  if(event_type=='fission'){
    event_loc <- as.data.frame(str_locate_all(seqs_str,'012')[[1]])
  } 
  if(event_type == 'fusion'){
    event_loc <- as.data.frame(str_locate_all(seqs_str,'210')[[1]])
  }
  
  #for seuqneces of hml or lmh (for fission and fusion respectively), get the time index when they start and end
  if(nrow(event_loc)>0){
    for(r in 1:nrow(event_loc)){
      event_loc$start_time[r] <- ti + sum(seqs$lengths[1:event_loc$start[r]])
      event_loc$end_time[r] <- ti + sum(seqs$lengths[1:event_loc$end[r]-1])
    }
  }
 
  print(nrow(event_loc))
  #get the before and after times
  
  if(plot == T){
    #quartz() #open a new plot for mac
    windows() #for windows
    par(mfrow=c(2,1))
    plot(ti:tf, dyad_dist[ti:tf],type='l', main = paste(event_type, datetime),xlab='Time (min)',ylab = 'Distance apart (m)')
    abline(v=t_event,col='black', lty = 2)
    abline(h = thresh_h, col = 'darkorange1')
    abline(h = thresh_l, col = 'magenta')
    for(r in 1:nrow(event_loc)){
      abline(v=event_loc$start_time[r], col = 'green')
      abline(v=event_loc$end_time[r], col = 'green')
    }
    
    xmin <- min(min(xA[,ti:tf],na.rm=T),min(xB[,ti:tf],na.rm=T))
    xmax <- max(max(xA[,ti:tf],na.rm=T),max(xB[,ti:tf],na.rm=T))
    ymin <- min(min(yA[,ti:tf],na.rm=T),min(yB[,ti:tf],na.rm=T))
    ymax <- max(max(yA[,ti:tf],na.rm=T),max(yB[,ti:tf],na.rm=T))
    plot(NULL, xlim=c(xmin,xmax),ylim=c(ymin,ymax),asp=1, xlab='Easting', ylab = 'Northing', main = paste('(Red =', group_A_names, '), (Blue =', group_B_names,')'))
    for(j in 1:nrow(xA)){
      lines(xA[j,ti:tf],yA[j,ti:tf],type='l',col='#FF000033')
    }
    for(j in 1:nrow(xB)){
      lines(xB[j,ti:tf],yB[j,ti:tf],type='l',col='#0000FF33')
    }
    lines(xcA[ti:tf],ycA[ti:tf], col = '#FF0000', lwd = 2)
    lines(xcB[ti:tf],ycB[ti:tf], col = '#0000FF', lwd = 2)
    
    points(xcA[t_event], ycA[t_event], pch = 8, col = 'black')
    points(xcB[t_event], ycB[t_event], pch = 8, col = 'black')
    points(xcA[ti], ycA[ti], pch = 1, col = 'black')
    points(xcB[ti], ycB[ti], pch = 1, col = 'black')
    points(xcA[tf], ycA[tf], pch = 4, col = 'black')
    points(xcB[tf], ycB[tf], pch = 4, col = 'black')
    
    #algorithm-identified start and end
    points(xcA[event_loc$start_time],ycA[event_loc$start_time],pch = 1, col = 'green')
    points(xcB[event_loc$start_time],ycB[event_loc$start_time],pch = 1, col = 'green')
    points(xcA[event_loc$end_time],ycA[event_loc$end_time],pch = 4, col = 'green')
    points(xcB[event_loc$end_time],ycB[event_loc$end_time],pch = 4, col = 'green')
    
  }
  
  #create an object to output the start and end times
  start_time <- event_loc$start_time
  end_time <- event_loc$end_time
  out <- list(start_time = start_time, end_time = end_time)
  
  #return the output object
  invisible(out)
}







 