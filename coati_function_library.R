#library of general functions

#LIBRARIES
library(dbscan)
library(rgdal)
library(lubridate)


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

#Converts a matrix of eastings and northings (eastings first column, northings second column) to UTM
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