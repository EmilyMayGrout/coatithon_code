#fission fusion analysis script
#first define groups for each time step (each 10 mins)

#--------priority list----------
#how the radius size influences the distribution of the number of sub groups over time
#who tends to be on their own - plot each individuals subgroup sizes  
#when they split, the distribution of sub group sizes 

#network of within group spatial associations and is that consistent over time
#network of fraction of time each pair is in the same sub group and is it consistent over time
#are these two things related? does the within group network correlate with the fission-fusion subgroup network


library(dbscan)

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
in_file <- "galaxy_xy_10min_level0.RData"

#radius for group membership
R <- 50

#-------MAIN-------
setwd(data_dir)

load(in_file)

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
ngroups <- apply(sub_groups, 2, max, na.rm = T)
#remove infinite values
ngroups[is.infinite(ngroups)] <- NA
hist(ngroups)
plot(ngroups, type = "l")


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



 