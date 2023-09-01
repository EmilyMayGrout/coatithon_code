#analyse merge events (opposite to the code for getting split events)

#NOTE: this code only works for up to 3 subgroups per split
#TODO: if more than 3 subgroups present, need to generalize

data_dir <- "/home/pranav/Personal/Temp/emily/Data/galaxy/"
code_dir <- '/home/pranav/Personal/Temp/emily/code/code_review/'
plot_dir <- '/home/pranav/Personal/Temp/emily/Figures/galaxy/'
gps_file <- "galaxy_xy_10min_level1.RData"
id_file <- 'galaxy_coati_ids.RData' 

library(fields)
library(viridis)
library(hms)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, 50)


#run through time by time and check if each time is a merge
merge <- c()

for(t in 1:(n_times-1)){
  
  #get subgroup membership now and previous
  merge_group <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(merge_group))==0 | sum(!is.na(subgroups_previously))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  merge_group[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(merge_group))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_merge_group <- length(unique(merge_group[!is.na(merge_group)]))
  n_subgroups_previously <- length(unique(subgroups_previously[!is.na(subgroups_previously)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(merge_group)==1)
  singletons_later <- sum(table(subgroups_previously)==1)
  
  #determine if this time step is a merge
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if(n_merge_group==1 
     & n_subgroups_previously >1
     & singletons_later==0
  ){
    merge <- c(merge, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_merge_group > 1 
     & (singletons_now + 1) == n_merge_group
     & n_subgroups_previously > n_merge_group
     & singletons_now == singletons_later
  ){
    merge <- c(merge, t)
  }
  
  
}

#this seems to work



#make a data frame of merges, with merged group and subgroups
merge_df <- data.frame(t = merge, merge_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)


i=5
for(i in 1:nrow(merge_df)){
  t <- merge_df$t[i]
  merge_group <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #transfer NAs from now to later and later to now
  merge_group[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(merge_group))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(merge_group)
  orig_subgroup_members <- which(merge_group == orig_subgroup_number)
  
  #store original group membership in data frame
  merge_df$merge_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_previously[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_previously==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
    #really hacky shit to get R to put lists into a data frame :(
    if(j==1){
      merge_df$sub1[i] <- list(orig_inds_in_group)
    }
    if(j==2){
      merge_df$sub2[i] <- list(orig_inds_in_group)
    }
    if(j==3){
      merge_df$sub3[i] <- list(orig_inds_in_group)
    }
    if(j==4){
      merge_df$sub4[i] <- list(orig_inds_in_group)
    }
    if(j==5){
      merge_df$sub5[i] <- list(orig_inds_in_group)
    }
  }
}

#number in each subgroup
merge_df$n_merge <- sapply(merge_df$merge_group, function(x){return(sum(!is.na(x)))})
merge_df$n_sub1 <- sapply(merge_df$sub1, function(x){return(sum(!is.na(x)))})
merge_df$n_sub2 <- sapply(merge_df$sub2, function(x){return(sum(!is.na(x)))})
merge_df$n_sub3 <- sapply(merge_df$sub3, function(x){return(sum(!is.na(x)))})

save(merge_df, file = paste0(data_dir, "merge_df.Rdata"))


#DONE WITH DATAFRAME!


#now look at time difference between merges and splits
#open splits_df
load(paste0(data_dir, "splits_df.Rdata"))

#luckily the number of merges and split events is 29
time_diff <- data.frame(splits_t = splits_df$t, merge_t = merge_df$t)
time_diff$diff <- time_diff$merge - time_diff$splits_t
time_diff$diff_time_hour <- (time_diff$diff*10)/60
time_diff$diff_time_hour <-  format(round(time_diff$diff_time_hour, 1), nsmall = 1)
time_diff$diff_time_hour <- as.numeric(time_diff$diff_time_hour)

#running to make a graph with presedente (code in merge_analysis_presedente)
time_diff_gal <- time_diff

png(file = paste0(plot_dir, "split_duration_hist.png"), width = 1000, height = 600, units = "px")

par(mar=c(8, 8, 8, 8))
hist(time_diff$diff_time_hour, breaks = 80, xlab = "Time between splits and merges (hours)", col = "lightblue3", main = "", cex.lab = 2,cex.axis = 1.5)
dev.off()

#------------------------------------------------------------------
#get latlon coords for merge events

#get the index of when there was a split
merge_indx <- merge_df$t

#find the x and y UTM coords of those splits
merge_xs <- xs[,merge_indx]
merge_ys <- ys[,merge_indx]

#make this into a dataframe to reshape it so its not a matrix and that we can join the x and y coords together for the mapping
merge_xs <- as.data.frame(merge_xs) 
merge_xs <- reshape(merge_xs, varying=1:ncol(merge_xs), v.names="xs", direction ="long", idvar = "ID")

merge_ys <- as.data.frame(merge_ys)
merge_ys <- reshape(merge_ys, varying=1:ncol(merge_ys), v.names="xs", direction ="long", idvar = "ID")

merge_xy <- cbind(merge_xs, merge_ys)
merge_xy <- merge_xy[,-c(3,4)]
colnames(merge_xy) <- c("event", "xs","ys", "ID")

#quick plot of merge locations 
plot(merge_xy$xs, merge_xy$ys, col = merge_xy$ID) 

#remove id column
merge_xy_1 <- merge_xy[,-c(1,4)]

#convert to latlon
merge_latlon <- as.dlibata.frame(utm.to.latlon(merge_xy_1, utm.zone = '17',southern_hemisphere=FALSE))
merges_utm_latlon <- cbind(merge_latlon, merge_xy)
merges_utm_latlon$event <- as.factor(merges_utm_latlon$event)

save(merges_utm_latlon, file = paste0(data_dir, "merges_utm_latlon.RData")) 






