#analyse merge events (opposite to the code for getting split events)

#NOTE: this code only works for up to 3 subgroups per split
#TODO: if more than 3 subgroups present, need to generalize

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/'
gps_file <- "presedente_xy_10min_level0.RData"
id_file <- 'coati_ids.RData' 

library(fields)
library(viridis)
library(hms)
library(dplyr)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

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

save(merge_df, file = "C:/Users/egrout/Dropbox/coatithon_notgithub/results/merge_results/Presedente/merge_df.Rdata")  




load("C:/Users/egrout/Dropbox/coatithon_notgithub/Presedente_fission_fusion/splits_df.Rdata")

#number of splits and merges is not the same, so should rbind them

merge_df$event <- "merge"
splits_df$event <- "split"

merge_df_subset <- merge_df[,c(1,12)]
splits_df_subset <- splits_df[,c(1,12)]

time_diff <- rbind(splits_df_subset,merge_df_subset)
time_diff <- time_diff %>% arrange(-desc(t))

#because the number of splits and merges are different, I'm filtering to times when its a merge then a split (so the duration of a split is from the first split to the next merge)
i = 1
time_diff$keep <- NA
time_diff$keep2 <- NA

for (i in 1:nrow(time_diff)){
  
  if (time_diff$event[i] == "merge" & time_diff$event[i+1] == "split"){
    time_diff$keep[i] <- "1"
    time_diff$keep2[i+1] <- "1"
  } else {
    time_diff$keep[i] <- "0"
    } 
  }
#replace NA's with 0 so I can sum the results 
time_diff$keep2[is.na(time_diff$keep2)] <- 0 
time_diff$keep3 <- as.numeric(time_diff$keep) + as.numeric(time_diff$keep2)
time_diff <- time_diff[time_diff$keep3 == 1,]
time_diff <- time_diff[,c(1,2)]
splits_time_diff <- time_diff[time_diff$event == "split",]
merge_time_diff <- time_diff[time_diff$event == "merge",]

time_diff <- data.frame(splits_t = splits_time_diff$t, merge_t = merge_time_diff$t)
time_diff$diff <- time_diff$splits_t - time_diff$merge_t
time_diff$diff_time_hour <- (time_diff$diff*10)/60
time_diff$diff_time_hour <-  format(round(time_diff$diff_time_hour, 1), nsmall = 1)
time_diff$diff_time_hour <- as.numeric(time_diff$diff_time_hour)

#to make a graph with galaxy distributions (running the merge_analysis of galaxy)
time_diff_pres <- time_diff

#time_diff_gal was made in the merge_analysis script
hist(time_diff_gal$diff_time_hour, col='green', add=TRUE)

png(file = "C:/Users/egrout/Dropbox/coatithon_notgithub/results/merge_results/Presedente/split_duration_hist.png", width = 1000, height = 600, units = "px")

par(mar=c(8, 8, 8, 8))
hist(time_diff_pres$diff_time_hour, breaks = 80, xlab = "Time between splits and merges (hours)", col = "lightblue3", main = "", cex.lab = 2,cex.axis = 1.5)
#time_diff_gal was made in the merge_analysis script
dev.off()



#making a plot for both galaxy and presedente overlapping (need to run the merge_analysis for galaxy for this to work)
png(file = "C:/Users/egrout/Dropbox/coatithon_notgithub/results/merge_results/split_duration_hist_both1.png", width = 1000, height = 600, units = "px")
par(mar=c(8, 8, 8, 8))
hist(time_diff_gal$diff_time_hour, col='lightblue4', breaks = 80, xlab = "Time between splits and merges (hours)", main = "", cex.lab = 2,cex.axis = 1.5, ylim = c(0, 8))
hist(time_diff_pres$diff_time_hour, breaks = 10 ,col=rgb(0.7,1,1,0.5),add=TRUE, alpha = 0.5)
legend(x = "topright",
       legend = c("Galaxy", "Presedente"),  
       lty = c(1, 1),          
       col = c('lightblue4', rgb(0.7,1,1)), lwd = 5, cex=2,  bty = "n")                
dev.off()







