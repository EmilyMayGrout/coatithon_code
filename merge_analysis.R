#analyse split events
#need to adapt this code to look at merge events
#Ari said to swap (t and t+1) for (t and t-1)
#I ran it with these changes but it doesn't seem to be correct, still need to play with it

#NOTE: this code only works for up to 3 subgroups per split
#TODO: if more than 3 subgroups present, need to generalize

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData' 

library(fields)
library(viridis)

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


#run through time by time and check if each time is a split
splits <- c()
for(t in 1:(n_times-1)){
  
  #get subgroup membership now and previous
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(subgroups_now))==0 | sum(!is.na(subgroups_previously))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(subgroups_now))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
  n_subgroups_previously <- length(unique(subgroups_previously[!is.na(subgroups_previously)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(subgroups_now)==1)
  singletons_later <- sum(table(subgroups_previously)==1)
  
  #determine if this time step is a merge
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if(n_subgroups_now==1 
     & n_subgroups_previously >1
     & singletons_later==0
  ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & (singletons_now+1) == n_subgroups_now
     & n_subgroups_previously > n_subgroups_now
     & singletons_now == singletons_later
  ){
    splits <- c(splits, t)
  }
  
  
}

#make a data frame of splits, with original group and subgroups
splits_df <- data.frame(t = splits, orig_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)
for(i in 1:nrow(splits_df)){
  t <- splits_df$t[i]
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_previously <- subgroup_data$ind_subgroup_membership[,t-1]
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_previously))] <- NA
  subgroups_previously[which(is.na(subgroups_now))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(subgroups_now)
  orig_subgroup_members <- which(subgroups_now == orig_subgroup_number)
  
  #store original group membership in data frame
  splits_df$orig_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_previously[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_previously==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
    #really hacky shit to get R to put lists into a data frame :(
    if(j==1){
      splits_df$sub1[i] <- list(orig_inds_in_group)
    }
    if(j==2){
      splits_df$sub2[i] <- list(orig_inds_in_group)
    }
    if(j==3){
      splits_df$sub3[i] <- list(orig_inds_in_group)
    }
    if(j==4){
      splits_df$sub4[i] <- list(orig_inds_in_group)
    }
    if(j==5){
      splits_df$sub5[i] <- list(orig_inds_in_group)
    }
  }
}

#number in each subgroup
splits_df$n_orig <- sapply(splits_df$orig_group, function(x){return(sum(!is.na(x)))})
splits_df$n_sub1 <- sapply(splits_df$sub1, function(x){return(sum(!is.na(x)))})
splits_df$n_sub2 <- sapply(splits_df$sub2, function(x){return(sum(!is.na(x)))})
splits_df$n_sub3 <- sapply(splits_df$sub3, function(x){return(sum(!is.na(x)))})

#save(splits_df, file = "C:/Users/egrout/Dropbox/coatithon/coatithon_code/splits_on_map/splits_df.Rdata")  

#DONE WITH DATAFRAME!

p_dyad_together <- get_p_dyad_together(splits_df_local = splits_df, n_inds_local = n_inds)

#Compute consistency based on consistency metric (see coati_function_library)
consistency_data <- get_consistency(p_dyad_together)

#randomize splits and recompute consistency metric
n_rands <- 1000
rando_consistencies <- rep(NA, n_rands)
for (i in 1:n_rands){
  
  rando <- randomise_splits(splits_df)
  rando_dyad <- get_p_dyad_together(splits_df_local = rando, n_inds_local = n_inds)
  consist <- get_consistency(rando_dyad)
  rando_consistencies[i] <- consist
  
}

png(height = 600, width = 900, units = 'px', filename = paste0(plot_dir,'consistency_splits_hist.png'))
par(mfrow=c(1,1), mar = c(6,6,3,3))#(bottom, left, top, right)
hist(rando_consistencies, breaks=seq(0,0.5,0.005), main = "", xlab = "Consistency of random sub-group allocations", col = "slategray4",  cex.lab = 2.5, cex.axis=2.5)
abline(v=consistency_data, col = 'orange2', lwd=4)
dev.off()


#this graph shows that the mean consistency of individuals splitting with others is more likely than randomly assigning individuals into sub-groups

#should look into repeatability of binary data - to see if there is a better metric for getting consistency values from binary data

#symmetrize matrix - copy entries from p_dyad_together[i,j] to p_dyad_together[j,i]
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    p_dyad_together[j,i] <- p_dyad_together[i,j]
  }
}


diag(p_dyad_together) <- NA
new_order <- c(5,1,11,4,10,2,3,6,7,8,9)
p_dyad_together_reorder <- p_dyad_together[new_order, new_order]

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network_splits.png'))
par(mfrow=c(1,1), mar = c(1,2,1,1))#(bottom, left, top, right)
visualize_network_matrix(p_dyad_together_reorder, coati_ids[new_order,])
dev.off()
#image.plot(p_dyad_together)

