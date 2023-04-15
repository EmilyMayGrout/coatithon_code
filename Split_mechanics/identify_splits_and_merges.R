#Identify splits and merges using algorith from Della Libera et al. submitted

library(lubridate)
library(dbscan)

#DIRECTORIES AND PARAMETERS
R_inner <- 15
R_outer <- 50

codedir <- '~/Dropbox/code_ari/coatithon_code/'
dir <- '~/Dropbox/coati/processed/' #directory where all data is stored
group <- 'galaxy' #subdirectory where the group data is stored

#get directory to group data
groupdir <- paste0(dir,group)

#LOAD DATA
#navigate into directory
setwd(codedir)

#read in events
events <- read.csv(paste0('Split_mechanics/',group,'_manual_split_merge_clean.csv'), sep=';')

#read in coati ids
setwd(groupdir)
load(file=paste0(group,'_coati_ids.RData'))

#read in timestamp data
load(file=paste0(group,'_xy_highres_level1.RData'))

#IDENTIFY SUBGROUPS AT EACH TIME POINT
#number of inds and times
n_inds <- nrow(xs)
n_times <- ncol(xs)

#day start indexes
days <- date(ts)
day_start_idxs <- c(1, which(diff(days)==1)+1)

#Get dyadic distances for each pair, then use double threhsold method to determine if htey are otgether at any moment
dyad_dists <- together <- array(NA, dim = c(n_inds, n_inds, n_times))
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    
    #dyadic distance
    dx <- xs[i,] - xs[j,]
    dy <- ys[i,] - ys[j,]
    dyad_dists[i,j,] <- sqrt(dx^2 + dy^2)
    
    #together or not
    #loop over days. for each day...
    for(d in 1:(length(day_start_idxs)-1)){
      
      #get times for that day
      t_day <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
      
      #dyadic distances on that day
      dyad_dists_ij <- dyad_dists[i,j,t_day]
      
      if(sum(!is.na(dyad_dists_ij))==0){
        next
      }
      
      #times when together within inner radius
      together_inner <- dyad_dists_ij <= R_inner
      
      #times when together within outer radius
      together_outer <- dyad_dists_ij <= R_outer
      
      together_ij <- together_inner
      
      if(sum(together_inner,na.rm=T)==0){
        together[i,j,t_day] <- together[j,i,t_day] <- together_ij
        next
      }
      
      #go backwards from crossing points into inner radius to find the 'starts' when crossed the outer radius
      inner_starts <- which(diff(together_inner)==1)+1
      if(length(inner_starts)==0){
        together[i,j,t_day] <- together[j,i,t_day] <- together_ij
        next
      }
      for(k in 1:length(inner_starts)){
        crossing <- inner_starts[k]
        curr_time <- crossing
        for(curr_time in seq(crossing,1,-1)){
          if(is.na(together_outer[curr_time])){
            start <- curr_time + 1
            break
          }
          if(curr_time == 1){
            start <- curr_time
            break
          }
          if(together_outer[curr_time]==F){
            start <- curr_time + 1
            break
          }
          
        }
        together_ij[start:crossing] <- T
      }
      
      #go forwards from crossing points out of outer radius to find the 'ends' when crossed the outer radius
      inner_ends <- which(diff(together_inner)==-1)
      if(length(inner_ends)==0){
        together[i,j,t_day] <- together[j,i,t_day] <- together_ij
        next
      }
      for(k in 1:length(inner_ends)){
        crossing <- inner_ends[k]
        curr_time <- crossing
        for(curr_time in seq(crossing,length(together_ij),1)){
          if(is.na(together_outer[curr_time])){
            end <- curr_time - 1
            break
          }
          if(curr_time == length(together_outer)){
            end <- curr_time
            break
          }
          if(together_outer[curr_time]==F){
            end <- curr_time - 1
            break
          }
          
        }
        together_ij[crossing:end] <- T
      }
      
      together[i,j,t_day] <- together[j,i,t_day] <- together_ij
      
    }
  }
}

# i <- 2
# j <- 8
# d <- 1
# quartz()
# par(mfrow=c(2,1))
# plot(dyad_dists[i,j,day_start_idxs[d]:(day_start_idxs[d+1]-1)],type='l',ylab='',xlab='')
# abline(h = R_inner, col = 'red')
# abline(h = R_outer, col = 'red')
# plot(together[i,j,day_start_idxs[d]:(day_start_idxs[d+1]-1)], type = 'l',ylim=c(0,1),ylab='',xlab='')

#Identify groups from together matrices
groups <- matrix(NA, nrow = n_inds, ncol = n_times)
for(t in 1:n_times){
  non.nas <- which(colSums(!is.na(together[,,t]))>0)
  if(length(non.nas)<=1){
    next
  }
  non.nas.together <- together[non.nas,non.nas,t]
  diag(non.nas.together) <- 1
  grps.non.nas <- dbscan(x = as.dist(1 - non.nas.together), eps = .1,minPts=1)$cluster
  groups[non.nas, t] <- grps.non.nas
}

#store groups as lists of lists
groups_list <- list()
for(t in 1:n_times){
  
  #create a list to hold the groups at that timestep
  groups_list[[t]] <- list()
  if(sum(!is.na(groups[,t]))==0){
    next
  }
  
  #
  max_group_id <- max(groups[,t],na.rm=T)
  for(i in seq(1,max_group_id)){
    groups_list[[t]][[i]] <- which(groups[,t]==i)
  }
  
}

#IDENTIFY EVENTS WHERE SOMETHING CHANGES
event_times <- c()
for(d in 1:(length(day_start_idxs)-1)){
  t_day <- day_start_idxs[d]:(day_start_idxs[d+1]-1)
  for(t in t_day[1:length(t_day)]){
    groups_curr <- groups_list[[t]]
    groups_next <- groups_list[[t+1]]
    if(sum(!is.na(groups_curr))==0 | sum(!is.na(groups_next))==0){
      next
    }
    n_groups_curr <- length(groups_curr)
    n_groups_next <- length(groups_next)
    for(i in 1:n_groups_curr){
      group_curr <- groups_curr[[i]]
      matched <- F
      for(j in 1:n_groups_next){
        group_next <- groups_next[[j]]
        if(setequal(group_curr,group_next)){
          matched <- T
        }
      }
    }
    if(!matched){
      event_times <- c(event_times, t)
    }
  }
}


#go through event times and classify events into types
changes <- data.frame(tidx = event_times)
changes$datetime <- ts[changes$tidx]
changes$event_type <- NA
changes$n_groups_curr <- unlist(lapply(groups_list, length))[event_times]
changes$n_groups_next <- unlist(lapply(groups_list, length))[event_times+1]
changes$n_inds_curr <- changes$n_inds_next <- NA
for(i in 1:nrow(changes)){
  t <- changes$tidx[i]
  groups_curr <- groups_list[[t]]
  groups_next <- groups_list[[t+1]]
  changes$n_inds_curr[i] <- sum(unlist(lapply(groups_curr, length)))
  changes$n_inds_next[i] <- sum(unlist(lapply(groups_next, length)))
}

changes$event_type[which(changes$n_groups_curr < changes$n_groups_next & changes$n_inds_curr==changes$n_inds_next)] <- 'fission'
changes$event_type[which(changes$n_groups_curr > changes$n_groups_next & changes$n_inds_curr==changes$n_inds_next)] <- 'fusion'

#get groups for fissions and fusions
changes$group_A_idxs <- changes$group_A <- changes$group_B_idxs <- changes$group_B <- list(c(0))
fissions <- which(changes$event_type=='fission')
fusions <- which(changes$event_type=='fusion')
for(i in fissions){
  
  groups_curr <- groups_list[[changes$tidx[i]]]
  groups_next <- groups_list[[changes$tidx[i]+1]]
  n_groups_next <- changes$n_groups_next[i]
  n_groups_curr <- changes$n_groups_curr[i]
  
  #if there are only 2 groups at the end, then name them A and B
  if(n_groups_next==2){
    changes$group_A_idxs[i] <- groups_next[1]
    changes$group_B_idxs[i] <- groups_next[2]
  }
  
  #if there are more than 2 groups, need to figure out which one changed
  if(n_groups_next > 2){
    matched_groups <- c()
    for(j in 1:n_groups_next){
      matched <- F
      for(k in 1:n_groups_curr){
        if(setequal(groups_next[[j]], groups_curr[[k]])){
          matched <- T
        }
      }
      if(matched){
        matched_groups <- c(matched_groups, j)
      }
    }
    unmatched_groups <- setdiff(1:n_groups_next, matched_groups)
    if(length(unmatched_groups)==2){
      changes$group_A_idxs[i] <- groups_next[unmatched_groups[1]]
      changes$group_B_idxs[i] <- groups_next[unmatched_groups[2]]
    }else{
      print(i)
    }
    
  }
}

for(i in fusions){
  
  groups_curr <- groups_list[[changes$tidx[i]]]
  groups_next <- groups_list[[changes$tidx[i]+1]]
  n_groups_next <- changes$n_groups_next[i]
  n_groups_curr <- changes$n_groups_curr[i]
  
  #if there are only 2 groups at the start, then name them A and B
  if(n_groups_curr==2){
    changes$group_A_idxs[i] <- groups_curr[1]
    changes$group_B_idxs[i] <- groups_curr[2]
  }
  
  #if there are more than 2 groups, need to figure out which one changed
  if(n_groups_curr > 2){
    matched_groups <- c()
    for(j in 1:n_groups_curr){
      matched <- F
      for(k in 1:n_groups_next){
        if(setequal(groups_curr[[j]], groups_next[[k]])){
          matched <- T
        }
      }
      if(matched){
        matched_groups <- c(matched_groups, j)
      }
    }
    unmatched_groups <- setdiff(1:n_groups_curr, matched_groups)
    if(length(unmatched_groups)==2){
      changes$group_A_idxs[i] <- groups_curr[unmatched_groups[1]]
      changes$group_B_idxs[i] <- groups_curr[unmatched_groups[2]]
    }else{
      print(i)
    }
    
  }
}

#get coati names
for(i in c(fissions,fusions)){
  changes$group_A[i] <- list(coati_ids$name_short[changes$group_A_idxs[i][[1]]])
  changes$group_B[i] <- list(coati_ids$name_short[changes$group_B_idxs[i][[1]]])
}

#final events dataframe (call it events_detected)
events_detected <- changes[c(fissions, fusions),c('tidx','datetime','event_type','group_A_idxs','group_B_idxs','group_A','group_B')]

#get number of individuals in each subgroup
events_detected$n_A <- sapply(events_detected$group_A_idxs, length)
events_detected$n_B <- sapply(events_detected$group_B_idxs, length)

#non-single-individual events
non_singles <- which(events_detected$n_A != 1 & events_detected$n_B!=1)

#visualize events
analyse_ff_event(14, events = events_detected, xs, ys)

