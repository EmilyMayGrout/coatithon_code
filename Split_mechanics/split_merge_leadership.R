#Script to look at leadership in split and merge events
#Q1: Who leaves and who stays
#Q2: Who leads the movement of the group when leaving or joining

#-----PARAMETERS-------

user <- 'emily'
group <- 'galaxy'
use_manual_events <- F
dist_moved_thresh <- 15 #minimum distance moved by a subgroup to count it as having moved (i.e. left or joined)
make_plots <- F
dist_frac_thresh <- 0.5
n_rands <- 1000
own_finish_line <- T

#bins for the normalized ranks in the entropy computation
#default is seq(0,1,.2) which is 5 bins of size 0.2
#if you want to look at front vs non-front, can set bins to c(0,.99,1) - this will give 2 bins, one running from 0 to .99 and the other giving only those with value 1
norm_rank_bins <- c(0,.99,1)

if(user=='ari'){
  groupdir <- paste0('~/Dropbox/coati/processed/', group)
} else{
  if(group == 'galaxy'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
  } else if(group == 'presedente'){
    groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
  }
}

#-----------GRAPHICS--------
#make compatible with windows OS
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

#FUNCTION
#Measure the displacement of each individual in a subgroup along that subgroup's movement path during a split / merge
#INPUTS:
#moving_inds: indexes of the individuals in the moving subgroup
#xs, ys: matrices of positions
#t0: start time index of event 
#tf: end time index of event - used for computing group direction of movement
#tmeas: time index at which you measure the individual relative to the group's heading (normally this would also be tf)
#OUTPUTS:
#proj_i: [vector] displacement of each individual along the group direction of motion (projection)
ind_disp_along_group_path <- function(moving_inds, xs, ys, t0, tf, tmeas){
  
  #get centroid initial location
  xc0 <- mean(xs[moving_inds,t0], na.rm = T)
  yc0 <- mean(ys[moving_inds,t0], na.rm = T)
  
  #centroid final location
  xcf <- mean(xs[moving_inds,tf], na.rm = T)
  ycf <- mean(ys[moving_inds,tf], na.rm = T)
  
  #centroid displacement vector
  dxc <- xcf - xc0
  dyc <- ycf - yc0
  
  #get individual end points
  xif <- xs[moving_inds, tmeas]
  yif <- ys[moving_inds, tmeas]
  
  #get individual displacement vectors (their end point to the group start point)
  dxi <- xif - xc0
  dyi <- yif - yc0
  
  #project individual displacement vectors onto the group displacement vector
  disp_i <- (dxi*dxc + dyi*dyc) / sqrt(dxc^2 + dyc^2)
  
  #get ranks
  if(sum(is.na(disp_i))==0){
    ranks <- rank(disp_i)
  } else{
    ranks <- rep(NA, length(disp_i))
  }
  
  #get normalized ranks
  if(length(disp_i)==1){
    norm_ranks <- c(NA)
  } else{
    min_rank <- min(ranks)
    max_rank <- max(ranks)
    norm_ranks <- (ranks - min_rank) / (max_rank - min_rank)
  }
  
  #return ind displacement along group vector, ranks, and normalized ranks
  out <- list()
  out$disp <- disp_i
  out$ranks <- ranks
  out$norm_ranks <- norm_ranks
  return(out) 
}

#Function to compute the order at which individuals cross a threshold distance along the fission/fusion direction
#INPUTS:
#moving_inds: indexes of the individuals in the moving subgroup
#xs, ys: matrices of positions
#t0: start time index of event 
#tf: end time index of event - used for computing group direction of movement
#dist_frac_thresh: threshold fractional distance along group trajectory from start point to compute passing time for each individual
#own_finish_line: defaults to F, in which case the individual ranks are computed in the normal way. If T, 
# then use the "everyone has their own finish line" version, where we check the times at which each individual has moved a distance of
# dist_frac_thresh*total group displacement along the group trajectory, relative to their initial starting location
ind_crossing_thresh_times_along_group_path <- function(moving_inds, xs, ys, t0, tf, dist_frac_thresh = 0.5, own_finish_line = F){
  
  #if times are missing, return NAs for ranks and crossing times
  if(is.na(t0) | is.na(tf)){
    out <- list()
    out$first_crossing_times <- rep(NA, length(moving_inds))
    out$ranks <- rep(NA, length(moving_inds))
    out$norm_ranks <- rep(NA, length(moving_inds))
    return(out)
  }
  
  #get centroid initial location
  xc0 <- mean(xs[moving_inds,t0], na.rm = T)
  yc0 <- mean(ys[moving_inds,t0], na.rm = T)
  
  #centroid final location
  xcf <- mean(xs[moving_inds,tf], na.rm = T)
  ycf <- mean(ys[moving_inds,tf], na.rm = T)
  
  #centroid displacement vector
  dxc <- xcf - xc0
  dyc <- ycf - yc0
  
  #centroid total distance traveled
  distc <- sqrt(dxc^2 + dyc^2)
  
  #if distance moved was zero, return NAs
  if(distc == 0){
    out <- list()
    out$first_crossing_times <- rep(NA, length(moving_inds))
    out$ranks <- rep(NA, length(moving_inds))
    out$norm_ranks <- rep(NA, length(moving_inds))
    return(out)
  }
  
  #get the distance threshold (along group trajectory) used to determine order of crossing
  dist_thresh <- dist_frac_thresh * distc
  
  #get individual end points
  xi <- xs[moving_inds, t0:ncol(xs)]
  yi <- ys[moving_inds, t0:ncol(xs)]
  
  #get individual displacement vectors (their end point to the group start point)
  dxi <- xi - xc0
  dyi <- yi - yc0
  
  #project individual displacement vectors onto the group displacement vector
  disp_i <- (dxi*dxc + dyi*dyc) / sqrt(dxc^2 + dyc^2)
  
  #if own_finish_line is F, use the normal metric - we get the order in which each individaul
  #crossed a (shared) threshold line along the group displacement trajectory
  #if own_finish_line is T, we get the order in which each individual
  #has moved a distance of dist_thresh along the group direction, from its original starting point
  if(own_finish_line){
    if(length(moving_inds) > 1){
      disp_i <- disp_i - disp_i[,1]
    } else{
      disp_i <- disp_i - disp_i[1]
    }
  }
  
  #get times of each individual crossing the threshold
  first_crossing_times <- rep(NA, length(moving_inds))
  for(ind in 1:length(moving_inds)){
    
    if(is.matrix(disp_i)){
      crossing_times <- which(disp_i[ind,] > dist_thresh)
    } else{
      crossing_times <- which(disp_i > dist_thresh)
    }
    
    #if it never crosses, then give infinity for the crossing time and throw a warning
    if(length(crossing_times)==0){
      first_crossing_times[ind] <- Inf
      
      warning('Individual never crossed the distance threshold - crossing time set to Inf')
    }else{
      first_crossing_times[ind] <- min(crossing_times)
    }
  }
  #get ranks
  if(sum(is.na(first_crossing_times))==0){
    ranks <- rank(-first_crossing_times) #rank negative value so that higher ranks are higher leadership
  } else{
    ranks <- rep(NA, length(first_crossing_times))
  }
  
  #get normalized ranks
  if(length(first_crossing_times)==1){
    norm_ranks <- c(NA)
  } else{
    min_rank <- min(ranks)
    max_rank <- max(ranks)
    norm_ranks <- (ranks - min_rank) / (max_rank - min_rank)
  }
  
  #return ind displacement along group vector, ranks, and normalized ranks
  out <- list()
  out$first_crossing_times <- first_crossing_times
  out$ranks <- ranks
  out$norm_ranks <- norm_ranks
  
  return(out) 
  
}

#Compute entropy from a set of measurements
#histo is a vector of frequencies or probabilities to take the entropy of
compute_entropy <- function(histo){
  
  #don't allow negative values
  if(sum(histo < 0)>0){
    stop('negative values found')
  }
  
  #if NAs, return NA
  if(sum(is.na(histo))>0){
    return(NA)
  }
  
  #if frequencies rather than probabilities, normalize to add up to 1
  histo <- histo / sum(histo) 
  
  #find nonzero elements
  nonzero <- which(histo!=0)
  
  #get entropy, excluding zero elements
  entropy <- -sum(histo[nonzero] * log(histo[nonzero], base = 2))
  
  return(entropy)
}

#Function to get information about leadership (determined by position) during fissions and fusions
#INPUTS:
# events data frame
# n_inds: number of individuals
# leadership_type: 'position' for using position along group movement vector at a particular time 
#     or 'crosstime' for using crossing times
#     or 'crosstime_ownfinishline' for using the crossing times where each individual has its own finish line
# meas_time: where to measure leadership, at the start, midpoint (mid) or end of the event
#OUTPUT:
# out: list containing
# out$fission_leaders: matrix where rows are events and columns are individuals, giving the normalized leadership ranks for fisisons
# out$fusion_leaders: same but for fusions
# out$fission_leadership_entropies: vector of entropies of the normalized rank distributions for fissions for each individual
# out$fusion_leadership_entropies: same but for fusions
get_fission_fusion_leadership <- function(events, n_inds, leadership_type = 'position', meas_time = 'end', norm_rank_bins = seq(0,1,.2)){
  
  #check if meas_time is one of the correct options
  if(!(meas_time %in% c('start','mid','end'))){
    stop('meas_time not specified as start, mid, or end')
  }
  
  if(!(leadership_type %in% c('position','crosstime','crosstime_ownfinishline'))){
    stop('leadership_type not specified as position or crosstime or crosstime_ownfinishline')
  }
  
  #get normalized ranks for individuals during fissions and fusions
  fission_leaders <- fusion_leaders <- matrix(NA, nrow = n_inds, ncol = nrow(events))
  for(i in 1:nrow(events)){
    
    if(!is.na(events$A_moved[i])){
      if(events$A_moved[i]){
        inds <- events$group_A_idxs[i][[1]]
        if(leadership_type=='position'){
          if(meas_time == 'start'){
            norm_ranks <- events$group_A_lead_pos_norm_rank_start[i][[1]] 
          } 
          if(meas_time == 'mid'){
            norm_ranks <- events$group_A_lead_pos_norm_rank_mid[i][[1]] 
          }
          if(meas_time == 'end'){
            norm_ranks <- events$group_A_lead_pos_norm_rank_end[i][[1]] 
          }
        } else if(leadership_type == 'crosstime'){
          norm_ranks <- events$group_A_lead_crosstime_norm_rank[i][[1]]
        } else if(leadership_type == 'crosstime_ownfinishline'){
          norm_ranks <- events$group_A_lead_crosstime_ownfinishline_norm_rank[i][[1]]
        } 
        if(events$event_type[i] == 'fission'){
          fission_leaders[inds,i] <- norm_ranks
        }
        if(events$event_type[i] == 'fusion'){
          fusion_leaders[inds,i] <- norm_ranks
        } 
      }
    }
    
    if(!is.na(events$B_moved[i])){
      if(events$B_moved[i]){
        inds <- events$group_B_idxs[i][[1]]
        if(leadership_type == 'position'){
          if(meas_time == 'start'){
            norm_ranks <- events$group_B_lead_pos_norm_rank_start[i][[1]] 
          } 
          if(meas_time == 'mid'){
            norm_ranks <- events$group_B_lead_pos_norm_rank_mid[i][[1]] 
          }
          if(meas_time == 'end'){
            norm_ranks <- events$group_B_lead_pos_norm_rank_end[i][[1]] 
          }
        } else if(leadership_type =='crosstime'){
          norm_ranks <- events$group_B_lead_crosstime_norm_rank[i][[1]]
        } else if(leadership_type =='crosstime_ownfinishline'){
          norm_ranks <- events$group_B_lead_crosstime_ownfinishline_norm_rank[i][[1]]
        }
        if(events$event_type[i] == 'fission'){
          fission_leaders[inds,i] <- norm_ranks
        }
        if(events$event_type[i] == 'fusion'){
          fusion_leaders[inds,i] <- norm_ranks
        } 
      }
    }
  }
  
  #For each individual get entropy of its leadership normalized rank distribution
  #for fusions
  fusion_leadership_entropies <- fission_leadership_entropies <- rep(NA, n_inds)
  for(i in 1:n_inds){
    if(sum(!is.na(fusion_leaders[i,]))>0){
      histo_fusion <- hist(fusion_leaders[i,], breaks= norm_rank_bins, plot = F)$counts
      fusion_leadership_entropies[i] <- compute_entropy(histo_fusion)
    }
    if(sum(!is.na(fission_leaders[i,]))>0){
      histo_fission <- hist(fission_leaders[i,], breaks= norm_rank_bins, plot = F)$counts
      fission_leadership_entropies[i] <- compute_entropy(histo_fission)
    }
  }
  
  out <- list()
  out$fission_leaders <- fission_leaders
  out$fusion_leaders <- fusion_leaders
  out$fission_leadership_entropies <- fission_leadership_entropies
  out$fusion_leadership_entropies <- fusion_leadership_entropies
  
  return(out)
}

#-------------------MAIN-------------------------

#LOAD EVENTS DATA
setwd(groupdir)

if(use_manual_events){
  load(paste0(group,'_manual_ff_events_characterized.RData'))
} else{
  load(paste0(group,'_auto_ff_events_characterized.RData'))
}

#LOAD GPS AND IDS DATA
#navigate into directory
setwd(groupdir)

#read in coati ids
load(file=paste0(group,'_coati_ids.RData'))

#modify coati ids to only include first 3 letters
coati_ids$name_short <- sapply(coati_ids$name, function(x){return(substr(x,1,3))})

#read in timestamp data
load(file=paste0(group,'_xy_highres_level1.RData'))

#number of individuals
n_inds <- nrow(coati_ids)

#Define whether groups moved or stayed during the split / merge using the dist_moved_thresh threshold
events$A_moved <- events$A_during_disp > dist_moved_thresh
events$B_moved <- events$B_during_disp > dist_moved_thresh

#For each individual compute the number of times they moved vs didn't move, for fissions and fusions separately
n_fission_moves <- n_fissions <- n_fusion_moves <- n_fusions <- rep(0, n_inds)
for(i in 1:nrow(events)){
  
  #get the info associated with that event
  row <- events[i,]
  
  #get the individuals in each subgroup
  group_A <- row$group_A_idxs[[1]]
  group_B <- row$group_B_idxs[[1]]
  
  if(row$event_type=='fission'){
    
    #add 1 to the elements associated with individuals in each group, for the total number of fission events
    n_fissions[group_A] <- n_fissions[group_A] + 1
    n_fissions[group_B] <- n_fissions[group_B] + 1
    
    if(!is.na(row$A_moved)){
      if(row$A_moved == T){
        n_fission_moves[group_A] <- n_fission_moves[group_A] + 1
      }
    }
    if(!is.na(row$B_moved)){
      if(row$B_moved == T){
        n_fission_moves[group_B] <- n_fission_moves[group_B] + 1
      }
    }
    
  } else if(row$event_type == 'fusion'){
    n_fusions[group_A] <- n_fusions[group_A] + 1
    n_fusions[group_B] <- n_fusions[group_B] + 1
    
    if(!is.na(row$A_moved)){
      if(row$A_moved == T){
        n_fusion_moves[group_A] <- n_fusion_moves[group_A] + 1
      }
    }
    if(!is.na(row$B_moved)){
      if(row$B_moved == T){
        n_fusion_moves[group_B] <- n_fusion_moves[group_B] + 1
      }
    }
  }
}

#Get the fraction of time each individual moved during a fission and fusion
fusion_move_fracs <- n_fusion_moves / n_fusions
fission_move_fracs <- n_fission_moves / n_fissions

#get error bars (Clopper Pearson intervals 95%)
fission_CIs <- fusion_CIs <- matrix(NA, nrow = n_inds, ncol = 2)
for(i in 1:n_inds){
  
  if(n_fissions[i]>0){
    test_out <- binom.test(n_fission_moves[i], n_fissions[i])
    fission_CIs[i,1] <- test_out$conf.int[1]
    fission_CIs[i,2] <- test_out$conf.int[2]
  }
  
  if(n_fusions[i]>0){
    test_out <- binom.test(n_fusion_moves[i], n_fusions[i])
    fusion_CIs[i,1] <- test_out$conf.int[1]
    fusion_CIs[i,2] <- test_out$conf.int[2]
  }
}

#LEADERSHIP DURING SPLITS
#TODO: add analysis of start vs. middle vs. end leadership
#get the displacement of each individual project along the subgroup vector, for each subgroup
for(i in 1:nrow(events)){
  
  #LEADERSIHP BASED ON POSITION ALONG THE GROUP TRAJECTORY AT SPECIFIC TIMES
  #get the leader info for each subgroup for each event - definition based on order at a sepcific time
  lead_pos_A_start <- ind_disp_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$start_time[i])
  lead_pos_A_mid <- ind_disp_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], floor((events$start_time[i]+events$end_time[i])/2))
  lead_pos_A_end <- ind_disp_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$end_time[i])
  lead_pos_B_start <- ind_disp_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$start_time[i])
  lead_pos_B_mid <- ind_disp_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], floor((events$start_time[i]+events$end_time[i])/2))
  lead_pos_B_end <- ind_disp_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$end_time[i])
  
  #save displacement info into events dataframe
  events$group_A_lead_disp_start[i] <- list(lead_pos_A_start$disp)
  events$group_B_lead_disp_start[i] <- list(lead_pos_B_start$disp)
  events$group_A_lead_disp_mid[i] <- list(lead_pos_A_mid$disp)
  events$group_B_lead_disp_mid[i] <- list(lead_pos_B_mid$disp)
  events$group_A_lead_disp_end[i] <- list(lead_pos_A_end$disp)
  events$group_B_lead_disp_end[i] <- list(lead_pos_B_end$disp)
  events$group_A_lead_pos_norm_rank_start[i] <- list(lead_pos_A_start$norm_ranks)
  events$group_B_lead_pos_norm_rank_start[i] <- list(lead_pos_B_start$norm_ranks)
  events$group_A_lead_pos_norm_rank_mid[i] <- list(lead_pos_A_mid$norm_ranks)
  events$group_B_lead_pos_norm_rank_mid[i] <- list(lead_pos_B_mid$norm_ranks)
  events$group_A_lead_pos_norm_rank_end[i] <- list(lead_pos_A_end$norm_ranks)
  events$group_B_lead_pos_norm_rank_end[i] <- list(lead_pos_B_end$norm_ranks)
  
  #LEADERSHIP BASED ON TIME OF CROSSING A DEFINED THRESHOLD DISTANCE ALONG GROUP TRAJECTORY
  #compute 
  lead_crosstime_A <- ind_crossing_thresh_times_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], dist_frac_thresh = dist_frac_thresh)
  lead_crosstime_B <- ind_crossing_thresh_times_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], dist_frac_thresh = dist_frac_thresh)
  
  #store in the data frame
  events$group_A_lead_crosstime_norm_rank[i] <- list(lead_crosstime_A$norm_ranks)
  events$group_B_lead_crosstime_norm_rank[i] <- list(lead_crosstime_B$norm_ranks)
  
  #LEADERSHIP BASED ON TIME OF CROSSING YOUR OWN FINISH LINE ALONG GROUP TRAJECTORY
  #compute 
  lead_crosstime_ownfinishline_A <- ind_crossing_thresh_times_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], dist_frac_thresh = dist_frac_thresh, own_finish_line = own_finish_line)
  lead_crosstime_ownfinishline_B <- ind_crossing_thresh_times_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], dist_frac_thresh = dist_frac_thresh, own_finish_line = own_finish_line)
  
  #store in the data frame
  events$group_A_lead_crosstime_ownfinishline_norm_rank[i] <- list(lead_crosstime_ownfinishline_A$norm_ranks)
  events$group_B_lead_crosstime_ownfinishline_norm_rank[i] <- list(lead_crosstime_ownfinishline_B$norm_ranks)
  
}

#Compute entropies for real data vs permuted data

#real data
out <- get_fission_fusion_leadership(events, n_inds, leadership_type = 'crosstime_ownfinishline', meas_time = 'end', norm_rank_bins = norm_rank_bins)
fission_entropies_data <- out$fission_leadership_entropies
fusion_entropies_data <- out$fusion_leadership_entropies
fission_leaders <- out$fission_leaders
fusion_leaders <- out$fusion_leaders

#permuted data - swap identities within subgroups
fusion_entropies_rand <- fission_entropies_rand <- matrix(NA, nrow = n_inds, ncol = n_rands)
for(n in 1:n_rands){
  events_rand <- events
  for(i in 1:nrow(events)){
    #get group A and gorup B indexes from the original data
    group_A_idxs <- events$group_A_idxs[i][[1]]
    group_B_idxs <- events$group_B_idxs[i][[1]]
    
    #randomize the order within group A and gruop B separately
    group_A_idxs_rand <- sample(group_A_idxs)
    group_B_idxs_rand <- sample(group_B_idxs)
    
    #store in the new events dataframe
    events_rand$group_A_idxs[i] <- list(group_A_idxs_rand)
    events_rand$group_B_idxs[i] <- list(group_B_idxs_rand)
    
  }
  
  out <- get_fission_fusion_leadership(events_rand, n_inds, leadership_type = 'crosstime_ownfinishline', meas_time = 'end', norm_rank_bins = norm_rank_bins)
  fission_entropies_rand[,n] <- out$fission_leadership_entropies
  fusion_entropies_rand[,n] <- out$fusion_leadership_entropies
  
}

fission_leaders_rand <- out$fission_leaders
fusion_leaders_rand <- out$fusion_leaders


#----------CORRELATIONS ACROSS LEADERSHIP METRICS---------
#position
out <- get_fission_fusion_leadership(events, n_inds, leadership_type = 'position', meas_time = 'mid', norm_rank_bins = norm_rank_bins)
fission_leaders_pos <- out$fission_leaders
fusion_leaders_pos <- out$fusion_leaders

filename <- paste0(groupdir, group, "_LeaderRank_position.RData")
save(out, file = filename)

#crosstime
out <- get_fission_fusion_leadership(events, n_inds, leadership_type = 'crosstime', meas_time = 'mid', norm_rank_bins = norm_rank_bins)
fission_leaders_crosstime <- out$fission_leaders
fusion_leaders_crosstime <- out$fusion_leaders

filename <- paste0(groupdir, group, "_LeaderRank_crosstime.RData")
save(out, file = filename)


out <- get_fission_fusion_leadership(events, n_inds, leadership_type = 'crosstime_ownfinishline', meas_time = 'mid', norm_rank_bins = norm_rank_bins)
fission_leaders_crosstime_ownfinishline <- out$fission_leaders
fusion_leaders_crosstime_ownfinishline <- out$fusion_leaders

filename <- paste0(groupdir, group, "_LeaderRank_crosstime_ownfinishline.RData")
save(out, file = filename)


#----------PLOTTING-------

if(make_plots){
  
  #Make a histogram of the distances moved by the groups during splits and merges
  #and plot the distance moved threshold on top of it as a sanity check
  quartz()
  hist(c(events$A_during_disp, events$B_during_disp), breaks = 20)
  abline(v = dist_moved_thresh, col = 'red', lwd = 3)
  
  #plot the % of moving out of all fission and fusion events for each individual w/ confidence intervals
  #fissions
  quartz()
  plot(NULL, xlim = c(0,100), ylim = c(0,n_inds), xlab = '% moved', yaxt = 'n', ylab = '', main = 'Fissions')
  arrows(fission_CIs[,1]*100,1:n_inds, fission_CIs[,2]*100, 1:n_inds, length = 0.1, code = 3, angle = 90, lwd = 2)
  points(fission_move_fracs*100, 1:n_inds, pch = 19, cex = 2)
  axis(2, at = 1:n_inds, labels = coati_ids$name, las =1)
  
  #fusions
  quartz()
  plot(NULL, xlim = c(0,100), ylim = c(0,n_inds), xlab = '% moved', yaxt = 'n', ylab = '', main = 'Fusions')
  arrows(fusion_CIs[,1]*100,1:n_inds, fusion_CIs[,2]*100, 1:n_inds, length = 0.1, code = 3, angle = 90, lwd = 2)
  points(fusion_move_fracs*100, 1:n_inds, pch = 19, cex = 2)
  axis(2, at = 1:n_inds, labels = coati_ids$name, las =1)
  
  #fissions vs fusions
  quartz()
  plot(fusion_move_fracs, fission_move_fracs)
  
  #leadership during fusions
  quartz()
  par(mfrow=c(5,5), mar = c(2,2,0,0))
  for(i in 1:n_inds){
    hist(fusion_leaders[i,], breaks= seq(0,1,.2), main = coati_ids$name[i])
  }
  
  #leadership during fissions
  quartz()
  par(mfrow=c(5,5), mar = c(2,2,0,0))
  for(i in 1:n_inds){
    hist(fission_leaders[i,], breaks= seq(0,1,.2), main = coati_ids$name[i])
  }
  
  
  #compare entropy between real and permuted data
  #test statistic = mean entropy
  quartz()
  fission_means_rand <- colMeans(fission_entropies_rand,na.rm=T)
  fission_mean_data <- mean(fission_entropies_data,na.rm=T)
  hist(fission_means_rand, breaks=20,main = 'Fission', xlab = 'Mean entropy')
  abline(v=fission_mean_data,col='red',lwd=2)
  p_fission <- sum((fission_means_rand < fission_mean_data),na.rm=T)/n_rands
  
  quartz()
  fusion_means_rand <- colMeans(fusion_entropies_rand, na.rm=T)
  fusion_mean_data <- mean(fusion_entropies_data, na.rm=T)
  hist(fusion_means_rand, breaks=20, main = 'Fusion', xlab = 'Mean entropy')
  abline(v=fusion_mean_data,col='red',lwd=2)
  p_fusion <- sum((fusion_means_rand < fusion_mean_data),na.rm=T)/n_rands
  
  
}


