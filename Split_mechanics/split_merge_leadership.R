#Script to look at leadership in split and merge events
#Q1: Who leaves and who stays
#Q2: Who leads the movement of the group when leaving or joining

#-----PARAMETERS-------

user <- 'ari'
group <- 'presedente'
use_manual_events <- F
dist_moved_thresh <- 15 #minimum distance moved by a subgroup to count it as having moved (i.e. left or joined)

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
  proj_i <- (dxi*dxc + dyi*dyc) / sqrt(dxc^2 + dyc^2)
  
  #return ind displacement along group vector
  return(proj_i) 
}

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


#Make a histogram of the distances moved by the groups during splits and merges
#and plot the distance moved threshold on top of it as a sanity check
quartz()
hist(c(events$A_during_disp, events$B_during_disp), breaks = 20)
abline(v = dist_moved_thresh, col = 'red', lwd = 3)

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

#LEADERSHIP DURING SPLITS

#get the displacement of each individual projectd along the subgroup vector, for each subgroup
events$group_A_lead_disp <- list(c(0,0,0))
events$group_B_lead_disp <- list(c(0,0,0))
events$group_A_lead_rank <- list(c(0,0,0))
events$group_B_lead_rank <- list(c(0,0,0))
for(i in 1:nrow(events)){
  events$group_A_lead_disp[i] <- list(ind_disp_along_group_path(events$group_A_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$end_time[i]))
  events$group_B_lead_disp[i] <- list(ind_disp_along_group_path(events$group_B_idxs[i][[1]], xs, ys, events$start_time[i], events$end_time[i], events$end_time[i]))
}

#get ranks (lower values = behind, higher values = in the lead)
for(i in 1:nrow(events)){
  
  #get ranks for group A and B
  ranks_A <- rank(events$group_A_lead_disp[i][[1]])
  ranks_B <- rank(events$group_B_lead_disp[i][[1]])
  
  #store raw ranks
  events$group_A_lead_rank[i] <- list(ranks_A)
  events$group_B_lead_rank[i] <- list(ranks_B)
  
  #get normalized ranks
  min_rank_A <- min(ranks_A)
  min_rank_B <- min(ranks_B)
  max_rank_A <- max(ranks_A)
  max_rank_B <- max(ranks_B)
  norm_ranks_A <- (ranks_A - min_rank_A) / (max_rank_A - min_rank_A)
  norm_ranks_B <- (ranks_B - min_rank_B) / (max_rank_B - min_rank_B)
  
  #store normalized ranks
  events$group_A_lead_rank_norm[i] <- list(norm_ranks_A)
  events$group_B_lead_rank_norm[i] <- list(norm_ranks_B)
  
}

#get normalized ranks for individuals during fissions and fusions
fission_leaders <- fusion_leaders <- matrix(NA, nrow = n_inds, ncol = nrow(events))
for(i in 1:nrow(events)){
  
  if(!is.na(events$A_moved[i])){
    if(events$A_moved[i]){
      inds <- events$group_A_idxs[i][[1]]
      norm_ranks <- events$group_A_lead_rank_norm[i][[1]]
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
      norm_ranks <- events$group_B_lead_rank_norm[i][[1]]
      if(events$event_type[i] == 'fission'){
        fission_leaders[inds,i] <- norm_ranks
      }
      if(events$event_type[i] == 'fusion'){
        fusion_leaders[inds,i] <- norm_ranks
      }
    }
  }
}

rowMeans(fission_leaders, na.rm=T)
apply(fission_leaders, 1, sd, na.rm=T)
rowMeans(fusion_leaders, na.rm=T)

quartz()
par(mfrow=c(5,5), mar = c(2,2,0,0))
for(i in 1:n_inds){
  hist(fusion_leaders[i,], breaks= seq(0,1,.2), main = coati_ids$name[i])
}

quartz()
par(mfrow=c(5,5), mar = c(2,2,0,0))
for(i in 1:n_inds){
  hist(fission_leaders[i,], breaks= seq(0,1,.2), main = coati_ids$name[i])
}
