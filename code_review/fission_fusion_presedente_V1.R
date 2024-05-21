#fission fusion analysis script with Presidente group
#first define groups for each time step (each 10 mins)

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level1.RData"
id_file <- 'presedente_coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 70
#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(vioplot)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)


#---------------------------------------------------------------
#calculate prop time males are with group

#n_inds <- nrow(xs)
#n_times <- ncol(xs)
# subgroup_data_males <- get_subgroup_data(xs, ys, 50)
# #males are 5,8,9,10,18,21 #ardern is 1, merkel is 13, Truss is 20 - used as measure for being in group
# subgroup_id <- subgroup_data_males$ind_subgroup_membership[c(1,5,8,9,10,13,18,20,21),]
# rownames(subgroup_id) <- coati_ids$name[c(1,5,8,9,10,13,18,20,21)]
# 
# #because I'm curious - number of times Wildflower is with males- only with Gendry
# sum(subgroup_id[9,] == subgroup_id[2,], na.rm = T) #32 #Wildflower and Gendry
# 
# #get sum of times the male is with either Ardern, Merkel, or Truss
# gendry_group <- sum((subgroup_id[1,] == subgroup_id[2,] | subgroup_id[6,] == subgroup_id[2,] | subgroup_id[8,] == subgroup_id[2,]), na.rm = TRUE) #Ardern, Merkel, Truss with Gendry #144
# kenyatta_group <- sum((subgroup_id[1,] == subgroup_id[3,] | subgroup_id[6,] == subgroup_id[3,] | subgroup_id[8,] == subgroup_id[3,]), na.rm = TRUE) #Ardern, Merkel, Truss with Kenyatta #253
# lula_group <- sum((subgroup_id[1,] == subgroup_id[4,] | subgroup_id[6,] == subgroup_id[4,] | subgroup_id[8,] == subgroup_id[4,]), na.rm = TRUE)#Ardern, Merkel, Truss with Lula #198
# mandela_group <- sum((subgroup_id[1,] == subgroup_id[5,] | subgroup_id[6,] == subgroup_id[5,] | subgroup_id[8,] == subgroup_id[5,]), na.rm = TRUE) #Ardern, Merkel, Truss with Mandela #144
# sam_group <- sum((subgroup_id[1,] == subgroup_id[7,] | subgroup_id[6,] == subgroup_id[7,] | subgroup_id[8,] == subgroup_id[7,]), na.rm = TRUE) #Ardern, Merkel with Sam #155
# 
# #sum gps points for each male
# gendry_sum <- sum(!is.na(subgroup_id[2,])) #Gendry 958
# kenyatta_sum <- sum(!is.na(subgroup_id[3,])) #Kenyatta 1046
# lula_sum <- sum(!is.na(subgroup_id[4,])) #Lula 1003
# mandela_sum <- sum(!is.na(subgroup_id[5,])) #Mandela 985
# sam_sum <- sum(!is.na(subgroup_id[7,])) #Sam 872
# 
# #prop time with group
# gendry_prop <- gendry_group/gendry_sum
# kenyatta_prop <- kenyatta_group/kenyatta_sum
# lula_prop <- lula_group/lula_sum
# mandela_prop <- mandela_group/mandela_sum
# sam_prop <- sam_group/sam_sum



#---------------------------------------------------------------


#removing males and wildflower
xs <- xs[c(1:4,6,7,11:17,19,20,22),]
ys <- ys[c(1:4,6,7,11:17,19,20,22),]
coati_ids <- coati_ids[-c(5,8,9,10,18,21),]

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

#------------------------------------------------------------------------------------

#Figure A1: number of sub groups when the radius is changed (graph put in dropbox results folder) 
png(height = 2200, width = 900, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists_onlyGroup_level1.png'))

par(mfrow=c(6,1), mar = c(8,9,12,1), mgp=c(6,1,0)) #(bottom, left, top, right)

for (i in 1:length(Rs)){
  
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  
  
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "aquamarine3", breaks = seq(.5,16,1), cex.lab = 4, cex.main = 4, cex.axis=3, freq = FALSE, ylim=c(0,.7), xlim=c(0.5,16), axes = F)
  axis(2, at = c(0,0.3,0.6), cex.axis = 4, las = 1, pos = 0.5)
  axis(1, at = c(5,10,15), cex.axis = 4, las = 1, pos = 0, padj = 0.6)
  
  
  # Add the title only for the first plot
  if(i == 1) {
    mtext("Presidente group", side = 3, line = 9, cex = 3)
  }
  
  
}

dev.off()

#-------------------------------------------------------------------------------------

R <- 70

# Figure 2d,e,f: Characterizing the subgroup patterns when group was split into 2 or 3 subgroups
png(height = 400, width = 1450, units = 'px', filename = paste0(plot_dir,'50m_charecterisations.png'))

blue <- rgb(69, 119, 116, maxColorValue = 255)

lightblue <- rgb(128, 218, 178, maxColorValue = 255)

par(mfrow=c(1,3), mar = c(6,11,2,1)) #(bottom, left, top, right)
par(mgp = c(4,0.4,-1.25)) #axis title, axis labels, axis line
subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs],main = "", xlab =  'Number of subgroups (radius = 70 m)', col = "aquamarine3", breaks = seq(.5,16,1), cex.lab = 2.9, cex.main = 3, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 6), border = lightblue, las = 1, yaxt = "n", xaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0:7), cex.axis = 3, las = 1, pos = 0, padj = 0.8)
text(x = 0.3, y = 0.55,label = "(d)", cex = 3)
mtext("Presidente group", side = 2, line = 8.5, cex = 2)

subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], breaks=seq(0.5,16,1), xlab = expression("Subgroup size (" * italic("N") * " subgroups = 2)"), main = '', col = blue, cex.lab = 3, cex.main = 2.9, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 16), border = "aquamarine4", las = 1,yaxt = "n", xaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,5,10,15), cex.axis = 3, las = 1, pos = 0, padj = 0.8)

text(x = 0.8, y = 0.55,label = "(e)", cex = 3)

hist(subgroup_counts[,s3], breaks=seq(0.5,16,1), xlab = expression("Subgroup size (" * italic("N") * " subgroups = 3)"), main = '', col = blue, cex.lab = 3, cex.main = 2.9, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 16), border = "aquamarine4", las = 1,yaxt = "n", xaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3, las = 1, pos = 0, hadj = 1.2)
axis(1, at = c(0,5,10,15), cex.axis = 3, las = 1, pos = 0, padj = 0.8)
text(x = 0.8, y = 0.55,label = "(f)", cex = 3)

dev.off()

#calculate proportion of times the group were either together or split into two, or three subgroups
all <- as.data.frame(table(subgroup_data$n_subgroups[all_tracked_idxs]))
sum_2.3 <- sum(all[c(2:3),2])
sum_4.5 <- sum(all[c(4:5),2])
sum_all <- sum(all[-1,2])
prop_2.3 <- (sum_2.3/ sum_all)*100




#-------------------------------------------------------------------------------------

#Figure 3c: which individuals tend to be in the same subgroup for only the group

subgroup_data <- get_subgroup_data(xs, ys, R=R)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA

#order - without Wildflower AND males: 
new_order <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14)

par(mgp=c(3, 1, 0), mar=c(11,10,4,2))

ffnet_reorder <- ff_net[new_order, new_order]
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_onlyGroup_level1.png'))

visualize_network_matrix_presedente(ffnet_reorder, coati_ids[new_order,])

dev.off()

#save matrix for mrqap analysis
write.table(ffnet_reorder,file="C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/presedente_matrix_10min_proptimeinsamesubgroup_50m.txt",row.names=FALSE)

#-------------------------------------------------------------------------------------

#Figure A3b,d and A6b,d: within full group individual associations to compare with the sub group memberships

#this is the probability of individuals being in the same sub-group using absolute dyadic distances 

subgroup_data <- get_subgroup_data(xs, ys, R)

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)
#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)
#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs[,full_group_index]
subset_y <- ys[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 10)

new_order <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14)

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'withingroup_network_onlyGroup_level1.png'))
visualize_network_matrix_presedente(within_group_data$proximity_net, coati_ids[new_order,])
dev.off()


#-------------------------------------------------------------------

#Plotting the mean sub-group size for each hour for all individuals

#get mean group size for each timepoint - so need to get the total number of individuals divided by the number of groups

# sub_counts <- as.data.frame(t(subgroup_data$subgroup_counts))
# colnames(sub_counts) <- c("sub1", "sub2", "sub3", "sub4")
# n_groups <- as.data.frame(subgroup_data$n_subgroups)
# #combine to make dataframe with the number of individuals in each sub group and the number of subgroups
# n_subs <- cbind(sub_counts, n_groups)
# colnames(n_subs)[colnames(n_subs) == 'subgroup_data$n_subgroups'] <- 'n_groups'
# #adding time to n_subs
# n_subs <- cbind(n_subs, ts)
# 
# #get the hour
# n_subs$hour <- hour(n_subs$ts)
# 
# #get total number of individuals
# n_subs$n_inds <- n_tracked
# 
# #calculate the mean group size
# n_subs$mean_group_size <- NA
# 
# n_subs$mean_group_size <- (n_subs$n_inds/n_subs$n_groups)
# 
# #change hour to Panama time
# n_subs$panama_time <- n_subs$hour-5
# n_subs$date <- as.Date(n_subs$ts)
# #n_subs_pres <- n_subs
# 
# #save n_subs as an RData object
# #save(n_subs_pres, file = "C:/Users/egrout/Dropbox/coatithon/processed/n_subs_pres.Rdata")
# #this dataframe will be used in the subgoups_vioplot_combined script to compare with the other group
# 
# #now plotting mean group size for each hour of the day
# png(height = 500, width = 700, units = 'px', filename = paste0(plot_dir, "mean_group_size_violin_level1_green.png"))
# vioplot(n_subs$mean_group_size ~ n_subs$panama_time,  xlab = "panama time", ylab = "mean subgroup size",col = "aquamarine3", cex.axis = 1.5, colMed = "black")
# dev.off()

#-------------------------------------------------------------------

#Figure A9c,d: make plot for the number of GPS points recorded for each individual

png(height = 600, width = 1400, units = 'px', filename = paste0(plot_dir, "number_tracked.png"))
par(mfrow=c(1,2), mar = c(10,9,2,1)) #c(bottom, left, top, right)
sum_tracked <- colSums(!is.na(xs))
hist(sum_tracked[ !(sum_tracked==0) ], ylim = c(0, 400), main = "", xlab = "Number of individuals tracked", ylab = "Frequency", col = "aquamarine4", cex.lab = 2, cex.axis = 2, yaxt = "n", xaxt = "n")
axis(2, at = c(0,100,200,300,400), cex.axis = 2, las = 1, pos = 1)
axis(1, cex.axis = 2, las = 1, pos = 0, padj = 0.6)
text(x = 2, y = 400,label = "(c)", cex = 3)
mtext("Presidente group", side = 2, line = 7, cex = 2)


each_sum <- data.frame(sum = rowSums(!is.na(xs)))
barplot(each_sum$sum, names.arg = coati_ids$name, las=2, col = "aquamarine3", ylab = "Number of GPS points",  cex.lab = 2, cex.axis = 2, cex.names=2, mgp=c(5,1,0), ylim = c(0, 1100))
text(x = 0.5, y = 1050,label = "(d)", cex = 3)

dev.off()

#------------------------------------------------------------------------------

#calculate the proportion of missing data
max(each_sum$sum)
min(each_sum$sum)
each_sum$missing <- max(each_sum$sum)- each_sum$sum
each_sum$prop <- (each_sum$missing/max(each_sum$sum))*100
mean(each_sum$prop)
sd(each_sum$prop)



#---------------------------------------------------------------------

#Figure A4,7,12: histogram for consistency 

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, R)

splits <- c()
for(t in 1:(n_times-1)){
  
  #get subgroup membership now and later
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(subgroups_now))==0 | sum(!is.na(subgroups_later))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
  n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(subgroups_now)==1)
  singletons_later <- sum(table(subgroups_later)==1)
  
  #determine if this time step is a split
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if((n_subgroups_now==1 )
     & n_subgroups_later >n_subgroups_now
     & singletons_later==0
  ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & ((singletons_now+1) == n_subgroups_now )
     & n_subgroups_later > n_subgroups_now
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
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(subgroups_now)
  orig_subgroup_members <- which(subgroups_now == orig_subgroup_number)
  
  #store original group membership in data frame
  splits_df$orig_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_later[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_later==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
    #put lists into a data frame
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

#this dataframe is used in merge_analysis_presedente_V1
save(splits_df, file = "C:/Users/egrout/Dropbox/coatithon_notgithub/Presedente_fission_fusion/splits_df.Rdata")  

#DONE WITH DATAFRAME!

p_dyad_together <- get_p_dyad_together(splits_df_local = splits_df, n_inds_local = n_inds)

#Compute consistency based on consistency metric (see coati_function_library)
consistency_data <- get_consistency(p_dyad_together)

#randomize splits and recompute consistency metric
n_rands <- 1000
set.seed(24)
rando_consistencies <- rep(NA, n_rands)
for (i in 1:n_rands){
  
  rando <- randomise_splits(splits_df)
  rando_dyad <- get_p_dyad_together(splits_df_local = rando, n_inds_local = n_inds)
  consist <- get_consistency(rando_dyad)
  rando_consistencies[i] <- consist
  
}

png(height = 600, width = 900, units = 'px', filename = paste0(plot_dir,'consistency_splits_hist.png'))
par(mfrow=c(1,1), mar = c(6,6,5,3))#(bottom, left, top, right)
hist(rando_consistencies, breaks=seq(0,0.5,0.005), main = "", xlab = "Consistency of random sub-group allocations", xlim = c(0,0.4), col = "slategray4",  cex.lab = 2.5, cex.axis=2.5, yaxt = "n", xaxt = "n") #for 70m, change xlim to 0.5
#abline(v=consistency_data, ylim = c(0, 50), col = 'orange2', lwd=4)
segments(x0 = consistency_data, x1 = consistency_data, y0 = 0.01, y1 = 300, col = 'orange2', lwd=4)

axis(2, at = c(0,100,200,300), cex.axis = 2.5, las = 1, pos = 0, padj = 0.2)
axis(1, at = c(0.1,0.2,0.3, 0.4), cex.axis = 2.5, las = 1, pos = 0, padj = 0.3) #for 70m, change to 0.5
text(x = 0.017, y = 250,label = "(b)", cex = 3) #50m y = 250, x = 0.017, ---- 30m y = 228, ---- 70m, y = 226 and x = 0.0.022
mtext("Presidente group", side = 3, line = 2, cex = 2.5)

dev.off()


#this graph shows that the mean consistency of individuals splitting with others is more likely than randomly assigning individuals into sub-groups

#should look into repeatability of binary data - to see if there is a better metric for getting consistency values from binary data

#----------------------------------------------------------------

#consistency matrix

#symmetrize matrix - copy entries from p_dyad_together[i,j] to p_dyad_together[j,i]
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    p_dyad_together[j,i] <- p_dyad_together[i,j]
  }
}

diag(p_dyad_together) <- NA

new_order <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14)

p_dyad_together_reorder <- p_dyad_together[new_order, new_order]

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_splits.png'))
par(mfrow=c(1,1), mar = c(1,2,1,1))#(bottom, left, top, right)
visualize_network_matrix_presedente(p_dyad_together_reorder, coati_ids[new_order,])
dev.off()

#--------------------------------------------------------------------------------------

#Figure A8: look at which age/sex classes tend to be on their own
inds_subgroup <- data.frame(subgroup_data$ind_subgroup_membership)
#make an empty dataframe to add the alone inds data to
df <- data.frame(matrix(nrow = 16, ncol = ncol(inds_subgroup)))

# Loop through columns
for (i in 1:ncol(inds_subgroup)) {
  # Find unique rows without considering NAs and excluding rows with all NAs
  unique_rows <- which(!duplicated(inds_subgroup[, i], na.rm = TRUE) & 
                         !duplicated(inds_subgroup[, i], fromLast = TRUE, na.rm = TRUE) &
                         !is.na(inds_subgroup[, i]))
  # Assign 1 to the cells where an individual is alone
  df[unique_rows, i] <- 1
}


coati_ids$inds_alone <- rowSums(df, na.rm = TRUE)
#get sum of times each ind has data
coati_ids$total_gps <- ncol(inds_subgroup) - rowSums(is.na(inds_subgroup))
#get proportion of time alone
coati_ids$prop_alone <- coati_ids$inds_alone/coati_ids$total_gps
coati_ids$prop_alone <- coati_ids$prop_alone*100

coati_ids$age_sex <- paste(coati_ids$age, coati_ids$sex, sep = " ")

colors <- c("orange3","orange2","orange","orange", "aquamarine4", "aquamarine3")

gg <- ggplot(aes(x = prop_alone, y = age_sex), data = coati_ids)+
  xlab("Proportion of time alone (%)")+
  ylab("Age class")+
  geom_point(aes(color = interaction(age, sex)), position = position_identity(), size = 3)+
  facet_grid(vars(age_sex), scales = "free", space = "free")+
  scale_color_manual(values = colors) +
  theme_classic()+
  guides(color = "none")+  # Remove the legend
  theme( strip.text.y = element_blank())

#change Ardera to Ardern (spelling error)

coati_ids$name[coati_ids$name == "Ardera"] <- "Ardern"

gg <- ggplot(aes(x = prop_alone, y = age_sex, color = name), data = coati_ids) +
  xlab("Proportion of time alone (%)") +
  ylab("Age class") + 
  scale_x_continuous(limits = c(0, 40))+
  geom_point(aes(color = name), position = position_jitter(width = 0, height = 0.3), size = 3, alpha = 0.7) +
  facet_grid(vars(age_sex), scales = "free", space = "free") +
  scale_color_manual(values = c("Ardern" = "#a6cee3", "Castro" = "#1f78b4", "Cleopatra" = "#b2df8a", "Gandhi" = "#33a02c", "Khan" = "#fb9a99", "Gillard" = "#e31a1c", "May" = "#fdbf6f", "Meir" = "#ff7f00", "Merkel" = "#cab2d6", "Moscoso" = "#6a3d9a", "Mujica" = "#b15928", "Obama"= "#EEEE00", "Peron" = "#feb24c", "Torrijos" = "#CD1076", "Truss" = "#bd0026" , "Zelenskyy" = "#006400")) + 
  theme_classic() +
  guides(color = guide_legend(title = "Individual"))+  # Add legend title
  theme( strip.text.y = element_blank())+
  NULL

gg


coati_ids_alone_pres <- coati_ids
#this dataframe is used in alone_inds_all_groups
#save(coati_ids_alone_pres, file = "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/pres_alone_inds_level1.Rdata")  

ggsave(filename = paste0(plot_dir, 'prop_time_alone.png'), plot = gg, width = 4, height = 5, dpi = 300)


