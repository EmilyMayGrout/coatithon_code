
#fission fusion analysis script for Galaxy group 

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
gps_file <- "galaxy_xy_10min_level1.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 50

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vioplot)
library(plotly)

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

#------------------------------------------------------------------------------
#find prop time Gus with the group

#n_inds <- nrow(xs)
#n_times <- ncol(xs)
# subgroup_data_males <- get_subgroup_data(xs, ys, 50)
# #group will be orbita is 6 and cometa is 11
# #gus is 5 
# subgroup_id <- subgroup_data_males$ind_subgroup_membership[c(5,6,11),]
# rownames(subgroup_id) <- coati_ids$name[c(5,6,11)]
# 
# #prop time the 2 subgroups are together (Orbita and Cometa) = 605
# sum((subgroup_id[1,] == subgroup_id[3,]), na.rm = TRUE) 
# #prop time Gus with either subgroup = 780
# gus_group <- sum((subgroup_id[1,] == subgroup_id[2,] | subgroup_id[3,] == subgroup_id[2,]), na.rm = TRUE)
# 
# #number of gps points for Gus
# gus_sum <- sum(!is.na(subgroup_id[2,]))
# gus_prop <- gus_group/gus_sum

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#----------------------------------------------------------------------

# Figure 2a,b,c: Characterizing the subgroup patterns when group was split into 2 or 3 subgroups
png(height = 400, width = 1280, units = 'px', filename = paste0(plot_dir,'50m_charecterisations.png'))

par(mfrow=c(1,3), mar = c(6,7,2,1)) #(bottom, left, top, right)
par(mgp = c(4,0.4,-1.2)) #axis title, axis labels, axis line
#for dark olivegreen3 match
green <- rgb(120, 170, 80, maxColorValue = 255)

darkgreen <- rgb(120, 160, 80, maxColorValue = 255)

subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs], main = "", xlab =  'Number of subgroups (radius = 50m)', col = green, breaks = seq(.5,11,1), cex.lab = 2.9, cex.main = 3, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 6), border = "darkolivegreen3", las = 1, yaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3, las = 1)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), main = "", xlab = 'Subgroup size (N subgroups = 2)',  col = "darkolivegreen", cex.lab = 2.9, cex.main = 3, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), border = darkgreen, las = 1, yaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3,las = 1)
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), main = "", xlab = 'Subgroup size (N subgroups = 3)',  col = "darkolivegreen", cex.lab = 2.9, cex.main = 3, cex.axis=3, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), border = darkgreen, las = 1, yaxt = "n")
axis(2, at = c(0,0.3,0.6), cex.axis = 3, las = 1)

dev.off()

#-----------------------------------------------------------

#Figure S1a: Number of sub groups when the radius is changed (graph put in dropbox results folder)

png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists_level1.png'))
par(mfrow=c(6,1), mar = c(8,7,1,1), mgp=c(4,1,0))

for (i in 1:length(Rs)){
  
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2.5, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7), yaxt = "n")
  axis(2, at = c(0,0.3,0.6), cex.axis = 2, las = 1)
  
}
dev.off()

#--------------------------------------------------------------
R <- 50

#Figure 3a: which individuals tend to be in the same subgroup

subgroup_data <- get_subgroup_data(xs, ys, R= R)

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
new_order <- c(5,1,11,4,10,2,3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]
#save matrix for mrqap analysis
write.table(ffnet_reorder,file="C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/gal_matrix_10min_proptimeinsamesubgroup_50m.txt",row.names=FALSE)


png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_level1_cut.png'))

visualize_network_matrix_galaxy(ffnet_reorder, coati_ids[new_order,])
dev.off()

#--------------------------------------------------------------------------

#Probability of individuals being in the same sub-group using absolute dyadic distances 
#removing Gus (adult male) from the full group to get proximity data within group

xs_nogus <- xs[-c(5), ]
ys_nogus <- ys[-c(5), ]
coati_ids_nogus <- coati_ids[-c(5),]

subgroup_data <- get_subgroup_data(xs_nogus, ys_nogus, R)

n_inds <- nrow(xs_nogus)
n_times <- ncol(xs_nogus)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs_nogus))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)
#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)
#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs_nogus[,full_group_index]
subset_y <- ys_nogus[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 10)

new_order <- c(1,10,4,9,2,3,5,6,7,8)

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'withingroup_network_withoutgus_level1.png'))
visualize_network_matrix_galaxy(within_group_data$proximity_net, coati_ids_nogus[new_order,])
dev.off()

#there's an additional 42 data points but the proximity values don't change much
#this gives more data points by not including Gus


#-------------------------------------------------------------------

#Plotting the mean sub-group size for each hour for all individuals

#get mean group size for each time point - so need to get the total number of individuals divided by the number of groups

# sub_counts <- as.data.frame(t(subgroup_data$subgroup_counts))
# colnames(sub_counts) <- c("sub1", "sub2", "sub3", "sub4", "sub5")
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
# #n_subs_gal <- n_subs
# 
# #save n_subs as an RData object
# #save(n_subs_gal, file = "C:/Users/egrout/Dropbox/coatithon/processed/n_subs_gal.Rdata")
# #this dataframe will be used in the subgoups_vioplot_combined script to compare with the other group
# 
# #now plotting mean group size for each hour of the day
# png(height = 500, width = 700, units = 'px', filename = paste0(plot_dir, "mean_group_size_violin_level1_red.png"))
# vioplot(n_subs$mean_group_size ~ n_subs$panama_time,  cex.axis = 1.5, cex.lab = 1, xlab = "panama time", ylab = "mean subgroup size", col = "darkolivegreen3", colMed = "black")
# dev.off()


#---------------------------------------------------------------------

#Figure S2ab: make plot for the number of GPS points recorded for each individual

png(height = 600, width = 1400, units = 'px', filename = paste0(plot_dir, "number_tracked.png"))
par(mfrow=c(1,2), mar = c(8,8,2,1)) #c(bottom, left, top, right)
sum_tracked <- colSums(!is.na(xs))
hist(sum_tracked[ !(sum_tracked==0) ], ylim = c(0, 800), main = "", xlab = "Number of individuals tracked", ylab = "Frequency", col = "darkolivegreen4", cex.lab = 2, cex.axis = 2)
each_sum <- data.frame(sum = rowSums(!is.na(xs)))
barplot(each_sum$sum, names.arg = coati_ids$name, las=2, col = "darkolivegreen3", ylab = "Number of GPS points",  cex.lab = 2, cex.axis = 2, cex.names=2, mgp=c(5,1,0))
dev.off()

#calculate the proportion of missing data
max(each_sum$sum)
min(each_sum$sum)
each_sum$missing <- max(each_sum$sum)- each_sum$sum
each_sum$prop <- (each_sum$missing/max(each_sum$sum))*100
mean(each_sum$prop)
sd(each_sum$prop)

#---------------------------------------------------------------------

#Figure S5a: histogram for consistency 

#need to run this again as Gus was removed for figure 6a
n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is R m
subgroup_data <- get_subgroup_data(xs, ys, R)


#run through time by time and check if each time is a split
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
  if((n_subgroups_now == 1 )
     & n_subgroups_later > n_subgroups_now
     & singletons_later==0
  ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & ((singletons_now + 1) == n_subgroups_now )
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
  
  #transfer NAs from now to later and later to now so errenous splits are not detected
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

#save dataframe as its needed for the figure1_fissionfusion_plot.R
save(splits_df, file = "C:/Users/egrout/Dropbox/coatithon_notgithub/Galaxy_fission_fusion/splits_df.Rdata")  

#DONE WITH DATAFRAME

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


#----------------------------------------------------------------

#Figure S8bd + S11bd: consistency matrix

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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_splits.png'))
par(mfrow=c(1,1), mar = c(1,2,1,1))#(bottom, left, top, right)
visualize_network_matrix_galaxy(p_dyad_together_reorder, coati_ids[new_order,])
dev.off()

#---------------------------------------------------------------------

#look at which age/sex classes tend to be on their own
inds_subgroup <- data.frame(subgroup_data$ind_subgroup_membership)
#make an empty dataframe to add the alone inds data to
df <- data.frame(matrix(nrow = 11, ncol = ncol(inds_subgroup)))

for(i in 1:ncol(inds_subgroup)){
  
  #find the individuals in each column who are on their own
  col_i <- inds_subgroup[,i]
  #get the indices of the rows which are on their own
  unique_rows <- which(!duplicated(inds_subgroup[,i]) & !duplicated(inds_subgroup[, i], fromLast = TRUE))
  #for the individual who is on its own, add 1 to that row in column i  
  df[unique_rows, i] <- 1
}

coati_ids$inds_alone <- rowSums(df, na.rm = TRUE)
#get sum of times each ind has data
coati_ids$total_gps <- ncol(inds_subgroup) - rowSums(is.na(inds_subgroup))
#get proportion of time alone
coati_ids$prop_alone <- coati_ids$inds_alone/coati_ids$total_gps

coati_ids$age_sex <- paste(coati_ids$age, coati_ids$sex, sep = " ")

colors <- c("orange3","orange2","orange","aquamarine4", "aquamarine3")

gg <- ggplot(aes(x = prop_alone, y = age_sex), data = coati_ids)+
  xlab("Proportion of time alone (%)")+
  ylab("Age class")+
  geom_point(aes(color = interaction(age, sex)), position = position_identity(), size = 3)+
  facet_grid(vars(age_sex), scales = "free", space = "free")+
  scale_color_manual(values = colors) +
  theme_classic()+
  guides(color = "none")+  # Remove the legend
  theme( strip.text.y = element_blank())

gg

coati_ids_alone_gal <- coati_ids
#this dataframe is used in alone_inds_all_groups
save(coati_ids_alone_gal, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/gal_alone_inds_level1.Rdata")  


ggsave(filename = paste0(plot_dir, 'prop_time_alone.png'), plot = gg, width = 6, height = 6, dpi = 300)

#-------------------------------------------------------------------------------------
#Figure 1e

subgroup_data_wide <- as.data.frame(rbind(subgroup_data$ind_subgroup_membership[,1:1633]))
#put id into dataframe
subgroup_data_wide <- cbind(coati_ids$name, subgroup_data_wide)
#changing column names into a list from 1 to 1633 - one day would be good to have as time
column_names <- c(1:1633)
colnames(subgroup_data_wide)[-1] <- column_names
#change ID's into numbers
subgroup_data_wide[subgroup_data_wide == 'Quasar'] <- '1'
subgroup_data_wide[subgroup_data_wide == 'Estrella'] <- '2'
subgroup_data_wide[subgroup_data_wide == 'Venus'] <- '3'
subgroup_data_wide[subgroup_data_wide == 'Lucero'] <- '4'
subgroup_data_wide[subgroup_data_wide == 'Gus'] <- '5'
subgroup_data_wide[subgroup_data_wide == 'Orbita'] <- '6'
subgroup_data_wide[subgroup_data_wide == 'Planeta'] <- '7'
subgroup_data_wide[subgroup_data_wide == 'Saturno'] <- '8'
subgroup_data_wide[subgroup_data_wide == 'Pluto'] <- '9'
subgroup_data_wide[subgroup_data_wide == 'Luna'] <- '10'
subgroup_data_wide[subgroup_data_wide == 'Cometa'] <- '11'

#pivot the data frame into a long format
subgroup_data_long <- subgroup_data_wide %>% pivot_longer(-c(`coati_ids$name`), names_to='time', values_to='sub-group')
subgroup_data_long$time <- as.numeric(subgroup_data_long$time)

#time df
time_df <- data.frame(ts, 1:length(ts))
names(time_df)[2] <- "time"
subgroup_data_long <- left_join(subgroup_data_long, time_df)

#filter test df to 05.01.22
start_time <- as.POSIXct("2022-01-05 11:00:00")
end_time <- as.POSIXct("2022-01-05 23:00:00")

# Subset the dataframe to include only rows within the date range
subset_df <- subgroup_data_long[subgroup_data_long$ts >= start_time & subgroup_data_long$ts <= end_time, ]

#convert time to Panama time
subset_df$ts <- subset_df$ts - hours(5)
colnames(subset_df) <- c("coati_ids.name", "time", "sub.group", "ts")

#prep the main subgrouping data for plotting
df_mod <- subset_df %>%
  fill(sub.group, .direction = "downup") %>% #get rid of NAs
  mutate(new_sg = sub.group) %>% #make a new subgroup column in case we mess with it (unnecessary step)
  arrange(ts, coati_ids.name, sub.group) %>% #sort by time, coati id, and subgroup
  group_by(ts, sub.group) %>% #group by time and subgroup
  #make new columns...
  mutate(n_in_sg = n(), #number of coatis in each subgroup 
         seq_in_sg = seq_along(coati_ids.name), #number each coati in each subgroup
         sg_jitter = new_sg - (0.07/2*n_in_sg) + (seq_in_sg*0.07)) # this is a trixie but 

df_mod2 <- df_mod #lets duplicated it before we mess with it, just in case we eff it up

#the following three little code chunks replace the df_mod2$new_sg for specific rows with 0
df_mod2$new_sg[which(df_mod2$coati_ids.name >=6 & 
                       df_mod2$coati_ids.name <= 9 &
                       df_mod2$ts == as.POSIXct("2022-01-05 12:20:00"))] <- 0

df_mod2$new_sg[which(df_mod2$coati_ids.name == 11 & 
                       (df_mod2$ts == as.POSIXct("2022-01-05 20:20:00") |
                          df_mod2$ts == as.POSIXct("2022-01-05 21:00:00") ))] <- 0

df_mod2$new_sg[which(df_mod2$sub.group != 1 & 
                       (df_mod2$ts >= as.POSIXct("2022-01-05 13:50:00") &
                          df_mod2$ts <= as.POSIXct("2022-01-05 15:40:00") ))] <- 0

df_mod2$new_sg[which(df_mod2$sub.group == 3 & 
                       (df_mod2$ts >= as.POSIXct("2022-01-06 09:10:00") &
                          df_mod2$ts <= as.POSIXct("2022-01-06 09:20:00") ))] <- 2

df_mod2$new_sg[which(df_mod2$sub.group == 2 & 
                       (df_mod2$ts >= as.POSIXct("2022-01-06 09:10:00") &
                          df_mod2$ts <= as.POSIXct("2022-01-06 09:20:00") ))] <- 3

df_mod2 <- df_mod2 %>%
  arrange(ts, coati_ids.name, sub.group) %>%
  group_by(ts, sub.group) %>%
  mutate(n_in_sg = n(),
         seq_in_sg = seq_along(coati_ids.name),
         sg_jitter = new_sg - (0.07/2*n_in_sg) + (seq_in_sg*0.07))

#make rectangles for each subgroup
df_group_polys <- df_mod2 %>%
  group_by(ts, sub.group) %>% # group by timestamp and subgroup
  summarize(top_subgroup = max(sg_jitter) + 0.1, #take the highest jitter value and add a little buffer
            bottom_subgroup = min(sg_jitter) - 0.1, #take the lowest jitter value and add a little buffer
            left_subgroup = mean(ts) - 180, #take that timestamp and subtract a little buffer
            right_subgroup = mean(ts) + 180) # take that time stamp and add a little buffer
#so, above, if you want more white space between the subgroup outlines and the points, use bigger buffers

#let's make some little lines that will connect subgroup boxes at the same timestamp (otherwise it's hard to see how they line up)
df_ts_lines <- df_mod2 %>%
  group_by(ts) %>%
  summarize(top_subgroup= max(sg_jitter), #take the highest point at each timestamp
            bottom_subgroup = min(sg_jitter)) # take the lowest point at each timestamp

p1 <- ggplot() +
  geom_rect(data = df_group_polys,
            mapping = aes(xmin = left_subgroup,
                          xmax = right_subgroup,
                          ymin = bottom_subgroup,
                          ymax = top_subgroup),
            fill = "white", 
            color = "transparent") + 
  geom_line(data = df_mod2, 
            mapping = aes(x = ts,
                          y = sg_jitter,
                          group = as.factor(coati_ids.name),
                          color = as.factor(coati_ids.name)),
            alpha = 0.3,
            size = 1) +
  geom_point(data = df_mod2,
             mapping = aes(x = ts,
                           y = sg_jitter,
                           color = as.factor(coati_ids.name)),
             size = 2.2,
             alpha = 1) +
  geom_rect(data = df_group_polys,
            mapping = aes(xmin = left_subgroup,
                          xmax = right_subgroup,
                          ymin = bottom_subgroup,
                          ymax = top_subgroup),
            fill = "transparent", #transparent fill so that they don't block otu the coati points and lines
            color = "grey40") + #dark grey outline
  scale_color_manual(values = sample(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#b15928'), #these are 11 distinct colors, they will be randomly assigned to the coatis
                                     11, replace = FALSE)) +
  theme_classic() +
  #get rid of all the extra stuff (note, this now is for a horizontal plot. if you want a vertical plot,
  # then replace the xs with ys and the ys with xs)
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y =element_blank(),
        axis.text.x = element_text(size = 20, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x.bottom=element_line(color="black"),
        axis.ticks.y = element_blank()) +
  NULL

p1

# save it
ggsave(filename = "C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/FF line and dot subgroupings_0501022_horz.png",
       width = 15,
       height = 3,
       units = "in",
       dpi = 350,
       scale = 1)



























