#this script is for making a dataframe with time, and list group members in original group and then the sub groups for each - this is to then be used for the ff over time graph (help from Alie)


data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

library(fields)
library(viridis)
library(tidyverse)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, 50)

subgroup_data$ind_subgroup_membership

t <- as.data.frame(rbind(subgroup_data$ind_subgroup_membership[,1:1633]))
#put id into dataframe
t <- cbind(coati_ids$name, t)
#changing column names into a list from 1 to 1633 - one day would be good to have as time
column_names <- c(1:1633)
colnames(t)[-1] <- column_names

library(tidyr)
#pivot the data frame into a long format
test <- t %>% pivot_longer(-c(`coati_ids$name`), names_to='time', values_to='sub-group')
test$time <- as.numeric(test$time)

#time df
time_df <- data.frame(ts, 1:length(ts))
names(time_df)[2] <- "time"
pivot_subgroups <- left_join(test, time_df)

#rename column names

colnames(pivot_subgroups) <- c("name", "time", "subgroup_ID", "ts")


#---------------------------------------------------------------------------------------
#Alie's example code:
library(tidyverse)


df_start <- data.frame(Time = rep(c(1:10), each = 5),
                       ID = rep(c("A", "B", "C", "D", "E"), times = 5),
                       SubGroup = sample(1:3, 50, replace = T))
boom <- df_start %>% group_by(Time, SubGroup) %>% nest(data = ID) %>% ungroup()
boom2 <- boom %>% pivot_wider(names_from = SubGroup, values_from = data)
all <- df_start %>% select(-SubGroup) %>% group_by(Time) %>% nest(data = ID) 
output <- full_join(boom2, all)

#--------------------------------------------------------------------------------------

## THIS WORKS!! ##

nesting_df <- pivot_subgroups %>% group_by(time, subgroup_ID) %>% nest(data = name) %>% ungroup()
pivot_nesting <- nesting_df %>% pivot_wider(names_from = subgroup_ID, values_from = data)
all_inds <- pivot_subgroups %>% select(-subgroup_ID) %>% group_by(time) %>% nest(data = name)
nest_df_all <- full_join(pivot_nesting, all_inds)

save(nest_df_all, file = paste0(data_dir,'galaxy_nest_subgroups_lowres_level0.RData'))

#---------------------------------------------------------------------
















