# this script is to make a plot with subgroup size and speed of travel

library(move)
library(data.table)
require(EMbC)
library(amt)
library(RColorBrewer)
library(zoo)
library(dplyr)
library(ggplot2)
library(lubridate)
library(arm)  #also loads library lme4
library(blmeco)

plot_dir <- "C:/Users/egrout/Dropbox/stats_Franzi/plots/"
data_dir <- "C:/Users/egrout/Dropbox/stats_Franzi/data/"

#this data is from lowres embc script in the animove folder in dropbox
load(file = paste0(data_dir,"resample_Galaxy.Rdata"))

#convert resample_Galaxy to a movestack
load(file = paste0(data_dir,"resample_Galaxy.Rdata"))
resamp_gal <- moveStack(resample_Galaxy)
resamp_gal <- spTransform(resamp_gal, CRS("+init=epsg:32617"))

#--------------------------------------------------------
# make a dataframe with ID, time, speed, age and sex

#first collapse the list

Galaxy_df <- as.data.frame(resamp_gal)
names(Galaxy_df)
Galaxy_df <- Galaxy_df[,c(35, 36, 27, 28, 29, 30, 31, 32)]
#read in the coati ids with the age and sex class
coati_id <- read.csv(paste0(data_dir, "coati_id.csv"), header = F)

#add the age and sex to each individual using the match function
Galaxy_df$age <- coati_id$V3[match(Galaxy_df$trackId, coati_id$V1)]

Galaxy_df$sex <- coati_id$V4[match(Galaxy_df$trackId, coati_id$V1)]

#define the age and sex classes as factors
Galaxy_df$age_logical <- factor(Galaxy_df$age, levels = c("Juvenile", "Sub-adult", "Adult"))
Galaxy_df$age <- factor(Galaxy_df$age)
Galaxy_df$sex <- as.factor(Galaxy_df$sex)

#add column for hour
Galaxy_df$hours <- factor(hour(Galaxy_df$timestamps))

#now open the dataframe made in coatithon/fission_fusion (around line 300) called n_subs

load(file = paste0(data_dir,"n_subs.Rdata"))

#merge n_subs and Galaxy_df
Galaxy_df$hours <- as.integer(Galaxy_df$hours)

#not possible to merge as Galaxy_df has diff individuals so it wouldn't make sense. Not sure what to do now...








