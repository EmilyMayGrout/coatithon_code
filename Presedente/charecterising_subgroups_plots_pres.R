#plots for the low res fission fusion paper

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level1.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 50

#-------SETUP-------
#making plots for characterizing the subgroups

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
#removing males and wildflower
xs <- xs[c(1:4,6,7,11:17,19,20,22),]
ys <- ys[c(1:4,6,7,11:17,19,20,22),]
coati_ids <- coati_ids[-c(5,8,9,10,18,21),]

n_inds <- nrow(xs)
n_times <- ncol(xs)



#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#make 3 plots next to each other, first is hist of number of subgroups at 50m, next it distribution of individuals in each group size with 2 subgroups, last is at 3 subgroups

png(height = 400, width = 1280, units = 'px', filename = paste0(plot_dir,'50m_charecterisations.png'))

par(mfrow=c(1,3), mar = c(6,5,2,1)) #(bottom, left, top, right)
subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs],main = "50m radii", xlab =  'Number of subgroups', col = "darkolivegreen3", breaks = seq(.5,22,1), cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.5), xlim = c(0, 6))

subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], breaks=seq(0.5,22,1), xlab = 'Subgroup size', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.4), xlim = c(0, 16), xaxp = c(1,16, 15))
hist(subgroup_counts[,s3], breaks=seq(0.5,22,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.4), xlim = c(0, 16), xaxp = c(1,16, 15))

dev.off()

#just for the 2 subgroups plot

png(height = 500, width = 600, units = 'px', filename = paste0(plot_dir,'50m_charecterisations_2groups.png'))
par( mar = c(6,5,2,1))
hist(subgroup_counts[,s2], breaks=seq(0.5,22,1), xlab = 'Subgroup size', main = '', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 16), xaxp = c(1,16, 15))
dev.off()





