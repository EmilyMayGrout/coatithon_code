#fission fusion analysis script
#first define groups for each time step (each 10 mins)

library(dbscan)

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
in_file <- "galaxy_xy_10min_level0.RData"

#radius for group membership
R <- 30

#-------MAIN-------
setwd(data_dir)

load(in_file)

t_idx <- 190
x <- xs[ , t_idx]
y <- ys[, t_idx]

plot(x, y)

groups <- dbscan(x= cbind(x, y), eps = R, minPts = 1)

str(groups)

plot(x, y, col = groups$cluster)

