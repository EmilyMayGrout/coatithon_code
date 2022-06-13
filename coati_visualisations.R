
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy_xy_10min_level0.RData")

plot(xs[,1001], ys[,1001])

plot(NULL, xlim = c(min(xs, na.rm=T), max(xs, na.rm=T)), ylim = c(min(ys, na.rm=T), max(ys, na.rm=T)), asp = 1)
points(xs[1,], ys[1,])
points(xs[5,], ys[5,], col = "blue")
points(xs[2,], ys[2,], col = "red")
