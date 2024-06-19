#running the same fission-fusion code with Trago group 

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/trago_results/'
gps_file <- "trago_xy_10min_level0.RData"
id_file <- 'trago_coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

t <- coati_ids

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


#------plot 1: number of sub groups when the radius is changed (graph put in dropbox results folder) -----------
png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists_trago.png'))
par(mfrow=c(6,1), mar = c(6,5,1,1))

for (i in 1:length(Rs)){
  
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7))
  
}
dev.off()



#------plot 2: number of individuals in each sub group when radius is 50m -----------

png(height = 540, width = 270, units = 'px', filename = paste0(plot_dir,'subgroup_size_hists_50m_trago.png'))
subgroup_data <- get_subgroup_data(xs, ys, R=50)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)


par(mfrow=c(2,1), mar = c(6,5,1,1))
hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))

dev.off()

#-------------------------------------------------------------------------
# plot 3
#Figure S1 in ff low res manuscript

subgroup_data <- get_subgroup_data(xs, ys, R=50)

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
new_order <- c(2,3,7:9,1,5,4,6) #order for matrices for microbiome
#"Amarulla","Limoncello","Whisky","Cerveza","Sake","Tiger","Bailey","Rum","Tequila" 

new_order <- c(4,1:3, 5:7)
ffnet_reorder <- ff_net[new_order, new_order]


png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_trago.png'))
par(mfrow=c(1,1), mar = c(6,6,1,1)) #bottom, left, top, and right
visualize_network_matrix_trago(ffnet_reorder, coati_ids[new_order,])
dev.off()

#save matrix for mrqap analysis
write.table(ffnet_reorder,file="C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/trago_matrix_50min_proptimeinsamesubgroup.txt",row.names=FALSE)


#---------------------------------------------------------------------------------------------
#look at which age/sex classes tend to be on their own
inds_subgroup <- data.frame(subgroup_data$ind_subgroup_membership)
#make an empty dataframe to add the alone inds data to
df <- data.frame(matrix(nrow = 9, ncol = ncol(inds_subgroup)))

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

gg

coati_ids_alone_trago <- coati_ids
#this dataframe is used in alone_inds_all_groups
save(coati_ids_alone_trago, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/trago_alone_inds_level1.Rdata")  

ggsave(filename = paste0(plot_dir, 'prop_time_alone.png'), plot = gg, width = 6, height = 6, dpi = 300)
