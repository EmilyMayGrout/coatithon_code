#this code is for trying out the STRAND package 

########################################
#
#   Binomial Analyses  
#
########################################

# Clear working space
rm(list = ls())


# Load libraries
#issue is that I can't download stan because there isn't enough memory 
library(ggplot2)

library(devtools)
install_github("ctross/STRAND", force = TRUE)

library(STRAND)
library(cmdstanr)

#with coati galaxy data
galaxy <- readRDS("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/array_mrqap.rds")

#make into the same format as the Baboon_Data
nets2 = list(Subgroup_membership = galaxy[,,1])

#is exposure needed for the coatis when the proportion of time each individual is in the same subgroup is accounting for ?
# need to ask Ari about this... and if its not accounted for, I need help working out how to calculate this
# for now having exposure as NA

#dyad for relatedness
dyad2 = list(Relatedness = t(galaxy[,,2]))

#getting sex into correct format for block2
gal_sex <- galaxy[,,3]
gal_sex <- gal_sex[1,]
gal_sex[is.na(gal_sex)] <- 1
gal_sex[gal_sex == 0] <- "female"
gal_sex[gal_sex == 1] <- "male"


block2 = data.frame(Sex = as.factor(gal_sex)) 

#keeping in same order as the galaxy matrix
indiv2 = data.frame(Age = c(3, 1, 3, 3, 3, 3, 3, 3, 2, 2, 2)) 

model_dat2 = make_strand_data(outcome = nets2,
                              individual_covariates = indiv2, 
                              block_covariates = block2,
                              dyadic_covariates = dyad2,
                              outcome_mode = "binomial",
                              exposure = NA
)

# model
fit2 =  fit_block_plus_social_relations_model(data=model_dat2,
                                              block_regression = ~ Sex,
                                              focal_regression = ~ Age,
                                              target_regression = ~ Age,
                                              dyad_regression = ~  Relatedness,
                                              mode="mcmc",
                                              stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                                          iter_warmup = 1000, iter_sampling = 1000,
                                                                          max_treedepth = NULL, adapt_delta = .98)
)


