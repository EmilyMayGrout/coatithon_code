
#this code is for trying out the STRAND package with the baboon data from ctross/STRAND github

########################################
#
#   Binomial Analyses  
#
########################################

# Clear working space
rm(list = ls())


# Load libraries
library(ggplot2)

#library(devtools)
#install_github("ctross/STRAND", force = TRUE)

library(STRAND)
library(cmdstanr)

#set_cmdstan_path(path = "C:/Users/egrout/Documents/R/win-library/.cmdstan/cmdstan-2.32.2")

# Import baboon data for practice
load("C:/Users/egrout/Dropbox/coatithon/processed/strand_test_data/Baboon_Data.RData")

# Number of grooming event and a sample-size measure
# Here, the term "exposure" relates to the number of trials for a binomial distribution
# for coatis this would be the subgroup membership
nets = list(Grooming = Baboon_Data$Grooming)
# number of observations (so in coatis, this is number of times both individuals have a gps point)
exposure_nets = list(Exposure = Baboon_Data$Exposure)

# Dyadic variable: transpose of Presenting 
# for coatis this would be relatedness
dyad = list(Presenting = t(Baboon_Data$Presenting),
            Threatening = t(Baboon_Data$Threatening))

#block-structuring variables are individual level data but must be factors as they are used to create a random intercept unique to the interaction of the focal/sender
block = data.frame(Sex = as.factor(Baboon_Data$Sex)) 

#individual level covariates data can be numeric variables, indicator variables, or categorical variables 
indiv =  data.frame(Age = Baboon_Data$Age) 

model_dat = make_strand_data(outcome = nets,
                             individual_covariates = indiv, 
                             block_covariates = block,
                             dyadic_covariates = dyad,
                             outcome_mode = "binomial",
                             exposure = exposure_nets
)

# model - this is where I get the error
fit =  fit_block_plus_social_relations_model(data=model_dat,
                                             block_regression = ~ Sex,
                                             focal_regression = ~ Age,
                                             target_regression = ~ Age,
                                             dyad_regression = ~  Presenting + Threatening,
                                             mode="mcmc",
                                             stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,iter_warmup = 1000, iter_sampling = 1000, max_treedepth = NULL, adapt_delta = .98)
)

res = summarize_strand_results(fit)


############################### Plots
vis_1 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=TRUE, site="HP", only_slopes=TRUE)
vis_1
#ggsave("Baboon_slopes.pdf", vis_1, width=7, height=2.5)

vis_2 = strand_caterpillar_plot(res, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects","Other estimates"), normalized=FALSE, site="HP", only_technicals=TRUE, only_slopes=FALSE)
vis_2
#ggsave("Baboon_corr.pdf", vis_2, width=6, height=2.5)


