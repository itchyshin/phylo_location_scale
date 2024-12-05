# model 3
# test

# packages

#library(pacman)

#remotes::install_github("stan-dev/cmdstanr")
#install_cmdstan(cores = 18)

pacman::p_load(tidyverse, 
               ggplot2,
               brms,
               rstan,
               rstanarm,
               loo,
               bayesplot,
               here,
               ape,
               cmdstanr
)

# load data

dat <- read.csv(here("data/Tyrannides/Tyrannides_1175spp.csv"))

str(dat)

# reading nexus file (use ape pacakge)
tree_all <- read.nexus(here("data/Tyrannides/tree_100/Tyrannides_1175spp_100.nex"))
tree <- tree_all[[1]]
# trun tree into correlation matrix using vcv function

A <- vcv.phylo(tree, corr = TRUE)

# center data of interest

dat$cbeak_length <- scale(log(dat$Beak.Length_Culmen), center = TRUE, scale = FALSE)
dat$cbeak_width <- scale(log(dat$Beak.Width), center = TRUE, scale = FALSE)
dat$cbeak_depth <- scale(log(dat$Beak.Depth), center = TRUE, scale = FALSE)
dat$cmass <- scale(log(dat$Mass), center = TRUE, scale = FALSE)
dat$cwing_length <- scale(log(dat$Wing.Length), center = TRUE, scale = FALSE)
dat$ctail_length <- scale(log(dat$Tail.Length), center = TRUE, scale = FALSE)
dat$tarsus_length <- scale(log(dat$Tarsus.Length), center = TRUE, scale = FALSE)
dat$crange_size <- scale(log(dat$Range.Size), center = TRUE, scale = FALSE)

# correlaiton check

cor(dat$cbeak_length, dat$cbeak_width)
cor(dat$cbeak_length, dat$cbeak_depth)
cor(dat$cbeak_width, dat$cbeak_depth)

#########################
# location-scale model 1 using brms with phylogenetic correlation
########################


# wing 

formula1 <- bf(cwing_length ~1 + cmass +  (1|p|gr(Phylo, cov = A)), 
               sigma ~ Family - 1 + cmass
)

prior1 <- default_prior(formula1, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_tf_1 <- brm(formula1, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  backend = "cmdstanr",
                  prior = prior1,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_tf_1)

# saving

saveRDS(mod_tf_1, here("Rdata", "mod_tf_1.rds"))


# beak length

formula2 <- bf(clength ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
               sigma ~ Family - 1 + cmass
)

prior2 <- default_prior(formula2, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_tf_2 <- brm(formula2, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  backend = "cmdstanr",
                  prior = prior2,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_tf_2)

# saving

saveRDS(mod_tf_2, here("Rdata", "mod_tf_2.rds"))


