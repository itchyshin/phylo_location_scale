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

dat <- read.csv(here("data", "Psittaciformes", "Psittaciformes_354spp.csv"))

str(dat)

# reading nexus file (use ape pacakge)

tree_all <- read.nexus(here("data", "Psittaciformes", "Psittaciformes_354spp_100.nex"))
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


# range size

formula1 <- bf(crange_size ~1 + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior1 <- default_prior(formula1, 
                         data = dat, 
                         data2 = list(A = A),
                         family = gaussian()
)

mod_psit_1 <- brm(formula1, 
             data = dat, 
             data2 = list(A = A),
             chains = 2, 
             cores = 2, 
             iter = 3000, 
             warmup = 2000,
             #backend = "cmdstanr",
             prior = prior1,
             threads = threading(9),
             control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_1)

# saving

saveRDS(mod_psit_1, here("Rdata", "mod_psit_1.rds"))


########################
# location-scale model 2 using brms with phylogenetic correlation
########################

formula3A <- bf(cbeak_width ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

formula3B <- bf(cbeak_depth ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

formula3 <- formula3A + formula3B + set_rescor(TRUE) #do we need to do correction??

# creat prior

prior3 <- default_prior(formula3, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

# fit model

mod_psit_2 <- brm(formula3, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 5000, 
            warmup = 3000,
            prior = prior3,
            #backend = "cmdstanr",
            threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_2)

# saving

saveRDS(mod_psit_2, here("Rdata", "mod_psit_2.rds"))

