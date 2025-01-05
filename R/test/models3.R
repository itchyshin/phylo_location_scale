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
               cmdstanr,
               bayesplot, 
               tidybayes,
               patchwork,
               phytools,
               coda,
               MASS
)


# pacman::p_load("coda", "tidyverse", "here", 
#                "MCMCglmm", "brms", "MASS", 
#                "phytools", "patchwork", "bayesplot", "tidybayes")

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

# create a variable which says whether they live in Forest or not

dat$forest <- ifelse(dat$Habitat == "Forest", 1, 0)


#########################
# location-scale model 1 using brms with phylogenetic correlation
########################


# range size

formula1 <- bf(crange_size ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
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

formula2A <- bf(cbeak_width ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

formula2B <- bf(cbeak_depth ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

formula2 <- formula2A + formula2B + set_rescor(TRUE) #do we need to do correction??

# creat prior

prior2 <- default_prior(formula2, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)



# fit model

mod_psit_2 <- brm(formula2, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 5000, 
            warmup = 3000,
            prior = prior2,
            #backend = "cmdstanr",
            threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 20)
)

summary(mod_psit_2)

# saving

saveRDS(mod_psit_2, here("Rdata", "mod_psit_2.rds"))

# loading and looking at it again

mod_psit_2 <- readRDS(here("Rdata", "mod_psit_2.rds"))

summary(mod_psit_2)

###########
# range

formula3 <- bf(crange_size ~ 1 + cmass + (1|gr(Phylo, cov = A)), 
               sigma ~ 1 + forest + cmass
)

prior3 <- default_prior(formula3, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_psit_3 <- brm(formula3, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  #backend = "cmdstanr",
                  prior = prior3,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_3)

# saving

saveRDS(mod_psit_3, here("Rdata", "mod_psit_3.rds"))

###############
# beak length
###############


formula4 <- bf(cbeak_length ~ 1 + cmass + (1|gr(Phylo, cov = A)), 
               sigma ~ 1 + forest + cmass
)

prior4 <- default_prior(formula4, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_psit_4 <- brm(formula4, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  #backend = "cmdstanr",
                  prior = prior4,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_4)

# saving

saveRDS(mod_psit_4, here("Rdata", "mod_psit_4.rds"))

# beak width

formula5 <- bf(cbeak_width ~ 1 + cmass + (1|gr(Phylo, cov = A)), 
               sigma ~ 1 + forest + cmass
)

prior5 <- default_prior(formula5, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_psit_5 <- brm(formula5, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  #backend = "cmdstanr",
                  prior = prior5,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_5)

# saving

saveRDS(mod_psit_5, here("Rdata", "mod_psit_5.rds"))

# beak depth

formula6 <- bf(cbeak_depth ~ 1 + cmass + (1|gr(Phylo, cov = A)), 
               sigma ~ 1 + forest + cmass
)

prior6 <- default_prior(formula6, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod_psit_6 <- brm(formula6, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 3000, 
                  warmup = 2000,
                  #backend = "cmdstanr",
                  prior = prior6,
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod_psit_6)

# saving

saveRDS(mod_psit_6, here("Rdata", "mod_psit_6.rds"))

