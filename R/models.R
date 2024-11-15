# test

# packages

#library(pacman)

pacman::p_load(tidyverse, 
       ggplot2,
       brms,
       rstan,
       rstanarm,
       loo,
       bayesplot,
       here,
       ape
)

# load data

dat <- read.csv(here("data", "corvidae_310spp.csv"))

str(dat)

# reading nexus file (use ape pacakge)

tree_all <- read.nexus(here("data", "corvidae_310spp_100.nex"))
tree <- tree_all[[1]]

# trun tree into correlation matrix using vcv function

A <- vcv.phylo(tree, corr = TRUE)

# center data of interest

dat$cbeak_length <- scale(dat$Beak.Length_Culmen, center = TRUE, scale = FALSE)
dat$cbeak_width <- scale(dat$Beak.Width, center = TRUE, scale = FALSE)
dat$cbeak_depth <- scale(dat$Beak.Depth, center = TRUE, scale = FALSE)

#########################
# location-scale model 1 using brms with phylogenetic correlation
########################

formula1 <- bf(cbeak_length ~1 + (1|gr(Phylo, cov = A)), 
                sigma ~ 1)


# creat prior

prior1 <- default_prior(formula1, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
                        )

# fit model

fit1 <- brm(formula1, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 4000, 
            warmup = 2000,
            prior = prior1,
            threads = threading(8)
            #control = list(adapt_delta = 0.95, max_treedepth = 15)
            )

summary(fit1)

########################
# location-scale model 2 using brms with phylogenetic correlation
########################

formula2 <- bf(cbeak_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)


# creat prior

prior2 <- default_prior(formula2, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

# fit model

fit2 <- brm(formula2, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 4000, 
            warmup = 2000,
            prior = prior2,
            threads = threading(8)
            #control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2)

########################
# location-scale model 2 using brms with phylogenetic correlation
########################

formula3A <- bf(cbeak_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula3B <- bf(cbeak_depth ~1 + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula3 <- formula3A + formula3B

# creat prior

prior3 <- default_prior(formula3, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

# fit model

fit3 <- brm(formula3, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 2000, 
            warmup = 1000,
            prior = prior3,
            threads = threading(9)
            #control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit3)

