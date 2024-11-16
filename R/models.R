# test

# packages

#library(pacman)

#remotes::install_github("stan-dev/cmdstanr")
install_cmdstan(cores = 18)

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

dat <- read.csv(here("data", "corvidae_310spp.csv"))

str(dat)

# reading nexus file (use ape pacakge)

tree_all <- read.nexus(here("data", "corvidae_310spp_100.nex"))
tree <- tree_all[[1]]

# trun tree into correlation matrix using vcv function

A <- vcv.phylo(tree, corr = TRUE)

# center data of interest

dat$cbeak_length <- scale(log(dat$Beak.Length_Culmen), center = TRUE, scale = FALSE)
dat$cbeak_width <- scale(log(dat$Beak.Width), center = TRUE, scale = FALSE)
dat$cbeak_depth <- scale(log(dat$Beak.Depth), center = TRUE, scale = FALSE)

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

#rebuild_cmdstan()

fit1 <- brm(formula1, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            prior = prior1,
            backend = "cmdstanr",
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
            )

summary(fit1)

########################
# location-scale model 2 using brms with phylogenetic correlation
########################

# length
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
            iter = 3000, 
            warmup = 2000,
            backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2)

# width

formula2b <- bf(cbeak_width ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

fit2b <- brm(formula2b, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2b)

# depth

formula2c <- bf(cbeak_depth ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

fit2c <- brm(formula2c, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2c)


########################
# location-scale model 2 using brms with phylogenetic correlation
########################

formula3A <- bf(cbeak_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula3B <- bf(cbeak_depth ~1 + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula3 <- formula3A + formula3B + set_rescor(TRUE) #do we need to do correlction??

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
            iter = 3000, 
            warmup = 2000,
            prior = prior3,
            backend = "cmdstanr",
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit3)

