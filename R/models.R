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
dat$cmass <- scale(log(dat$Mass), center = TRUE, scale = FALSE)
dat$cwing_length <- scale(log(dat$Wing.Length), center = TRUE, scale = FALSE)
dat$ctail_length <- scale(log(dat$Tail.Length), center = TRUE, scale = FALSE)
dat$tarsus_length <- scale(log(dat$Tarsus.Length), center = TRUE, scale = FALSE)
dat$crange_size <- scale(log(dat$Range.Size), center = TRUE, scale = FALSE)

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
            iter = 5000, 
            warmup = 3000,
            prior = prior1,
            #backend = "cmdstanr",
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
            )

summary(fit1)


# putting body mass into the model

formula1b <- bf(cbeak_length ~1 + cmass (1|gr(Phylo, cov = A)), 
                sigma ~ 1)

prior1b <- default_prior(formula1b, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit1b <- brm(formula1b, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 2000, 
            warmup = 1000,
            prior = prior1b,
            #backend = "cmdstanr",
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

########################
# location-scale model 2 using brms with phylogenetic correlation
########################

# length
formula2 <- bf(cbeak_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)


# create prior

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
            #backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2)

# width

formula2b <- bf(cbeak_width ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2b <- default_prior(formula2b, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)


fit2b <- brm(formula2b, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2b,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2b)

# depth

formula2c <- bf(cbeak_depth ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2c <- default_prior(formula2c, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)


fit2c <- brm(formula2c, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2c)

# mass

formula2d <- bf(cmass ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2d <- default_prior(formula2d, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit2d <- brm(formula2d, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2d,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2d)

# wing length

formula2e <- bf(cwing_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2e <- default_prior(formula2e, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit2e <- brm(formula2e, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2e,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2e)

# tail length

formula2f <- bf(ctail_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2f <- default_prior(formula2f, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit2f <- brm(formula2f, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2f,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2f)

# tarsus length

formula2g <- bf(tarsus_length ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2g <- default_prior(formula2g, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit2g <- brm(formula2g, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2g,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2g)

# range size

formula2h <- bf(crange_size ~1 + (1|p|gr(Phylo, cov = A)), 
               sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

prior2h <- default_prior(formula2h, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

fit2h <- brm(formula2h, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 3000, 
            warmup = 2000,
            #backend = "cmdstanr",
            prior = prior2h,
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit2h)




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
            iter = 5000, 
            warmup = 3000,
            prior = prior3,
            backend = "cmdstanr",
            threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit3)

# saving

saveRDS(fit1, here("Rdata", "fit3.rds"))


# fit 4

formula4A <- bf(cbeak_width ~1 + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula4B <- bf(cbeak_depth ~1 + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + (1|p|gr(Phylo, cov = A))
)

formula4 <- formula4A + formula4B + set_rescor(TRUE) #do we need to do correlction??

# creat prior

prior4 <- default_prior(formula4, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

# fit model

fit4 <- brm(formula4, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 5000, 
            warmup = 3000,
            prior = prior4,
            #backend = "cmdstanr",
            #threads = threading(9),
            control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(fit4)

# save models rds

saveRDS(fit1, here("Rdata", "fit4.rds"))




