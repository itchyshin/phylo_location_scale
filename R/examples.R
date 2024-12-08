# key models - examples

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

# create a variable which says whether they live in Forest or not

dat$forest <- ifelse(dat$Habitat == "Forest", 1, 0)

# Different varainces for different groups (Model 3)

form1 <- bf(cbeak_length ~ 1 + cmass + (1|p|gr(Phylo, cov = A)), 
            sigma ~ 1 + forest 
)

prior1 <- default_prior(form1, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod1 <- brm(form1, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 20000, 
            warmup = 5000,
            #backend = "cmdstanr",
            prior = prior1,
            threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod1)

# saving

saveRDS(mod1, here("Rdata", "mod1.rds"))


# loading and checking

mod1 <- readRDS(here("Rdata", "mod1.rds"))

summary(mod1)


# Mean-variance correlation (Model 4)

form2 <- bf(crange_size ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
            sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

prior2 <- default_prior(form1, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)

mod2 <- brm(form2, 
            data = dat, 
            data2 = list(A = A),
            chains = 2, 
            cores = 2, 
            iter = 6000, 
            warmup = 3000,
            #backend = "cmdstanr",
            prior = prior2,
            threads = threading(9),
            control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod2)

# saving

saveRDS(mod2, here("Rdata", "mod2.rds"))

# loading and checking

mod2 <- readRDS(here("Rdata", "mod2.rds"))

summary(mod2)


# Co-evolution of means and variances (Model 5)

form3A <- bf(cbeak_width ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

form3B <- bf(cbeak_depth ~1 + cmass + (1|p|gr(Phylo, cov = A)), 
                sigma ~ 1 + cmass + (1|p|gr(Phylo, cov = A))
)

form3 <- form3A + form3B + set_rescor(TRUE) #do we need to do correction??

# creat prior

prior3 <- default_prior(form3, 
                        data = dat, 
                        data2 = list(A = A),
                        family = gaussian()
)



# fit model

mod3 <- brm(form3, 
                  data = dat, 
                  data2 = list(A = A),
                  chains = 2, 
                  cores = 2, 
                  iter = 6000, 
                  warmup = 4000,
                  prior = prior3,
                  #backend = "cmdstanr",
                  threads = threading(9),
                  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(mod3)

# saving

saveRDS(mod3, here("Rdata", "mod3.rds"))

# loading and looking at it again

mod3 <- readRDS(here("Rdata", "mod3.rds"))

summary(mod3)

