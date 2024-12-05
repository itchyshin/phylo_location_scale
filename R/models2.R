# brain, body size etc model

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
# reading xlsx file

dat <- readxl::read_xlsx(here("data", "201126_all_species_corrected.xlsx"))

str(dat)


dat <- select(dat, tip_label, brain, weight, devo_mode)

# reading nexus file (use ape pacakge)

tree <- read.nexus(here("data", "Jetztree1.nex"))

# trimming tree

#tree <- drop.tip(tree, tree$tip.label[dat$tip_label %in% tree$tip.label])

retain <- match(tree$tip.label, dat$tip_label)
keep <- is.na(retain)

tree <- drop.tip(tree, tree$tip.label[keep])               
# trun tree into correlation matrix using vcv function

A <- vcv.phylo(tree, corr = TRUE)

exclude <- setdiff(dat$tip_label, row.names(A))

# centering the data

dat$cbrain <- scale(log(dat$brain), center = TRUE, scale = FALSE)
dat$cweight <- scale(log(dat$weight), center = TRUE, scale = FALSE)

dat_full <- dat %>% 
  filter(!tip_label %in% exclude) %>%
  filter(!is.na(cweight)) 

# modeling

from1 <- bf(cbrain ~1 + devo_mode + cweight + (1|a|gr(tip_label, cov = A)), 
               sigma ~ 1 + devo_mode + cweight)


# create prior

prior1 <- default_prior(from1, 
                        data = dat_full, 
                        data2 = list(A = A),
                        family = gaussian()
)


mod1 <- brm(from1, 
               data = dat_full, 
               data2 = list(A = A),
               chains = 2, 
               cores = 2, 
               iter = 3000, 
               warmup = 2000,
               #backend = "cmdstanr",
               prior = prior1,
               threads = threading(9),
               control = list(adapt_delta = 0.95, max_treedepth = 15)
)

summary(mod1)

# save mod1

save(mod1, file = here("Rdata", "mod1.RData"))

