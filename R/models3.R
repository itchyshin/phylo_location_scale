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


# range

formula3 <- bf(crange_size ~ 1 + cmass + (1|p|gr(Phylo, cov = A)), 
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

# beak length

formula4 <- bf(cbeak_length ~ 1 + cmass + (1|p|gr(Phylo, cov = A)), 
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

formula5 <- bf(cbeak_width ~ 1 + cmass + (1|p|gr(Phylo, cov = A)), 
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

formula6 <- bf(cbeak_depth ~ 1 + cmass + (1|p|gr(Phylo, cov = A)), 
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


# chatGPT code

# load all the models rds

# Load models
mod_psit_1 <- readRDS(here("Rdata", "mod_psit_1.rds"))
mod_psit_2 <- readRDS(here("Rdata", "mod_psit_2.rds"))
mod_psit_3 <- readRDS(here("Rdata", "mod_psit_3.rds"))
mod_psit_4 <- readRDS(here("Rdata", "mod_psit_4.rds"))
mod_psit_5 <- readRDS(here("Rdata", "mod_psit_5.rds"))
mod_psit_6 <- readRDS(here("Rdata", "mod_psit_6.rds"))

library(tidyverse)
library(tidybayes)
library(ggplot2)
library(patchwork)

# Function to get variable names dynamically
get_variables_dynamic <- function(model, pattern) {
  variables <- get_variables(model)
  variables[grep(pattern, variables)]
}


rename_vars <- function(variable) {
  variable <- gsub("b_", "", variable)                # Remove "b_"
  variable <- gsub("sd_", "SD_", variable)            # Replace "sd_" with "SD_"
  variable <- gsub("cor_", "Correlation_", variable)  # Replace "cor_" with "Correlation_"
  variable <- gsub("_", " ", variable)                # Replace "_" with space
  return(variable)
}

# Function to visualize fixed effects
visualize_fixed_effects <- function(model, title) {
  fixed_effect_vars <- get_variables_dynamic(model, "^b_")
  if (length(fixed_effect_vars) == 0) {
    message("No fixed effects found for ", title)
    return(NULL)
  }
  
  tryCatch({
    fixed_effects_samples <- model %>%
      spread_draws(!!!syms(fixed_effect_vars)) %>%
      pivot_longer(cols = all_of(fixed_effect_vars), names_to = ".variable", values_to = ".value")
    
    ggplot(fixed_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "lightcyan3", 
        color = "lightcyan4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(title = title, y = "Fixed Effects", x = "Posterior Values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects: ", e$message)
    return(NULL)
  })
}

# Function to visualize fixed effects
visualize_fixed_effects <- function(model, title) {
  fixed_effect_vars <- get_variables_dynamic(model, "^b_")
  if (length(fixed_effect_vars) == 0) {
    message("No fixed effects found for ", title)
    return(NULL)
  }
  
  tryCatch({
    fixed_effects_samples <- model %>%
      spread_draws(!!!syms(fixed_effect_vars)) %>%
      pivot_longer(cols = all_of(fixed_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable))
    
    ggplot(fixed_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "lightcyan3", 
        color = "lightcyan4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(title = title, y = "Fixed Effects", x = "Posterior Values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects: ", e$message)
    return(NULL)
  })
}

# Function to visualize random effects
visualize_random_effects <- function(model, title) {
  random_effect_vars <- get_variables_dynamic(model, "^sd_")
  if (length(random_effect_vars) == 0) {
    message("No random effects found for ", title)
    return(NULL)
  }
  
  tryCatch({
    random_effects_samples <- model %>%
      spread_draws(!!!syms(random_effect_vars)) %>%
      pivot_longer(cols = all_of(random_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable)) %>%
      mutate(.value = .value^2)  # Convert SD to variance for interpretation
    
    ggplot(random_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "olivedrab3", 
        color = "olivedrab4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(title = title, y = "Random Effects (Variance)", x = "Posterior Values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_random_effects: ", e$message)
    return(NULL)
  })
}

# Function to visualize correlations
visualize_correlations <- function(model, title) {
  correlation_vars <- get_variables_dynamic(model, "^cor_")
  if (length(correlation_vars) == 0) {
    message("No correlations found for ", title)
    return(NULL)
  }
  
  tryCatch({
    correlation_samples <- model %>%
      spread_draws(!!!syms(correlation_vars)) %>%
      pivot_longer(cols = all_of(correlation_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable))
    
    ggplot(correlation_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        fill = "#FF6347", 
        color = "#8B3626"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(title = title, y = "Correlations", x = "Posterior Values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_correlations: ", e$message)
    return(NULL)
  })
}

# Visualizing each model
plots_model1 <- list(
  visualize_fixed_effects(mod_psit_1, "Model 1: Fixed Effects"),
  visualize_random_effects(mod_psit_1, "Model 1: Random Effects"),
  visualize_correlations(mod_psit_1, "Model 1: Correlations")
)

plots_model2 <- list(
  visualize_fixed_effects(mod_psit_2, "Model 2: Fixed Effects"),
  visualize_random_effects(mod_psit_2, "Model 2: Random Effects"),
  visualize_correlations(mod_psit_2, "Model 2: Residual Correlations")
)

plots_model3 <- list(
  visualize_fixed_effects(mod_psit_3, "Model 3: Fixed Effects"),
  visualize_random_effects(mod_psit_3, "Model 3: Random Effects")
)

plots_model4 <- list(
  visualize_fixed_effects(mod_psit_4, "Model 4: Fixed Effects"),
  visualize_random_effects(mod_psit_4, "Model 4: Random Effects")
)

plots_model5 <- list(
  visualize_fixed_effects(mod_psit_5, "Model 5: Fixed Effects"),
  visualize_random_effects(mod_psit_5, "Model 5: Random Effects")
)

plots_model6 <- list(
  visualize_fixed_effects(mod_psit_6, "Model 6: Fixed Effects"),
  visualize_random_effects(mod_psit_6, "Model 6: Random Effects")
)

# Create a placeholder plot for models without correlations
placeholder_plot <- ggplot() +
  theme_void() +
  labs(title = "No Correlations Available")

# Arrange plots in 6 rows (one for each model) and 3 columns (Fixed, Random, Correlations)
all_plots <- patchwork::wrap_plots(
  list(
    # Row 1 (Model 1)
    plots_model1[[1]], plots_model1[[2]], plots_model1[[3]],
    # Row 2 (Model 2)
    plots_model2[[1]], plots_model2[[2]], plots_model2[[3]],
    # Row 3 (Model 3)
    plots_model3[[1]], plots_model3[[2]], placeholder_plot,
    # Row 4 (Model 4)
    plots_model4[[1]], plots_model4[[2]], placeholder_plot,
    # Row 5 (Model 5)
    plots_model5[[1]], plots_model5[[2]], placeholder_plot,
    # Row 6 (Model 6)
    plots_model6[[1]], plots_model6[[2]], placeholder_plot
  ),
  ncol = 3,  # Specify 3 columns: Fixed Effects, Random Effects, Correlations
  byrow = TRUE
)

# Save the combined plot
#ggsave("model_visualizations_6x3.png", all_plots, height = 25, width = 15)
all_plots


#########################
# drawing

# brms
get_variables(brm_mo3)

## fixed effects
fixed_effects_samples_brms <- brm_mo3 %>%
  spread_draws(b_logMass, b_Habitat.Densityopen, b_Habitat.DensitysemiMopen)
fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "lightcyan3", 
    color = "lightcyan4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-3.0, 3.0, 1), limits = c(-3.0, 4.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effects"
  ) +
  theme_classic()
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- brm_mo3 %>%
  spread_draws(sd_Phylo__Intercept)
random_effects_samples_brms <- random_effects_samples_brms %>%
  pivot_longer(cols = starts_with("sd_"), 
               names_to = ".variable", 
               values_to = ".value")　%>% 
  mutate(.value = .value^2)

head(random_effects_samples_brms)

brms_p2 <- ggplot(random_effects_samples_brms, aes(x = .value, y = "sd(Intercept)")) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "olivedrab3", 
    color = "olivedrab4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0, 5.0, 2), limits = c(0, 7.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects"
  ) +
  theme_classic() 


## cutpoint
cutpoint_samples_brms <- brm_mo3 %>%
  spread_draws(b_Intercept[cutpoint]) 

brms_p3 <- ggplot(cutpoint_samples_brms, aes(x = b_Intercept, y = factor(cutpoint)))+  
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "pink2", 
    color = "pink4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-5.0, 5.0, 2), limits = c(-5.0, 7.0)) + 
  labs(
    title = "Posterior distributions of cutpoints - brms",
    y = "Cutpoints"
  ) +
  theme_classic()

# brms_plot <- brms_p1 + brms_p2 + brms_p3 + plot_layout(ncol = 1)
# brms_plot


##brms

get_variables(brm_mn3)

## fixed effects
fixed_effects_samples_brms <- brm_mn3 %>%
  spread_draws(b_muInsessorial_Intercept, 
               b_muTerrestrial_Intercept, 
               b_muInsessorial_log_Tail_Length_centered, 
               b_muInsessorial_IsOmnivore, 
               b_muTerrestrial_log_Tail_Length_centered, 
               b_muTerrestrial_IsOmnivore)

fixed_effects_samples_brms <- fixed_effects_samples_brms %>%
  pivot_longer(cols = starts_with("b_"), 
               names_to = ".variable", 
               values_to = ".value")

brms_p1 <- ggplot(fixed_effects_samples_brms, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy", 
    point_interval = "mean_qi", 
    fill = "lightcyan3", 
    color = "lightcyan4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-10.0, 10.0, 2), limits = c(-10.0, 10.0)) + 
  labs(title = "Posterior distributions of fixed effects - brms",
       y = "Fixed effects"
  ) +
  theme_classic()
brms_p1
head(fixed_effects_samples_brms)

## random effects
random_effects_samples_brms <- brm_mn3 %>%
  spread_draws(sd_Phylo__muInsessorial_Intercept, 
               sd_Phylo__muTerrestrial_Intercept
  )

random_effects_samples_brms_long <- random_effects_samples_brms %>%
  pivot_longer(cols = contains("phylo"), 
               names_to = ".variable", 
               values_to = ".value") %>% 
  mutate(.value = .value^2)

head(random_effects_samples_brms_long)

brms_p2 <- ggplot(random_effects_samples_brms_long, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    point_interval = "mean_qi", 
    fill = "olivedrab3", 
    color = "olivedrab4"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(0.0, 80.0, 20), limits = c(0.0, 80.0)) + 
  labs(
    title = "Posterior distributions of random effects - brms",
    y = "Random effects"
  ) +
  theme_classic()

brms_p2

## correlation

correlation_effects_samples_brms <- brm_mn3 %>%
  spread_draws(cor_Phylo__muInsessorial_Intercept__muTerrestrial_Intercept)

correlation_effects_samples_brms_long <- correlation_effects_samples_brms %>%
  pivot_longer(cols = contains("phylo"), 
               names_to = ".variable", 
               values_to = ".value")

head(correlation_effects_samples_brms_long)

brms_p3 <- ggplot(correlation_effects_samples_brms_long, aes(x = .value, y = .variable)) +
  stat_halfeye(
    normalize = "xy",
    fill = "#FF6347", 
    color = "#8B3626"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
  scale_x_continuous(breaks = seq(-1.0, 1.0, 0.2), limits = c(-1.0, 1.0)) + 
  labs(
    title = "Posterior distributions of phylogenetic correlation - brms",
    y = "Correlation"
  ) +
  theme_classic()

brms_p3
# compare
#p_nominal <- mcmcglmm_p1 + mcmcglmm_p2 + mcmcglmm_p3 + brms_p1 + brms_p2 + brms_p3 







