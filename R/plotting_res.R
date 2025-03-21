# figs

library(tidyverse)
library(tidybayes)
library(ggplot2)
library(patchwork)
library(here)
# Load models
mod1 <- readRDS(here("Rdata", "mod1.rds"))
mod2 <- readRDS(here("Rdata", "mod2.rds"))
mod3 <- readRDS(here("Rdata", "mod3.rds"))


# Function to get variable names dynamically
get_variables_dynamic <- function(model, pattern) {
  variables <- get_variables(model)
  variables[grep(pattern, variables)]
}

# to change

rename_vars <- function(variable) {
  variable <- gsub("b_sigma_intercept", "b_s", variable)            
  variable <- gsub("b_intercept", "b_l", variable)        
  variable <- gsub("cor_", "r_", variable)  
  variable <- gsub("_", " ", variable)    
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
  visualize_fixed_effects(mod1, "Model 1: Fixed Effects"),
  visualize_random_effects(mod1, "Model 1: Random Effects")
)

plots_model2 <- list(
  visualize_fixed_effects(mod2, "Model 2: Fixed Effects"),
  visualize_random_effects(mod2, "Model 2: Random Effects"),
  visualize_correlations(mod2, "Model 2: Residual Correlations")
)

plots_model3 <- list(
  visualize_fixed_effects(mod3, "Model 3: Fixed Effects"),
  visualize_random_effects(mod3, "Model 3: Random Effects"),
  visualize_correlations(mod3, "Model 3: Correlations")
)


# Create a placeholder plot for models without correlations
placeholder_plot <- ggplot() +
  theme_void() +
  labs(title = "No Correlations Available")

# Arrange plots in 6 rows (one for each model) and 3 columns (Fixed, Random, Correlations)
all_plots <- patchwork::wrap_plots(
  list(
    # Row 1 (Model 1)
    plots_model1[[1]], plots_model1[[2]], placeholder_plot,
    # Row 2 (Model 2)
    plots_model2[[1]], plots_model2[[2]], plots_model2[[3]],
    # Row 3 (Model 3)
    plots_model3[[1]], plots_model3[[2]], plots_model3[[3]]),
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







