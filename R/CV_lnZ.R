# --- CV vs SD[log(z)] simulation demo ---------------------------------------
set.seed(123)

# helper: theoretical SD[log(z)] from CV (holds exactly for log-normal)
sdlog_from_cv <- function(cv) sqrt(log(1 + cv^2))

# grid of CVs to test
cv_grid <- seq(0.02, 1.0, by = 0.02)
n_per_cv <- 10000L

# ---- (A) Log-normal: theory vs simulation ----------------------------------
sim_lnorm <- lapply(cv_grid, function(cv) {
  sdlog <- sdlog_from_cv(cv)        # exact for log-normal
  mu_log <- 0                       # arbitrary; CV doesn't depend on mu
  x <- rlnorm(n_per_cv, meanlog = mu_log, sdlog = sdlog)
  data.frame(
    dist = "lognormal",
    cv_target = cv,
    cv_sample = sd(x) / mean(x),
    sdlog_sample = sd(log(x)),
    sdlog_theory = sdlog
  )
})
sim_lnorm <- do.call(rbind, sim_lnorm)

# ---- (B) Gamma: small-CV approximation test --------------------------------
# For Gamma(shape=k, scale=theta), CV = 1/sqrt(k) (independent of scale).
# We'll set theta = 1 and choose shape to hit a target CV.
sim_gamma <- lapply(cv_grid, function(cv) {
  shape <- 1 / (cv^2)
  scale <- 1
  x <- rgamma(n_per_cv, shape = shape, scale = scale)
  data.frame(
    dist = "gamma",
    cv_target = cv,
    cv_sample = sd(x) / mean(x),
    sdlog_sample = sd(log(x)),
    sdlog_theory = sdlog_from_cv(cv) # not exact; approximation to compare
  )
})
sim_gamma <- do.call(rbind, sim_gamma)

sim_all <- rbind(sim_lnorm, sim_gamma)

# error metrics
sim_all$abs_err_theory_vs_sample_sdlog <- with(sim_all, abs(sdlog_theory - sdlog_sample))
sim_all$rel_err_theory_vs_sample_sdlog <- with(sim_all, abs(sdlog_theory - sdlog_sample) / sdlog_theory)
sim_all$ratio_sdlog_to_cv <- with(sim_all, sdlog_sample / cv_sample)

# ---- Quick summary table at CV = {0.10, 0.30, 0.50} ------------------------
pick <- c(0.10, 0.30, 0.50)
summary_tab <- subset(sim_all, round(cv_target, 2) %in% pick, 
                      select = c(dist, cv_target, cv_sample, sdlog_sample, sdlog_theory, 
                                 abs_err_theory_vs_sample_sdlog, rel_err_theory_vs_sample_sdlog))
summary_tab <- summary_tab[order(summary_tab$dist, summary_tab$cv_target), ]
print(summary_tab, row.names = FALSE, digits = 4)

# ---- Plots ------------------------------------------------------------------
# 1) SD[log(z)] vs CV: theory line and simulation points
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  # Theory line (same for both panels)
  theory_df <- data.frame(
    cv = cv_grid,
    sdlog = sdlog_from_cv(cv_grid)
  )
  
  p1 <- ggplot(sim_all, aes(x = cv_sample, y = sdlog_sample, colour = dist)) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    geom_line(data = theory_df, aes(x = cv, y = sdlog), inherit.aes = FALSE) +
    geom_point(alpha = 0.4, size = 1) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Sample CV (original scale)",
         y = "Sample SD of log(z)",
         colour = "Distribution",
         title = "SD[log(z)] vs CV: simulation (points) and theory for log-normal (solid line)",
         subtitle = "Dashed line: y = x (small-CV approximation SD[log(z)] ≈ CV)")
  
  # 2) Approximation quality: SD[log]/CV
  p2 <- ggplot(sim_all, aes(x = cv_sample, y = ratio_sdlog_to_cv, colour = dist)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(alpha = 0.5, size = 1) +
    labs(x = "Sample CV", y = "SD[log(z)] / CV",
         colour = "Distribution",
         title = "How close is SD[log(z)] to CV?") +
    ylim(0, 1.6)
  
  print(p1)
  print(p2)
} else {
  message("Install ggplot2 for plots; the numeric summary has been printed above.")
}

# ---- Takeaways printed to console ------------------------------------------
cat("\nKey takeaways:\n",
    "- Log-normal: sd(log(z)) matches sqrt(log(1+CV^2)) exactly in theory; simulations should align closely.\n",
    "- Gamma: sd(log(z)) ≈ CV when CV is small (e.g., ≤ 0.3), but departs as CV grows, as expected.\n")