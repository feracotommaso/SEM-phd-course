library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)

library(lavaan)
library(modsem)

# -----------------------------
# 1) Set simulation parameters
# -----------------------------
R     <- 500   # number of replications
N     <- 400   # sample size

rho_xz <- 0.30

beta_x  <- 0.30
beta_z  <- 0.30
beta_xz <- 0.20  # <-- ground truth interaction (what we'll compare against)

lambda_int <- c(.80, .70, .75)
lambda_anx <- c(.80, .70, .75)
lambda_ach <- c(.80, .70, .75)

# -----------------------------
# 2) Model syntax (lavaan/modsem)
# -----------------------------
model_int <- '
  intelligence =~ int1 + int2 + int3
  anxiety      =~ anx1 + anx2 + anx3
  achievement  =~ ach1 + ach2 + ach3

  achievement ~ intelligence + anxiety +
                intelligence:anxiety
'

# -----------------------------
# 3) One replication: simulate data explicitly (rnorm) + indicators
# -----------------------------
simulate_one <- function(seed) {
  set.seed(seed)
  
  # Latent predictors with correlation rho_xz
  e1 <- rnorm(N)
  e2 <- rnorm(N)
  intelligence <- e1
  anxiety      <- rho_xz * e1 + sqrt(1 - rho_xz^2) * e2
  
  int_x_anx_lat <- intelligence * anxiety
  
  # Latent outcome with chosen residual variance (approx var(Y) ~ 1)
  linpart <- beta_x * intelligence + beta_z * anxiety + beta_xz * int_x_anx_lat
  var_e   <- max(1e-6, 1 - var(linpart))
  achievement <- linpart + rnorm(N, 0, sqrt(var_e))
  
  # Helper: indicators with approx unit variance
  make_inds <- function(f, lambda, prefix) {
    theta <- 1 - lambda^2
    out <- sapply(seq_along(lambda), function(j) {
      lambda[j] * f + rnorm(N, 0, sqrt(theta[j]))
    })
    out <- as.data.frame(out)
    names(out) <- paste0(prefix, 1:3)
    out
  }
  
  dat <- bind_cols(
    make_inds(intelligence, lambda_int, "int"),
    make_inds(anxiety,      lambda_anx, "anx"),
    make_inds(achievement,  lambda_ach, "ach")
  )
  
  dat
}

# -----------------------------
# 4) Fit models + extract interaction term each replication
# -----------------------------

# Extract interaction estimate from modsem objects (PI / LMS / QML)
extract_modsem_int <- function(fit) {
  pe <- modsem::parameter_estimates(fit)
  
  se_col <- intersect(c("std.error", "se"), names(pe))[1]
  p_col  <- intersect(c("p.value", "pvalue"), names(pe))[1]
  
  row <- pe %>%
    filter(op == "~", lhs == "achievement") %>%
    filter(grepl("intelligence", rhs) & grepl("anxiety", rhs)) %>%
    slice(1)
  
  if (nrow(row) == 0) {
    return(tibble(est = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  tibble(
    est = row$est,
    se  = if (!is.na(se_col)) row[[se_col]] else NA_real_,
    p   = if (!is.na(p_col))  row[[p_col]]  else NA_real_
  )
}

# Extract interaction from observed sum-score regression (lm)
extract_lm_int <- function(dat) {
  dat2 <- dat %>%
    mutate(
      intelligence_total = rowMeans(across(c(int1, int2, int3))),
      anxiety_total      = rowMeans(across(c(anx1, anx2, anx3))),
      achievement_total  = rowMeans(across(c(ach1, ach2, ach3))),
      int_c = as.numeric(scale(intelligence_total, scale = FALSE)),
      anx_c = as.numeric(scale(anxiety_total,      scale = FALSE)),
      int_x_anx = int_c * anx_c
    )
  
  m <- lm(achievement_total ~ int_c + anx_c + int_x_anx, data = dat2)
  co <- coef(summary(m))
  
  tibble(
    est = unname(co["int_x_anx", "Estimate"]),
    se  = unname(co["int_x_anx", "Std. Error"]),
    p   = unname(co["int_x_anx", "Pr(>|t|)"])
  )
}

# Extract interaction from observed sum-score lavaan::sem
extract_lavaan_obs_int <- function(dat) {
  dat2 <- dat %>%
    mutate(
      intelligence_total = rowMeans(across(c(int1, int2, int3))),
      anxiety_total      = rowMeans(across(c(anx1, anx2, anx3))),
      achievement_total  = rowMeans(across(c(ach1, ach2, ach3))),
      int_c = as.numeric(scale(intelligence_total, scale = FALSE)),
      anx_c = as.numeric(scale(anxiety_total,      scale = FALSE)),
      int_x_anx = int_c * anx_c
    )
  
  mod_obs <- '
    achievement_total ~ int_c + anx_c + int_x_anx
  '
  
  fit <- sem(mod_obs, data = dat2, meanstructure = TRUE)
  
  pe <- parameterEstimates(fit)
  row <- pe %>% filter(op == "~", lhs == "achievement_total", rhs == "int_x_anx") %>% slice(1)
  
  tibble(
    est = if (nrow(row) == 1) row$est else NA_real_,
    se  = if (nrow(row) == 1) row$se  else NA_real_,
    p   = if (nrow(row) == 1) row$pvalue else NA_real_
  )
}


# ---------------------------------------------------------------------------

extract_modsem_int_std <- function(fit) {
  pe <- modsem::standardized_estimates(fit)
  
  # Same structure as parameter_estimates(); 'est' is now standardized
  se_col <- intersect(c("std.error", "se"), names(pe))[1]
  p_col  <- intersect(c("p.value", "pvalue"), names(pe))[1]
  
  row <- pe |>
    dplyr::filter(op == "~", lhs == "achievement") |>
    dplyr::filter(grepl("intelligence", rhs) & grepl("anxiety", rhs)) |>
    dplyr::slice(1)
  
  if (nrow(row) == 0) {
    return(tibble::tibble(est = NA_real_, se = NA_real_, p = NA_real_))
  }
  
  tibble::tibble(
    est = row$est,  # standardized estimate
    se  = if (!is.na(se_col)) row[[se_col]] else NA_real_,
    p   = if (!is.na(p_col))  row[[p_col]]  else NA_real_
  )
}



extract_lm_int_std <- function(dat) {
  dat2 <- dat |>
    dplyr::mutate(
      intelligence_total = rowMeans(dplyr::across(c(int1, int2, int3))),
      anxiety_total      = rowMeans(dplyr::across(c(anx1, anx2, anx3))),
      achievement_total  = rowMeans(dplyr::across(c(ach1, ach2, ach3)))
    ) |>
    dplyr::mutate(
      int_z = as.numeric(scale(intelligence_total)),
      anx_z = as.numeric(scale(anxiety_total)),
      ach_z = as.numeric(scale(achievement_total)),
      int_x_anx_z = int_z * anx_z
    )
  
  m <- lm(ach_z ~ int_z + anx_z + int_x_anx_z, data = dat2)
  co <- coef(summary(m))
  
  tibble::tibble(
    est = unname(co["int_x_anx_z", "Estimate"]),
    se  = unname(co["int_x_anx_z", "Std. Error"]),
    p   = unname(co["int_x_anx_z", "Pr(>|t|)"])
  )
}


extract_lavaan_obs_int_std <- function(dat) {
  dat2 <- dat |>
    dplyr::mutate(
      intelligence_total = rowMeans(dplyr::across(c(int1, int2, int3))),
      anxiety_total      = rowMeans(dplyr::across(c(anx1, anx2, anx3))),
      achievement_total  = rowMeans(dplyr::across(c(ach1, ach2, ach3))),
      int_z = as.numeric(scale(intelligence_total)),
      anx_z = as.numeric(scale(anxiety_total)),
      ach_z = as.numeric(scale(achievement_total)),
      int_x_anx_z = int_z * anx_z
    )
  
  mod_obs <- "ach_z ~ int_z + anx_z + int_x_anx_z"
  fit <- lavaan::sem(mod_obs, data = dat2, meanstructure = TRUE)
  
  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)
  row <- pe |>
    dplyr::filter(op == "~", lhs == "ach_z", rhs == "int_x_anx_z") |>
    dplyr::slice(1)
  
  tibble::tibble(
    est = row$std.all,   # standardized coefficient
    se  = row$se,
    p   = row$pvalue
  )
}

fit_all_models <- function(dat) {
  
  # Observed models
  out_lm <- tryCatch(
    extract_lm_int_std(dat) %>% mutate(model = "Observed sum-scores (lm)", ok = TRUE),
    error = function(e) tibble(est = NA_real_, se = NA_real_, p = NA_real_, model = "Observed sum-scores (lm)", ok = FALSE)
  )
  
  out_obs_sem <- tryCatch(
    extract_lavaan_obs_int_std(dat) %>% mutate(model = "Observed sum-scores (lavaan)", ok = TRUE),
    error = function(e) tibble(est = NA_real_, se = NA_real_, p = NA_real_, model = "Observed sum-scores (lavaan)", ok = FALSE)
  )
  
  # Latent models (modsem)
  out_pi <- tryCatch(
    extract_modsem_int_std(modsem(model_int, data = dat, method = "dblcent")) %>%
      mutate(model = "Latent PI dblcent", ok = TRUE),
    error = function(e) tibble(est = NA_real_, se = NA_real_, p = NA_real_, model = "Latent PI dblcent", ok = FALSE)
  )
  
  out_lms <- tryCatch(
    extract_modsem_int_std(modsem(model_int, data = dat, method = "lms")) %>%
      mutate(model = "Latent LMS", ok = TRUE),
    error = function(e) tibble(est = NA_real_, se = NA_real_, p = NA_real_, model = "Latent LMS", ok = FALSE)
  )
  
  out_qml <- tryCatch(
    extract_modsem_int_std(modsem(model_int, data = dat, method = "qml")) %>%
      mutate(model = "Latent QML", ok = TRUE),
    error = function(e) tibble(est = NA_real_, se = NA_real_, p = NA_real_, model = "Latent QML", ok = FALSE)
  )
  
  bind_rows(out_lm, out_obs_sem, out_pi, out_lms, out_qml)
}

# -----------------------------
# 5) Run the simulation
# -----------------------------
results <- map_dfr(1:R, function(r) {
  dat_r <- simulate_one(seed = 1000 + r)
  fit_all_models(dat_r) %>%
    mutate(rep = r, beta_true = beta_xz)
})

# -----------------------------
# 6) Summaries: bias, RMSE, convergence, "power" (p<.05)
# -----------------------------
summ <- results %>%
  group_by(model) %>%
  summarise(
    n_total = n(),
    n_ok    = sum(ok, na.rm = TRUE),
    conv_rate = mean(ok, na.rm = TRUE),
    
    mean_est = mean(est, na.rm = TRUE),
    sd_est   = sd(est, na.rm = TRUE),
    bias     = mean(est - beta_true, na.rm = TRUE),
    rmse     = sqrt(mean((est - beta_true)^2, na.rm = TRUE)),
    
    mean_se  = mean(se, na.rm = TRUE),
    power_p  = mean(p < .05, na.rm = TRUE),
    .groups = "drop"
  )

summ

# -----------------------------
# 7) Plot: sampling distribution of interaction estimates vs truth
# -----------------------------
ggplot(results %>% filter(!is.na(est)),
       aes(x = est)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = beta_xz, linetype = 2, linewidth = 0.8) +
  facet_wrap(~ model, scales = "free_y") +
  labs(
    title = "Sampling distribution of the interaction estimate across methods",
    subtitle = paste0("Dashed line = ground truth beta_xz = ", beta_xz),
    x = "Estimated interaction (b_XZ)",
    y = "Count"
  )


