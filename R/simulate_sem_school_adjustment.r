# R/simulate_sem_school_adjustment.R
# ------------------------------------------------------------
# Simulate a psychologically realistic SEM teaching dataset.
#
# Structure:
# - Students nested in classes nested in schools
# - Demographics
# - Observed academic and behavioural outcomes
# - Ordinal Likert indicators
# - Binary correct/incorrect indicators
# - Continuous performance indicators
# - Missing data: item missingness, task missingness, attrition, planned missingness
# - Mild local misfit, cross-loadings, longitudinal residual stability, non-invariance
#
# Dependencies: base R only.
# ------------------------------------------------------------

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

zscore <- function(x) {
  as.numeric(scale(x))
}

clamp <- function(x, lower, upper) {
  pmin(upper, pmax(lower, x))
}

make_ordinal <- function(eta, thresholds) {
  as.integer(cut(
    eta,
    breaks = c(-Inf, thresholds, Inf),
    labels = FALSE,
    right = TRUE
  ))
}

generate_likert <- function(theta,
                            loading = 0.70,
                            intercept = 0,
                            thresholds = c(-1.3, -0.4, 0.4, 1.2),
                            extra = 0,
                            method = 0) {
  n <- length(theta)
  loading_vec <- rep(loading, length.out = n)
  err_sd <- sqrt(pmax(0.30, 1 - loading_vec^2))
  
  eta <- intercept +
    loading_vec * theta +
    extra +
    method +
    rnorm(n, 0, err_sd)
  
  make_ordinal(eta, thresholds)
}

generate_binary_item <- function(theta, difficulty, discrimination = 1.20) {
  p <- inv_logit(discrimination * theta - difficulty)
  rbinom(length(theta), size = 1, prob = p)
}

apply_missing <- function(dat, vars, logit_p) {
  vars <- intersect(vars, names(dat))
  p <- inv_logit(logit_p)
  
  for (v in vars) {
    miss <- rbinom(nrow(dat), size = 1, prob = p) == 1
    miss <- miss & !is.na(dat[[v]])
    dat[[v]][miss] <- NA
  }
  
  dat
}

infer_construct <- function(v) {
  if (grepl("^se[1-5]_", v)) return("academic_self_efficacy")
  if (grepl("^anx[1-6]_", v)) return("school_anxiety")
  if (grepl("^bel[1-5]_", v)) return("school_belonging")
  if (grepl("^mast[1-5]_", v)) return("mastery_goal_orientation")
  if (grepl("^ts[1-5]_", v)) return("teacher_support")
  if (grepl("^abil[0-9]+$", v)) return("academic_reasoning_ability")
  if (grepl("^att_", v)) return("attention_regulation")
  if (grepl("^gpa_|^math_grade|^language_grade|^absences|^study_hours|^screen_time|^homework|^teacher_rating", v)) {
    return("observed_academic_behavioural")
  }
  if (grepl("^missing_|^dropout_|^planned_", v)) return("missing_data_design")
  if (v %in% c(
    "school_id", "class_id", "region", "school_type",
    "class_size", "school_resources", "class_climate_mean"
  )) return("cluster_context")
  if (v %in% c(
    "student_id", "age", "grade", "gender", "ses", "parent_edu",
    "books_home", "language_home", "migration_background",
    "special_education_support"
  )) return("demographics")
  
  "other"
}

infer_scale <- function(v) {
  if (grepl("^(se|anx|bel|mast|ts)[0-9]+_t[123]$", v)) return("ordinal_1_5")
  if (v %in% c("parent_edu", "books_home", "homework_completion")) return("ordinal_1_5")
  if (grepl("^abil[0-9]+$", v)) return("binary_0_1")
  if (grepl("^dropout_|^missing_any_", v)) return("binary_0_1")
  if (v %in% c("migration_background", "special_education_support")) return("binary_0_1")
  if (grepl("^absences_", v)) return("count")
  if (v %in% c(
    "gender", "region", "school_type", "language_home",
    "planned_missing_form", "school_id", "class_id"
  )) return("categorical")
  
  "continuous"
}

infer_wave <- function(v) {
  if (grepl("_t1$", v)) return("t1")
  if (grepl("_t2$", v)) return("t2")
  if (grepl("_t3$", v)) return("t3")
  if (grepl("^abil|^att_|^mast|^ts|math_grade_t1|language_grade_t1", v)) return("t1")
  ""
}

simulate_sem_school_adjustment <- function(seed = 20260529,
                                           n_schools = 30,
                                           classes_per_school = 4,
                                           mean_class_size = 20,
                                           sd_class_size = 4) {
  set.seed(seed)
  
  # ------------------------------------------------------------
  # 1. Schools
  # ------------------------------------------------------------
  
  schools <- data.frame(
    school_id = sprintf("S%02d", seq_len(n_schools)),
    region = sample(
      c("North", "Centre", "South"),
      size = n_schools,
      replace = TRUE,
      prob = c(0.45, 0.25, 0.30)
    ),
    stringsAsFactors = FALSE
  )
  
  region_effect <- ifelse(
    schools$region == "North", 0.25,
    ifelse(schools$region == "Centre", 0.05, -0.25)
  )
  
  schools$school_resources <- zscore(
    rnorm(n_schools, 0, 0.80) + region_effect
  )
  
  schools$school_ses_context <- zscore(
    0.55 * schools$school_resources +
      region_effect +
      rnorm(n_schools, 0, 0.85)
  )
  
  school_type_score <- schools$school_resources + rnorm(n_schools, 0, 0.70)
  
  schools$school_type <- ifelse(
    school_type_score > 0.55, "general",
    ifelse(school_type_score < -0.75, "vocational", "technical")
  )
  
  # ------------------------------------------------------------
  # 2. Classes
  # ------------------------------------------------------------
  
  class_list <- list()
  class_counter <- 1
  
  for (s in seq_len(n_schools)) {
    for (cc in seq_len(classes_per_school)) {
      class_size <- max(12, round(rnorm(1, mean_class_size, sd_class_size)))
      
      class_list[[class_counter]] <- data.frame(
        school_id = schools$school_id[s],
        class_id = sprintf("C%03d", class_counter),
        class_size = class_size,
        grade = sample(9:12, size = 1),
        class_climate_raw = 0.50 * schools$school_resources[s] + rnorm(1, 0, 0.80),
        class_ability_u = rnorm(1),
        class_anxiety_u = rnorm(1),
        class_gpa_u = rnorm(1),
        stringsAsFactors = FALSE
      )
      
      class_counter <- class_counter + 1
    }
  }
  
  classes <- do.call(rbind, class_list)
  classes$class_climate_mean <- zscore(classes$class_climate_raw)
  
  # ------------------------------------------------------------
  # 3. Student rows
  # ------------------------------------------------------------
  
  student_list <- vector("list", nrow(classes))
  student_counter <- 1
  
  for (i in seq_len(nrow(classes))) {
    ids <- seq(from = student_counter, length.out = classes$class_size[i])
    
    student_list[[i]] <- data.frame(
      student_id = ids,
      school_id = classes$school_id[i],
      class_id = classes$class_id[i],
      stringsAsFactors = FALSE
    )
    
    student_counter <- max(ids) + 1
  }
  
  dat <- do.call(rbind, student_list)
  n <- nrow(dat)
  
  class_index <- match(dat$class_id, classes$class_id)
  school_index <- match(dat$school_id, schools$school_id)
  
  dat$region <- schools$region[school_index]
  dat$school_type <- schools$school_type[school_index]
  dat$class_size <- classes$class_size[class_index]
  dat$school_resources <- round(schools$school_resources[school_index], 3)
  dat$class_climate_mean <- round(classes$class_climate_mean[class_index], 3)
  dat$grade <- classes$grade[class_index]
  
  school_ses_context <- schools$school_ses_context[school_index]
  class_ability_u <- classes$class_ability_u[class_index]
  class_anxiety_u <- classes$class_anxiety_u[class_index]
  class_gpa_u <- classes$class_gpa_u[class_index]
  
  # ------------------------------------------------------------
  # 4. Demographics
  # ------------------------------------------------------------
  
  dat$age <- round(dat$grade + 5.5 + rnorm(n, 0, 0.35), 1)
  
  dat$gender <- sample(
    c("girl", "boy", "another_or_no_answer"),
    size = n,
    replace = TRUE,
    prob = c(0.50, 0.48, 0.02)
  )
  
  dat$ses <- round(zscore(
    0.55 * school_ses_context +
      0.20 * dat$class_climate_mean +
      rnorm(n, 0, 0.90)
  ), 3)
  
  dat$parent_edu <- make_ordinal(
    dat$ses + rnorm(n, 0, 0.85),
    thresholds = c(-1.20, -0.40, 0.40, 1.15)
  )
  
  dat$books_home <- make_ordinal(
    dat$ses + rnorm(n, 0, 0.95),
    thresholds = c(-1.10, -0.25, 0.55, 1.30)
  )
  
  p_other_language <- inv_logit(
    -1.70 -
      0.45 * dat$ses +
      0.10 * (dat$region == "North")
  )
  
  dat$language_home <- ifelse(
    rbinom(n, 1, p_other_language) == 1,
    "other_language",
    "majority_language"
  )
  
  p_migration <- inv_logit(
    -2.20 +
      1.10 * (dat$language_home == "other_language") -
      0.25 * dat$ses
  )
  
  dat$migration_background <- rbinom(n, 1, p_migration)
  
  # ------------------------------------------------------------
  # 5. True latent variables
  # ------------------------------------------------------------
  
  theta_ability_t1 <- zscore(
    0.45 * dat$ses +
      0.25 * dat$school_resources +
      0.25 * class_ability_u +
      rnorm(n, 0, 0.85)
  )
  
  theta_attention_t1 <- zscore(
    0.50 * theta_ability_t1 +
      0.25 * dat$ses +
      rnorm(n, 0, 0.90)
  )
  
  theta_teachsup_t1 <- zscore(
    0.65 * dat$class_climate_mean +
      0.25 * dat$school_resources +
      rnorm(n, 0, 0.85)
  )
  
  theta_mastery_t1 <- zscore(
    0.30 * dat$ses +
      0.25 * theta_ability_t1 +
      0.25 * theta_teachsup_t1 +
      rnorm(n, 0, 0.90)
  )
  
  theta_belong_t1 <- zscore(
    0.45 * theta_teachsup_t1 +
      0.25 * dat$class_climate_mean +
      0.15 * dat$ses +
      rnorm(n, 0, 0.85)
  )
  
  theta_selfeff_t1 <- zscore(
    0.35 * theta_ability_t1 +
      0.25 * theta_mastery_t1 +
      0.20 * theta_teachsup_t1 +
      0.15 * dat$ses +
      rnorm(n, 0, 0.80)
  )
  
  theta_anxiety_t1 <- zscore(
    -0.35 * theta_selfeff_t1 -
      0.30 * theta_belong_t1 -
      0.15 * dat$ses +
      0.15 * class_anxiety_u +
      rnorm(n, 0, 0.85)
  )
  
  p_support <- inv_logit(
    -3.10 -
      0.35 * dat$ses -
      0.55 * theta_ability_t1 +
      0.30 * (dat$language_home == "other_language")
  )
  
  dat$special_education_support <- rbinom(n, 1, p_support)
  
  # ------------------------------------------------------------
  # 6. Observed T1 outcomes
  # ------------------------------------------------------------
  
  gpa_latent_t1 <- zscore(
    0.40 * theta_ability_t1 +
      0.25 * theta_selfeff_t1 -
      0.25 * theta_anxiety_t1 +
      0.20 * theta_attention_t1 +
      0.20 * dat$ses +
      0.15 * dat$school_resources +
      0.20 * class_gpa_u +
      rnorm(n, 0, 0.75)
  )
  
  dat$gpa_t1 <- round(clamp(7.0 + 0.90 * gpa_latent_t1, 4.0, 10.0), 1)
  
  math_latent <- zscore(
    0.55 * theta_ability_t1 +
      0.20 * theta_attention_t1 +
      0.15 * theta_selfeff_t1 +
      rnorm(n, 0, 0.85)
  )
  
  language_latent <- zscore(
    0.35 * theta_ability_t1 +
      0.20 * theta_belong_t1 +
      0.20 * theta_selfeff_t1 -
      0.10 * (dat$language_home == "other_language") +
      rnorm(n, 0, 0.90)
  )
  
  dat$math_grade_t1 <- round(clamp(7.0 + 0.85 * math_latent, 4.0, 10.0), 1)
  dat$language_grade_t1 <- round(clamp(7.0 + 0.85 * language_latent, 4.0, 10.0), 1)
  
  mu_abs_t1 <- exp(
    1.35 +
      0.25 * theta_anxiety_t1 -
      0.20 * theta_belong_t1 -
      0.18 * dat$ses -
      0.15 * gpa_latent_t1
  )
  
  dat$absences_t1 <- pmin(rnbinom(n, size = 2.2, mu = mu_abs_t1), 60)
  
  dat$study_hours <- round(clamp(
    6.0 +
      0.85 * theta_selfeff_t1 +
      0.45 * theta_mastery_t1 -
      0.45 * theta_anxiety_t1 +
      0.30 * dat$ses +
      rnorm(n, 0, 2.0),
    0, 25
  ), 1)
  
  dat$screen_time <- round(clamp(
    rlnorm(
      n,
      meanlog = log(2.6) +
        0.12 * theta_anxiety_t1 -
        0.12 * dat$ses,
      sdlog = 0.45
    ),
    0.2, 10
  ), 1)
  
  dat$homework_completion <- make_ordinal(
    0.45 * theta_selfeff_t1 +
      0.25 * theta_mastery_t1 -
      0.25 * theta_anxiety_t1 +
      rnorm(n, 0, 0.90),
    thresholds = c(-1.20, -0.30, 0.45, 1.20)
  )
  
  effort_latent <- zscore(
    0.40 * theta_selfeff_t1 +
      0.30 * theta_mastery_t1 -
      0.20 * theta_anxiety_t1 +
      0.20 * gpa_latent_t1 +
      rnorm(n, 0, 0.85)
  )
  
  dat$teacher_rating_effort <- round(clamp(50 + 8 * effort_latent, 25, 75), 1)
  
  # ------------------------------------------------------------
  # 7. Longitudinal latent variables and outcomes
  # ------------------------------------------------------------
  
  abs_z_t1 <- zscore(dat$absences_t1)
  gpa_z_t1 <- zscore(dat$gpa_t1)
  
  theta_selfeff_t2 <- zscore(
    0.65 * theta_selfeff_t1 +
      0.15 * gpa_z_t1 +
      0.10 * theta_teachsup_t1 +
      0.10 * dat$school_resources +
      rnorm(n, 0, 0.75)
  )
  
  theta_belong_t2 <- zscore(
    0.65 * theta_belong_t1 +
      0.15 * theta_teachsup_t1 +
      0.10 * dat$class_climate_mean -
      0.10 * theta_anxiety_t1 +
      rnorm(n, 0, 0.75)
  )
  
  theta_anxiety_t2 <- zscore(
    0.65 * theta_anxiety_t1 -
      0.15 * theta_selfeff_t1 -
      0.15 * theta_belong_t1 +
      0.12 * abs_z_t1 +
      rnorm(n, 0, 0.75)
  )
  
  gpa_latent_t2 <- zscore(
    0.60 * gpa_z_t1 +
      0.20 * theta_selfeff_t1 -
      0.22 * theta_anxiety_t1 +
      0.15 * theta_ability_t1 +
      0.10 * theta_attention_t1 +
      0.15 * dat$school_resources +
      0.15 * class_gpa_u +
      rnorm(n, 0, 0.80)
  )
  
  dat$gpa_t2 <- round(clamp(7.0 + 0.90 * gpa_latent_t2, 4.0, 10.0), 1)
  
  mu_abs_t2 <- exp(
    1.30 +
      0.25 * theta_anxiety_t2 -
      0.20 * theta_belong_t2 -
      0.15 * dat$ses -
      0.12 * gpa_latent_t2
  )
  
  dat$absences_t2 <- pmin(rnbinom(n, size = 2.1, mu = mu_abs_t2), 70)
  
  abs_z_t2 <- zscore(dat$absences_t2)
  gpa_z_t2 <- zscore(dat$gpa_t2)
  
  theta_selfeff_t3 <- zscore(
    0.65 * theta_selfeff_t2 +
      0.15 * gpa_z_t2 +
      0.10 * theta_teachsup_t1 +
      0.08 * dat$school_resources +
      rnorm(n, 0, 0.75)
  )
  
  theta_belong_t3 <- zscore(
    0.65 * theta_belong_t2 +
      0.15 * theta_teachsup_t1 +
      0.10 * dat$class_climate_mean -
      0.10 * theta_anxiety_t2 +
      rnorm(n, 0, 0.75)
  )
  
  theta_anxiety_t3 <- zscore(
    0.65 * theta_anxiety_t2 -
      0.15 * theta_selfeff_t2 -
      0.15 * theta_belong_t2 +
      0.12 * abs_z_t2 +
      rnorm(n, 0, 0.75)
  )
  
  gpa_latent_t3 <- zscore(
    0.60 * gpa_z_t2 +
      0.20 * theta_selfeff_t2 -
      0.22 * theta_anxiety_t2 +
      0.12 * theta_ability_t1 +
      0.10 * theta_attention_t1 +
      0.12 * dat$school_resources +
      0.15 * class_gpa_u +
      rnorm(n, 0, 0.80)
  )
  
  dat$gpa_t3 <- round(clamp(7.0 + 0.90 * gpa_latent_t3, 4.0, 10.0), 1)
  
  mu_abs_t3 <- exp(
    1.30 +
      0.25 * theta_anxiety_t3 -
      0.20 * theta_belong_t3 -
      0.15 * dat$ses -
      0.12 * gpa_latent_t3
  )
  
  dat$absences_t3 <- pmin(rnbinom(n, size = 2.1, mu = mu_abs_t3), 70)
  
  # ------------------------------------------------------------
  # 8. Ordinal Likert indicators
  # ------------------------------------------------------------
  
  theta_by_wave <- list(
    selfeff = list(t1 = theta_selfeff_t1, t2 = theta_selfeff_t2, t3 = theta_selfeff_t3),
    anxiety = list(t1 = theta_anxiety_t1, t2 = theta_anxiety_t2, t3 = theta_anxiety_t3),
    belong = list(t1 = theta_belong_t1, t2 = theta_belong_t2, t3 = theta_belong_t3)
  )
  
  method_selfeff <- matrix(rnorm(n * 5, 0, 0.25), nrow = n, ncol = 5)
  method_anxiety <- matrix(rnorm(n * 6, 0, 0.25), nrow = n, ncol = 6)
  method_belong <- matrix(rnorm(n * 5, 0, 0.25), nrow = n, ncol = 5)
  
  anx12_content_method <- rnorm(n, 0, 0.22)
  
  se_loadings <- c(0.76, 0.72, 0.70, 0.70, 0.66)
  anx_loadings <- c(0.74, 0.72, 0.70, 0.68, 0.66, 0.63)
  bel_loadings <- c(0.75, 0.72, 0.65, 0.68, 0.66)
  
  for (w in c("t1", "t2", "t3")) {
    for (j in 1:5) {
      loading_j <- se_loadings[j]
      
      if (j == 4) {
        loading_j <- ifelse(
          dat$school_type == "vocational", 0.55,
          ifelse(dat$school_type == "technical", 0.65, 0.75)
        )
      }
      
      dat[[paste0("se", j, "_", w)]] <- generate_likert(
        theta = theta_by_wave$selfeff[[w]],
        loading = loading_j,
        intercept = 0.35,
        thresholds = c(-1.55, -0.65, 0.10, 0.95),
        method = method_selfeff[, j]
      )
    }
    
    for (j in 1:6) {
      extra_j <- 0
      
      if (j == 3) {
        extra_j <- extra_j - 0.25 * (dat$gender == "boy")
      }
      
      if (j == 5) {
        extra_j <- extra_j - 0.22 * theta_by_wave$selfeff[[w]]
      }
      
      if (j %in% c(1, 2)) {
        extra_j <- extra_j + anx12_content_method
      }
      
      dat[[paste0("anx", j, "_", w)]] <- generate_likert(
        theta = theta_by_wave$anxiety[[w]],
        loading = anx_loadings[j],
        intercept = -0.45,
        thresholds = c(-1.00, -0.15, 0.75, 1.60),
        extra = extra_j,
        method = method_anxiety[, j]
      )
    }
    
    for (j in 1:5) {
      extra_j <- 0
      
      if (j == 2 && w == "t2") extra_j <- extra_j + 0.15
      if (j == 2 && w == "t3") extra_j <- extra_j + 0.30
      
      if (j == 3) {
        extra_j <- extra_j + 0.25 * theta_teachsup_t1
      }
      
      dat[[paste0("bel", j, "_", w)]] <- generate_likert(
        theta = theta_by_wave$belong[[w]],
        loading = bel_loadings[j],
        intercept = 0.30,
        thresholds = c(-1.40, -0.50, 0.25, 1.10),
        extra = extra_j,
        method = method_belong[, j]
      )
    }
  }
  
  for (j in 1:5) {
    dat[[paste0("mast", j, "_t1")]] <- generate_likert(
      theta = theta_mastery_t1,
      loading = c(0.75, 0.72, 0.70, 0.66, 0.64)[j],
      intercept = 0.35,
      thresholds = c(-1.45, -0.55, 0.20, 1.05),
      method = rnorm(n, 0, 0.10)
    )
  }
  
  for (j in 1:5) {
    dat[[paste0("ts", j, "_t1")]] <- generate_likert(
      theta = theta_teachsup_t1,
      loading = c(0.78, 0.74, 0.70, 0.68, 0.65)[j],
      intercept = 0.25,
      thresholds = c(-1.35, -0.45, 0.30, 1.10),
      method = rnorm(n, 0, 0.12)
    )
  }
  
  # ------------------------------------------------------------
  # 9. Binary academic reasoning indicators
  # ------------------------------------------------------------
  
  item_difficulty <- c(
    -1.80, -1.45, -1.10, -0.75, -0.45, -0.15,
    0.15, 0.45, 0.75, 1.05, 1.35, 1.65
  )
  
  item_discrimination <- c(
    1.05, 1.05, 1.10, 1.15, 1.20, 1.20,
    1.25, 1.25, 1.20, 1.15, 1.10, 1.05
  )
  
  for (j in 1:12) {
    dat[[paste0("abil", j)]] <- generate_binary_item(
      theta = theta_ability_t1,
      difficulty = item_difficulty[j],
      discrimination = item_discrimination[j]
    )
  }
  
  # ------------------------------------------------------------
  # 10. Continuous attention indicators
  # ------------------------------------------------------------
  
  dat$att_acc <- round(clamp(
    inv_logit(1.20 + 0.65 * theta_attention_t1 + rnorm(n, 0, 0.55)),
    0.35, 0.99
  ), 3)
  
  dat$att_rt_inv <- round(zscore(
    0.75 * theta_attention_t1 + rnorm(n, 0, 0.70)
  ), 3)
  
  dat$att_inhibit <- round(clamp(
    50 + 8.5 * theta_attention_t1 + rnorm(n, 0, 7.0),
    20, 80
  ), 1)
  
  dat$att_switch <- round(clamp(
    50 + 8.0 * theta_attention_t1 + rnorm(n, 0, 7.5),
    20, 80
  ), 1)
  
  # ------------------------------------------------------------
  # 11. Missingness
  # ------------------------------------------------------------
  
  vars_t1 <- c(
    "gpa_t1", "math_grade_t1", "language_grade_t1", "absences_t1",
    "study_hours", "screen_time", "homework_completion", "teacher_rating_effort",
    paste0("se", 1:5, "_t1"),
    paste0("anx", 1:6, "_t1"),
    paste0("bel", 1:5, "_t1"),
    paste0("mast", 1:5, "_t1"),
    paste0("ts", 1:5, "_t1"),
    paste0("abil", 1:12),
    "att_acc", "att_rt_inv", "att_inhibit", "att_switch"
  )
  
  vars_t2 <- c(
    "gpa_t2", "absences_t2",
    paste0("se", 1:5, "_t2"),
    paste0("anx", 1:6, "_t2"),
    paste0("bel", 1:5, "_t2")
  )
  
  vars_t3 <- c(
    "gpa_t3", "absences_t3",
    paste0("se", 1:5, "_t3"),
    paste0("anx", 1:6, "_t3"),
    paste0("bel", 1:5, "_t3")
  )
  
  ordinal_item_vars <- grep(
    "^(se|anx|bel|mast|ts)[0-9]+_t[123]$",
    names(dat),
    value = TRUE
  )
  
  attention_vars <- c("att_acc", "att_rt_inv", "att_inhibit", "att_switch")
  
  dat$planned_missing_form <- sample(
    c("A", "B", "C"),
    size = n,
    replace = TRUE,
    prob = c(0.34, 0.33, 0.33)
  )
  
  form_b_vars <- c(
    paste0("anx", 4:6, "_t2"),
    paste0("bel", 4:5, "_t3")
  )
  
  form_c_vars <- c(
    paste0("se", 4:5, "_t2"),
    paste0("anx", 4:6, "_t3")
  )
  
  dat[dat$planned_missing_form == "B", form_b_vars] <- NA
  dat[dat$planned_missing_form == "C", form_c_vars] <- NA
  
  p_drop_t2 <- inv_logit(
    -2.25 +
      0.40 * theta_anxiety_t1 -
      0.35 * gpa_latent_t1 -
      0.30 * dat$ses +
      0.25 * zscore(dat$absences_t1)
  )
  
  dat$dropout_t2 <- rbinom(n, 1, p_drop_t2)
  
  p_drop_t3 <- inv_logit(
    -1.75 +
      0.45 * theta_anxiety_t1 -
      0.40 * gpa_latent_t1 -
      0.30 * dat$ses +
      0.25 * zscore(dat$absences_t1)
  )
  
  dropout_t3_new <- rbinom(n, 1, p_drop_t3)
  dat$dropout_t3 <- ifelse(dat$dropout_t2 == 1, 1, dropout_t3_new)
  
  dat[dat$dropout_t2 == 1, vars_t2] <- NA
  dat[dat$dropout_t3 == 1, vars_t3] <- NA
  
  item_missing_logit <- -3.10 +
    0.25 * theta_anxiety_t1 -
    0.20 * dat$ses +
    0.15 * zscore(dat$absences_t1)
  
  dat <- apply_missing(dat, ordinal_item_vars, item_missing_logit)
  
  ability_missing_logit <- -3.35 -
    0.15 * dat$ses +
    0.12 * theta_anxiety_t1
  
  dat <- apply_missing(dat, paste0("abil", 1:12), ability_missing_logit)
  
  task_missing_logit <- -2.10 +
    0.30 * zscore(dat$absences_t1) -
    0.20 * dat$ses
  
  task_missing <- rbinom(n, 1, inv_logit(task_missing_logit)) == 1
  dat[task_missing, attention_vars] <- NA
  
  demo_missing_logit <- -4.00 +
    0.15 * theta_anxiety_t1 -
    0.20 * dat$ses
  
  dat <- apply_missing(dat, c("parent_edu", "books_home"), demo_missing_logit)
  
  dat$missing_any_t1 <- as.integer(rowSums(is.na(dat[, intersect(vars_t1, names(dat))])) > 0)
  dat$missing_any_t2 <- as.integer(rowSums(is.na(dat[, intersect(vars_t2, names(dat))])) > 0)
  dat$missing_any_t3 <- as.integer(rowSums(is.na(dat[, intersect(vars_t3, names(dat))])) > 0)
  
  # ------------------------------------------------------------
  # 12. Reorder columns
  # ------------------------------------------------------------
  
  id_vars <- c(
    "student_id", "school_id", "class_id", "region", "school_type",
    "class_size", "school_resources", "class_climate_mean"
  )
  
  demo_vars <- c(
    "age", "grade", "gender", "ses", "parent_edu", "books_home",
    "language_home", "migration_background", "special_education_support"
  )
  
  observed_vars <- c(
    "gpa_t1", "gpa_t2", "gpa_t3",
    "math_grade_t1", "language_grade_t1",
    "absences_t1", "absences_t2", "absences_t3",
    "study_hours", "screen_time", "homework_completion",
    "teacher_rating_effort"
  )
  
  likert_vars <- c(
    paste0("se", 1:5, "_t1"), paste0("se", 1:5, "_t2"), paste0("se", 1:5, "_t3"),
    paste0("anx", 1:6, "_t1"), paste0("anx", 1:6, "_t2"), paste0("anx", 1:6, "_t3"),
    paste0("bel", 1:5, "_t1"), paste0("bel", 1:5, "_t2"), paste0("bel", 1:5, "_t3"),
    paste0("mast", 1:5, "_t1"),
    paste0("ts", 1:5, "_t1")
  )
  
  binary_vars <- paste0("abil", 1:12)
  
  continuous_indicator_vars <- attention_vars
  
  missing_vars <- c(
    "missing_any_t1", "missing_any_t2", "missing_any_t3",
    "dropout_t2", "dropout_t3", "planned_missing_form"
  )
  
  ordered_cols <- c(
    id_vars,
    demo_vars,
    observed_vars,
    likert_vars,
    binary_vars,
    continuous_indicator_vars,
    missing_vars
  )
  
  student_data <- dat[, ordered_cols]
  
  # ------------------------------------------------------------
  # 13. Instructor truth file
  # ------------------------------------------------------------
  
  true_scores <- data.frame(
    student_id = dat$student_id,
    school_id = dat$school_id,
    class_id = dat$class_id,
    
    theta_selfeff_t1 = round(theta_selfeff_t1, 4),
    theta_selfeff_t2 = round(theta_selfeff_t2, 4),
    theta_selfeff_t3 = round(theta_selfeff_t3, 4),
    
    theta_anxiety_t1 = round(theta_anxiety_t1, 4),
    theta_anxiety_t2 = round(theta_anxiety_t2, 4),
    theta_anxiety_t3 = round(theta_anxiety_t3, 4),
    
    theta_belong_t1 = round(theta_belong_t1, 4),
    theta_belong_t2 = round(theta_belong_t2, 4),
    theta_belong_t3 = round(theta_belong_t3, 4),
    
    theta_mastery_t1 = round(theta_mastery_t1, 4),
    theta_teachsup_t1 = round(theta_teachsup_t1, 4),
    theta_ability_t1 = round(theta_ability_t1, 4),
    theta_attention_t1 = round(theta_attention_t1, 4),
    
    gpa_latent_t1 = round(gpa_latent_t1, 4),
    gpa_latent_t2 = round(gpa_latent_t2, 4),
    gpa_latent_t3 = round(gpa_latent_t3, 4),
    
    stringsAsFactors = FALSE
  )
  
  parameters <- list(
    seed = seed,
    n_students = nrow(student_data),
    n_schools = n_schools,
    n_classes = nrow(classes),
    classes_per_school = classes_per_school,
    
    item_difficulty = item_difficulty,
    item_discrimination = item_discrimination,
    
    likert_thresholds = list(
      selfeff = c(-1.55, -0.65, 0.10, 0.95),
      anxiety = c(-1.00, -0.15, 0.75, 1.60),
      belonging = c(-1.40, -0.50, 0.25, 1.10),
      mastery = c(-1.45, -0.55, 0.20, 1.05),
      teacher_support = c(-1.35, -0.45, 0.30, 1.10)
    ),
    
    intentional_imperfections = c(
      "anx5 weakly cross-loads on low self-efficacy",
      "bel3 weakly cross-loads on teacher support",
      "anx1 and anx2 share a content residual",
      "same items share person-specific residual components over time",
      "se4 has a lower loading in vocational schools",
      "anx3 has mild gender-related threshold/DIF",
      "bel2 becomes easier to endorse over time",
      "attrition is MAR-like and related to anxiety, GPA, SES, and absences",
      "attention indicators have whole-block task missingness"
    )
  )
  
  codebook <- data.frame(
    variable = names(student_data),
    type = vapply(student_data, function(x) class(x)[1], character(1)),
    scale = vapply(names(student_data), infer_scale, character(1)),
    wave = vapply(names(student_data), infer_wave, character(1)),
    construct = vapply(names(student_data), infer_construct, character(1)),
    missing_rate = round(vapply(student_data, function(x) mean(is.na(x)), numeric(1)), 3),
    stringsAsFactors = FALSE
  )
  
  list(
    data = student_data,
    codebook = codebook,
    truth = list(
      true_scores = true_scores,
      parameters = parameters
    )
  )
}

write_sem_school_adjustment <- function(output_dir = "data/processed",
                                        seed = 20260529,
                                        n_schools = 30,
                                        classes_per_school = 4,
                                        mean_class_size = 20,
                                        sd_class_size = 4) {
  sim <- simulate_sem_school_adjustment(
    seed = seed,
    n_schools = n_schools,
    classes_per_school = classes_per_school,
    mean_class_size = mean_class_size,
    sd_class_size = sd_class_size
  )
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  data_path <- file.path(output_dir, "sem_school_adjustment.csv")
  codebook_path <- file.path(output_dir, "sem_school_adjustment_codebook.csv")
  truth_path <- file.path(output_dir, "sem_school_adjustment_truth.rds")
  
  write.csv(sim$data, data_path, row.names = FALSE, na = "")
  write.csv(sim$codebook, codebook_path, row.names = FALSE, na = "")
  saveRDS(sim$truth, truth_path)
  
  invisible(list(
    data_path = data_path,
    codebook_path = codebook_path,
    truth_path = truth_path,
    summary = list(
      n_students = nrow(sim$data),
      n_variables = ncol(sim$data),
      n_schools = length(unique(sim$data$school_id)),
      n_classes = length(unique(sim$data$class_id)),
      dropout_t2_rate = mean(sim$data$dropout_t2),
      dropout_t3_rate = mean(sim$data$dropout_t3),
      mean_item_missing = mean(sim$codebook$missing_rate[sim$codebook$scale %in% c("ordinal_1_5", "binary_0_1")])
    )
  ))
}

# Run when called as a script, e.g. Rscript R/simulate_sem_school_adjustment.R
if (sys.nframe() == 0) {
  out <- write_sem_school_adjustment()
  print(out$summary)
  message("Wrote: ", out$data_path)
  message("Wrote: ", out$codebook_path)
  message("Wrote: ", out$truth_path)
}