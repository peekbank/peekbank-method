# source me for all analyses

library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(tibble)
library(agreement)
library(multidplyr)
library(boot)
# because we don't have full on tidyverse on the cluster, we specify the parts we actually need

# Seed for random number generation
set.seed(42)


# key ICC function
get_icc <- function(x, column = "accuracy", object = "stimulus") {
  if (object == "stimulus") {
    iccs <- dim_icc(x,
      model = "2A",
      type = "consistency",
      unit = "average",
      object = target_label,
      rater = administration_id,
      trial = repetition,
      score = {{ column }},
      bootstrap = 0
    )
  } else {
    iccs <- dim_icc(x,
      model = "2A",
      type = "consistency",
      unit = "average",
      object = administration_id,
      rater = target_label,
      trial = repetition,
      score = {{ column }},
      bootstrap = 0
    )
  }

  return(iccs$Inter_ICC)
}
bootstrap_icc <- function(x, column = "accuracy", bootstrap = 2000) {
  iccs <- dim_icc(x,
    model = "2A",
    type = "consistency",
    unit = "average",
    object = administration_id,
    rater = target_label,
    trial = repetition,
    score = {{ column }},
    bootstrap = bootstrap
  )
  # View(iccs$boot_results[["t"]])
  t <- iccs$boot_results[["t"]][, 6]
  t_no_na <- t[!is.na(t)]
  # print(t_no_na)
  # View(iccs$boot_results[["t0"]])
  t0 <- iccs$boot_results[["t0"]][6][1]
  # print(t0)
  if (is.na(t0)) {
    return(tibble(est = NA, lower = NA, upper = NA))
  }
  if (length(t_no_na) == 0) {
    return(tibble(est = iccs$Inter_ICC, lower = NA, upper = NA))
  }
  iccs$boot_results$R <- length(t_no_na)
  ci <- boot::boot.ci(boot.out = iccs$boot_results, t = t_no_na, t0 = t0, type = "basic")
  return(tibble(est = iccs$Inter_ICC, lower = ci$basic[4], upper = ci$basic[5]))
}

do_cdi <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    left_join(cdi_data) |>
    ungroup() |>
    summarise(
      cor_comp = ifelse(sum(!is.na(comp) & !is.na(mean_var)) > 2, cor.test(mean_var, comp)$estimate, as.numeric(NA)),
      cor_prod = ifelse(sum(!is.na(prod) & !is.na(mean_var)) > 2, cor.test(mean_var, prod)$estimate, as.numeric(NA)),
      cor_age = ifelse(sum(!is.na(age) & !is.na(mean_var)) > 2, cor.test(mean_var, age)$estimate, as.numeric(NA))
    )
  cors <- c(summ$cor_comp[1], summ$cor_prod[1], summ$cor_age[1])
  names(cors) <- c("cor_comp", "cor_prod", "cor_age")

  return(cors)
}

boot_cdi <- function(data) {
  data |>
    group_by(dataset_name) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi, 2000)
      if (is.na(b$t0[1])) {
        comp_lower <- NA
        comp_upper <- NA
      } else {
        ci_comp <- boot::boot.ci(b, index = 1, type = "basic")
        comp_lower <- ci_comp$basic[4]
        comp_upper <- ci_comp$basic[5]
      }
      if (is.na(b$t0[2])) {
        prod_lower <- NA
        prod_upper <- NA
      } else {
        ci_prod <- boot::boot.ci(b, index = 2, type = "basic")
        prod_lower <- ci_prod$basic[4]
        prod_upper <- ci_prod$basic[5]
      }
      if (is.na(b$t0[3])) {
        age_lower <- NA
        age_upper <- NA
      } else {
        ci_age <- boot::boot.ci(b, index = 3, type = "basic")
        age_lower <- ci_age$basic[4]
        age_upper <- ci_age$basic[5]
      }
      tibble(
        comp_est = b$t0[1], comp_lower = comp_lower, comp_upper = comp_upper,
        prod_est = b$t0[2], prod_lower = prod_lower, prod_upper = prod_upper,
        age_est = b$t0[3], age_lower = age_lower, age_upper = age_upper,
      )
    })) |>
    select(-data) |>
    unnest(corr)
}


boot_cdi_age <- function(data) {
  data |>
    group_by(dataset_name, age_bin) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi, 2000)
      if (is.na(b$t0[1])) {
        comp_lower <- NA
        comp_upper <- NA
      } else {
        ci_comp <- boot::boot.ci(b, index = 1, type = "basic")
        comp_lower <- ci_comp$basic[4]
        comp_upper <- ci_comp$basic[5]
      }
      if (is.na(b$t0[2])) {
        prod_lower <- NA
        prod_upper <- NA
      } else {
        ci_prod <- boot::boot.ci(b, index = 2, type = "basic")
        prod_lower <- ci_prod$basic[4]
        prod_upper <- ci_prod$basic[5]
      }
      if (is.na(b$t0[3])) {
        age_lower <- NA
        age_upper <- NA
      } else {
        ci_age <- boot::boot.ci(b, index = 3, type = "basic")
        age_lower <- ci_age$basic[4]
        age_upper <- ci_age$basic[5]
      }
      tibble(
        comp_est = b$t0[1], comp_lower = comp_lower, comp_upper = comp_upper,
        prod_est = b$t0[2], prod_lower = prod_lower, prod_upper = prod_upper,
        age_est = b$t0[3], age_lower = age_lower, age_upper = age_upper,
      )
    })) |>
    select(-data) |>
    unnest(corr)
}

test_retest_corr <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    summarise(cor_test_retest = ifelse(sum(!is.na(first_admin) & !is.na(second_admin)) > 2, cor.test(first_admin, second_admin)$estimate, NA))
  return(summ$cor_test_retest[1])
}

boot_test_retest <- function(data) {
  data |>
    select(-administration_id) |>
    filter(!is.na(mean_var)) |>
    group_by(dataset_name, subject_id, pair_number) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    select(-count) |>
    pivot_wider(names_from = session_num, values_from = mean_var) |>
    group_by(dataset_name) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, test_retest_corr, 2000)
      ci <- boot::boot.ci(b, type = "basic")
      print(ci)
    if (!is.null(ci$basic) && length(ci$basic) >= 5) {
      tibble(est = b$t0, lower = ci$basic[4], upper = ci$basic[5])
    } else {
      tibble(est = b$t0, lower = NA, upper = NA)}})) |>
    select(-data) |>
    unnest(corr)
}
options(dplyr.summarise.inform = FALSE)
