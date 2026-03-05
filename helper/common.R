# source me for all analyses

library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(agreement)
library(multidplyr)
library(boot)
# because we don't have full on tidyverse on the cluster, we specify the parts we actually need

# Seed for random number generation
set.seed(42)

get_age_bin_cutoff <- function(d_aoi, min_per_bin = 5) {
  d_aoi |>
    filter(!is.na(correct)) |>
    distinct(administration_id, age, dataset_name) |>
    mutate(age_bin = case_when(
      age < 18 ~ "<18",
      age < 24 ~ "18-24",
      age < 36 ~ "24-36",
      age >= 36 ~ ">=36"
    )) |>
    group_by(dataset_name, age_bin) |>
    mutate(count = n()) |>
    filter(count >= min_per_bin) |>
    ungroup()
}

make_age_bins <- function(d_aoi, min_per_bin = 5) {
  d_aoi |> inner_join(get_age_bin_cutoff(d_aoi, min_per_bin))
}

make_test_retest_pairs <- function(d_aoi) {
  admins <- d_aoi |>
    select(dataset_name, subject_id, administration_id, age) |>
    distinct()

  repeated <- admins |>
    group_by(dataset_name, subject_id) |>
    tally() |>
    filter(n > 1)

  repeated_subjects <- admins |> inner_join(repeated)

  pairs <- repeated_subjects |>
    group_by(dataset_name, subject_id) |>
    mutate(
      forward_age = lead(age),
      forward_diff = forward_age - age,
      test_num = case_when(
        forward_diff < 1.5 ~ 1,
      ),
      mean_age = case_when(
        test_num == 1 ~ (age + forward_age) / 2,
      ),
      second_admin = case_when(
        test_num == 1 ~ lead(administration_id)
      )
    ) |>
    filter(!is.na(test_num)) |>
    rename(first_admin = administration_id) |>
    select(-n, -age) |>
    left_join(repeated_subjects |> select(-age, -n), by = c("dataset_name", "subject_id", "second_admin" = "administration_id")) |>
    ungroup() |>
    mutate(pair_number = row_number()) |>
    select(-forward_age, -forward_diff, -test_num)

  pairs |> pivot_longer(c("first_admin", "second_admin"), names_to = "session_num", values_to = "administration_id")
}

make_baseline_corrected <- function(d_aoi) {
  baseline_lengths <- d_aoi |>
    group_by(dataset_name, trial_id) |>
    summarise(t_min = min(t_norm), .groups = "drop")

  d_aoi |>
    left_join(baseline_lengths, by = c("dataset_name", "trial_id")) |>
    filter(t_min < 0)
}

# key ICC function
get_icc <- function(x, column = "accuracy", object = "stimulus") {
  x <- x |> filter(!is.na(.data[[column]]))
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
  x <- x |> filter(!is.na(.data[[column]]))
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

safe_boot_ci <- function(b, index) {
  t0 <- b$t0[index]
  if (is.na(t0)) {
    return(list(lower = NA_real_, upper = NA_real_))
  }
  t_vals <- b$t[, index]
  t_no_na <- t_vals[!is.na(t_vals)]
  if (length(t_no_na) == 0) {
    return(list(lower = NA_real_, upper = NA_real_))
  }
  b$R <- length(t_no_na)
  ci <- boot::boot.ci(boot.out = b, t = t_no_na, t0 = t0, type = "basic")
  if (!is.null(ci$basic) && length(ci$basic) >= 5) {
    list(lower = ci$basic[4], upper = ci$basic[5])
  } else {
    list(lower = NA_real_, upper = NA_real_)
  }
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

boot_cdi <- function(data, by_age = FALSE) {
  grp <- if (by_age) c("dataset_name", "age_bin") else "dataset_name"
  data |>
    group_by(across(all_of(grp))) |>
    nest() |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, do_cdi, 2000)
      ci_comp <- safe_boot_ci(b, 1)
      ci_prod <- safe_boot_ci(b, 2)
      ci_age <- safe_boot_ci(b, 3)
      tibble(
        comp_est = b$t0[1], comp_lower = ci_comp$lower, comp_upper = ci_comp$upper,
        prod_est = b$t0[2], prod_lower = ci_prod$lower, prod_upper = ci_prod$upper,
        age_est = b$t0[3], age_lower = ci_age$lower, age_upper = ci_age$upper,
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
      ci <- safe_boot_ci(b, 1)
      tibble(est = b$t0, lower = ci$lower, upper = ci$upper)
    })) |>
    select(-data) |>
    unnest(corr)
}

setup_cluster <- function(libs, copy_names = character(), envir = parent.frame(), n_workers = 16) {
  cluster <- new_cluster(n_workers)
  for (lib in libs) cluster_library(cluster, lib)
  for (name in copy_names) cluster_copy(cluster, name, env = envir)
  cluster
}

options(dplyr.summarise.inform = FALSE)
