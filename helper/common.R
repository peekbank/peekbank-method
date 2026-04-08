# source me for all analyses

library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(agreement)

# Seed for random number generation
set.seed(42)

get_age_bin_cutoff <- function(d_aoi, min_per_bin = 5) {
  d_aoi |>
    # filter(!is.na(correct)) |>
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
    filter(!(dataset_name == "adams_marchman_2018" & mean_age > 28)) |> # these do have multiple sessions but with very different items banks for the two sessions!
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
get_icc <- function(x, column = "accuracy") {
  x <- x |> filter(!is.na(.data[[column]]))
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

  return(iccs$Inter_ICC)
}
safe_cor <- function(x, y) {
  complete <- !is.na(x) & !is.na(y)
  if (sum(complete) > 2) {
    tryCatch(
      suppressWarnings(cor.test(x[complete], y[complete])$estimate[[1]]),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
}

calc_cdi <- function(data, by_age = FALSE) {
  data <- ungroup(data) |> mutate(dataset_name = as.character(dataset_name))
  if (nrow(data) == 0) {
    cols <- tibble(dataset_name = character(), comp_est = numeric(), prod_est = numeric(), age_est = numeric())
    if (by_age) cols <- mutate(cols, age_bin = character())
    return(cols)
  }
  grp <- if (by_age) c("dataset_name", "age_bin") else "dataset_name"
  data |>
    group_by(across(all_of(grp))) |>
    summarize(
      comp_est = safe_cor(mean_var, comp),
      prod_est = safe_cor(mean_var, prod),
      age_est = safe_cor(mean_var, age),
      .groups = "drop"
    )
}

calc_test_retest <- function(data) {
  na_tr_row <- tibble(dataset_name = NA_character_, est = NA_real_)
  tryCatch(
    {
      wide_data <- data |>
        select(-administration_id) |>
        filter(!is.na(mean_var)) |>
        group_by(dataset_name, subject_id, pair_number) |>
        mutate(count = n()) |>
        filter(count == 2) |>
        select(-count)

      if (nrow(wide_data) == 0) {
        return(na_tr_row)
      }

      wide_data |>
        pivot_wider(names_from = session_num, values_from = mean_var) |>
        group_by(dataset_name) |>
        summarize(est = safe_cor(first_admin, second_admin), .groups = "drop")
    },
    error = function(e) {
      warning("calc_test_retest failed: ", conditionMessage(e))
      na_tr_row
    }
  )
}

options(dplyr.summarise.inform = FALSE)
