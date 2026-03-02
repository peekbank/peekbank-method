source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- d_aoi |>
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
  filter(count >= 5) |>
  ungroup()

d_aoi_age <- d_aoi |> inner_join(age_bin_cutoff)

baseline_lengths <- d_aoi |>
  group_by(dataset_name, trial_id) |>
  summarise(t_min = min(t_norm))

d_aoi_bc <- d_aoi |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

d_aoi_bc_age <- d_aoi_age |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

icc_bc_window_bootstrap_sim <- function(b_start = -2000, b_end = 0,
                                        t_start = 500, t_end = 4000,
                                        object) {
  # get baseline corrected accuracies for all trials with ANY baseline info
  df <- d_aoi_bc |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    filter(!is.na(bc_accuracy)) |>
    group_by(dataset_id, dataset_name, administration_id, target_label) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name) |>
    nest() |>
    mutate(icc = map(data, \(d) bootstrap_icc(d, "bc_accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}


icc_bc_window_age_bootstrap_sim <- function(b_start = -2000, b_end = 0,
                                            t_start = 500, t_end = 4000,
                                            object) {
  # get baseline corrected accuracies for all trials with ANY baseline info
  df <- d_aoi_bc_age |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id, age_bin) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    filter(!is.na(bc_accuracy)) |>
    group_by(dataset_id, dataset_name, administration_id, target_label, age_bin) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name, age_bin) |>
    nest() |>
    mutate(icc = map(data, \(d) bootstrap_icc(d, "bc_accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}

bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
  object = c("administration")
)

library(multidplyr)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "d_aoi_bc")
cluster_copy(cluster, "d_aoi_bc_age")
cluster_copy(cluster, "icc_bc_window_bootstrap_sim")
cluster_copy(cluster, "icc_bc_window_age_bootstrap_sim")


bc_accs <- bc_acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(b_start, b_end, t_start, t_end, object), icc_bc_window_bootstrap_sim)) |>
  collect() |>
  unnest(col = icc)

saveRDS(bc_accs, "../cached_intermediates/2_bc_icc_boot.rds")

bc_accs_age <- bc_acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(b_start, b_end, t_start, t_end, object), icc_bc_window_age_bootstrap_sim)) |>
  collect() |>
  unnest(col = icc)

saveRDS(bc_accs_age, "../cached_intermediates/2_bc_icc_boot_byage.rds")
