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

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

baseline_lengths <- d_aoi |>
  group_by(dataset_name, trial_id) |>
  summarise(t_min = min(t_norm))

d_aoi_bc <- d_aoi |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

d_aoi_bc_age <- d_aoi_age |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

bc_acc_cdi <- function(b_start = -2000, b_end = 0,
                       t_start = 500, t_end = 4000) {
  d_aoi_bc |>
    group_by(dataset_name, dataset_id, administration_id, trial_id) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    group_by(administration_id, dataset_name) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T)) |>
    boot_cdi()
}



bc_acc_cdi <- function(b_start = -2000, b_end = 0,
                       t_start = 500, t_end = 4000) {
  d_aoi_bc_age |>
    group_by(dataset_name, dataset_id, administration_id, trial_id, age_bin) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
        na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
        na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    group_by(administration_id, dataset_name, age_bin) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T)) |>
    boot_cdi()
}

library(boot)
cluster <- new_cluster(36)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_library(cluster, "boot")
cluster_copy(cluster, "do_cdi")
cluster_copy(cluster, "d_aoi_bc")
cluster_copy(cluster, "d_aoi_bc_age")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "boot_cdi")
cluster_copy(cluster, "bc_acc_cdi")
cluster_copy(cluster, "bc_acc_cdi_age")


bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

bc_boot_cdi <- bc_acc_params |>
  partition(cluster) |>
  mutate(cdi = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_acc_cdi(b_s, b_e, t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(bc_boot_cdi, "../cached_intermediates/2_bc_cdi_boot.rds")

bc_boot_cdi_age <- bc_acc_params |>
  partition(cluster) |>
  mutate(cdi = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_acc_cdi_age(b_s, b_e, t_s, t_e))) |>
  collect() |>
  unnest(cdi)

saveRDS(bc_boot_cdi_age, "../cached_intermediates/2_bc_cdi_boot_byage.rds")
