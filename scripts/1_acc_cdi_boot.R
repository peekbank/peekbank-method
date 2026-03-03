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


acc_cdi <- function(t_start = -500, t_end = 4000) {
  d_aoi |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    group_by(dataset_name, dataset_id, administration_id, target_label) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name) |>
    summarize(mean_var = mean(accuracy, na.rm = T))
}



acc_cdi_age <- function(t_start = -500, t_end = 4000) {
  d_aoi_age |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, age_bin, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, age_bin) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, age_bin) |>
    summarize(mean_var = mean(accuracy, na.rm = T))
}

library(boot)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_library(cluster, "boot")
cluster_copy(cluster, "do_cdi")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "boot_cdi")
cluster_copy(cluster, "boot_cdi_age")


acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)

accs_boot_cdi <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi(t_s, t_e))) |>
  partition() |>
  mutate(cdi = map(summary_data, boot_cdi)) |>
  collect() |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "../cached_intermediates/1_acc_cdi_boot.rds")

accs_boot_cdi_byage <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) acc_cdi_age(t_s, t_e))) |>
  partition() |>
  mutate(cdi = map(summary_data, boot_cdi_age)) |>
  collect() |>
  unnest(cdi)

saveRDS(accs_boot_cdi_byage, "../cached_intermediates/1_acc_cdi_boot_byage.rds")
