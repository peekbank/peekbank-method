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

icc_window_sim_age_bootstrap <- function(t_start = -500, t_end = 4000, object) {
  print(paste(t_start, t_end))

  df <- d_aoi_age |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id, age_bin) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    filter(!is.na(accuracy)) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, age_bin) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name, age_bin) |>
    nest() |>
    mutate(icc = map(data, \(d) bootstrap_icc(d, "accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}

icc_window_sim_bootstrap <- function(t_start = -500, t_end = 4000, object) {
  print(paste(t_start, t_end))

  df <- d_aoi |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    filter(!is.na(accuracy)) |>
    group_by(dataset_name, dataset_id, administration_id, target_label) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name) |>
    nest() |>
    mutate(icc = map(data, \(d) bootstrap_icc(d, "accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}

acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
  object = c("administration")
)

cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "icc_window_sim_age_bootstrap")
cluster_copy(cluster, "icc_window_sim_bootstrap")
cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "d_aoi")
cluster_copy(cluster, "d_aoi_age")


accs_boot <- acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(t_start, t_end, object), \(t_s, t_e, o) icc_window_sim_bootstrap(t_s, t_e, o))) |>
  collect() |>
  unnest(col = icc)


saveRDS(accs_boot, "../cached_intermediates/1_acc_icc_boot.rds")

accs_boot_age <- acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(t_start, t_end, object), \(t_s, t_e, o) icc_window_sim_age_bootstrap(t_s, t_e, o))) |>
  collect() |>
  unnest(col = icc)

saveRDS(accs_boot_age, "../cached_intermediates/1_acc_icc_boot_byage.rds")
