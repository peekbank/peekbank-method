source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_sim <- d_aoi |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE)
  ) |>
  left_join(d_aoi)

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

d_sim_age <- d_aoi_age |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label, age_bin
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE)
  ) |>
  left_join(d_aoi_age)

icc_trial_exclusion_age_bootstrap <- function(t_start, t_end, exclude_less_than, look_both) {
  df_temp <- d_sim_age
  if (look_both == "ever") {
    df_temp <- filter(df_temp, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    df_temp <- filter(df_temp, pre_looking > 0, pre_looking < 1)
  }

  df <- df_temp |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id, age_bin) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, age_bin) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name, age_bin) |>
    nest() |>
    mutate(
      icc = map(data, \(d) bootstrap_icc(d, "accuracy", 2000)),
      num_trials = map(data, \(d) nrow(d))
    ) |>
    select(-data) |>
    unnest(icc)
}

icc_trial_exclusion_bootstrap <- function(t_start, t_end, exclude_less_than, look_both) {
  df_temp <- d_sim
  if (look_both == "ever") {
    df_temp <- filter(df_temp, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    df_temp <- filter(df_temp, pre_looking > 0, pre_looking < 1)
  }

  df <- df_temp |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(dataset_name, dataset_id, administration_id, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct))
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(dataset_name, dataset_id, administration_id, target_label) |>
    mutate(repetition = row_number())

  # compute ICCs
  df |>
    group_by(dataset_name) |>
    nest() |>
    mutate(
      icc = map(data, \(d) bootstrap_icc(d, "accuracy", 2000)),
      num_trials = map(data, \(d) nrow(d))
    ) |>
    select(-data) |>
    unnest(icc)
}

acc_params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
  look_both = c("before", "ever", "no_need")
)

cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "icc_trial_exclusion_age_bootstrap")
cluster_copy(cluster, "icc_trial_exclusion_bootstrap")
cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "d_sim")
cluster_copy(cluster, "d_sim_age")


accs_boot <- acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(t_start, t_end, exclude_less_than, look_both), \(t_s, t_e, e, l) icc_trial_exclusion_bootstrap(t_s, t_e, e, l))) |>
  collect() |>
  unnest(col = icc)


saveRDS(accs_boot, "../cached_intermediates/4_acc_trial_icc_boot.rds")

accs_boot_age <- acc_params |>
  partition(cluster) |>
  mutate(icc = pmap(list(t_start, t_end, exclude_less_than, look_both), \(t_s, t_e, e, l) icc_trial_exclusion_bootstrap(t_s, t_e, e, l))) |>
  collect() |>
  unnest(col = icc)

saveRDS(accs_boot_age, "../cached_intermediates/4_acc_trial_icc_boot_byage.rds")
