source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

# Trial-level exclusion flags (small lookup, not joined back to d_aoi)
trial_flags <- d_aoi |>
  group_by(dataset_name, trial_id, dataset_id, subject_id, administration_id, target_label) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE),
    .groups = "drop"
  )

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

# Summarize to trial-level accuracy after applying exclusion filters.
# Uses trial_flags for look_both filtering via semi_join (avoids full left_join).
# Includes age_bin in grouping if present in the data.
summarize_trial_exclusion <- function(d, flags, t_start, t_end, exclude_less_than, look_both) {
  if (look_both == "ever") {
    flags <- filter(flags, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    flags <- filter(flags, pre_looking > 0, pre_looking < 1)
  }

  join_cols <- c("dataset_name", "trial_id", "dataset_id", "administration_id", "target_label")

  d |>
    semi_join(flags, by = join_cols) |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
    ))), across(any_of("age_bin"))) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

# Bootstrap ICCs on pre-computed trial-level summaries.
run_trial_icc_bootstrap <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(
      icc = map(data, \(x) bootstrap_icc(x, "accuracy", 2000)),
      num_trials = map(data, \(x) nrow(x))
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

# Pre-compute trial-level summaries on main process
accs_summarized <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end, exclude_less_than, look_both),
    \(t_s, t_e, e, l) summarize_trial_exclusion(d_aoi, trial_flags, t_s, t_e, e, l)))

accs_summarized_age <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end, exclude_less_than, look_both),
    \(t_s, t_e, e, l) summarize_trial_exclusion(d_aoi_age, trial_flags, t_s, t_e, e, l)))

rm(d_aoi, d_aoi_age, age_bin_cutoff, trial_flags)
gc()

# Only bootstrap_icc and run_trial_icc_bootstrap go to workers (NOT d_aoi or d_sim)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "run_trial_icc_bootstrap")

accs_boot <- accs_summarized |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_trial_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_boot, "../cached_intermediates/4_acc_trial_icc_boot.rds")

accs_boot_age <- accs_summarized_age |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_trial_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_boot_age, "../cached_intermediates/4_acc_trial_icc_boot_byage.rds")
