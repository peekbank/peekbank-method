source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

# Trial-level exclusion flags (small lookup, not joined back to d_aoi)
trial_flags <- d_aoi |>
  group_by(dataset_name, trial_id, dataset_id, subject_id, administration_id, target_label) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE),
    .groups = "drop"
  )

d_aoi_age <- make_age_bins(d_aoi)

# Summarize to trial-level accuracy after applying exclusion filters.
# Uses trial_flags for look_both filtering via semi_join (avoids full left_join).
# Includes age_bin in grouping if present in the data.
summarize_trial_exclusion <- function(d, flags, t_start, t_end, exclude_less_than, look_both, min_trial) {
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
    group_by(administration_id, dataset_name, across(any_of("age_bin"))) |>
    mutate(
      count = n(),
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

run_trial_icc <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(
      est = map_dbl(data, \(x) get_icc(x, "accuracy")),
      num_trials = map(data, \(x) nrow(x))
    ) |>
    select(-data)
}

acc_params <- acc_trial_params_main

# Pre-compute trial-level summaries on main process
accs_summarized <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) summarize_trial_exclusion(d_aoi, trial_flags, t_s, t_e, e, l, m)
  ))

accs_summarized_age <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) summarize_trial_exclusion(d_aoi_age, trial_flags, t_s, t_e, e, l, m)
  ))

gc()

accs_icc <- accs_summarized |>
  mutate(icc = map(summary_data, run_trial_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_icc, "../cached_intermediates/4_acc_trial_icc.rds")

accs_icc_age <- accs_summarized_age |>
  mutate(icc = map(summary_data, run_trial_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_icc_age, "../cached_intermediates/4_acc_trial_icc_byage.rds")

acc_params_kid <- acc_trial_params_kid

kid_accs_summarized <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) summarize_trial_exclusion(d_aoi, trial_flags, t_s, t_e, e, l, m)
  ))

kid_accs_summarized_age <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) summarize_trial_exclusion(d_aoi_age, trial_flags, t_s, t_e, e, l, m)
  ))


kid_accs_icc <- kid_accs_summarized |>
  mutate(icc = map(summary_data, run_trial_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(kid_accs_icc, "../cached_intermediates/4_acc_kid_icc.rds")

kid_accs_icc_age <- kid_accs_summarized_age |>
  mutate(icc = map(summary_data, run_trial_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(kid_accs_icc_age, "../cached_intermediates/4_acc_kid_icc_byage.rds")
