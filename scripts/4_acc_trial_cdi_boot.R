source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

trial_totals <- get_total_trials(d_aoi) |>
  distinct(dataset_name, administration_id, max_trials)

trial_flags <- d_aoi |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label
  ) |>
  summarise(
    has_correct_at_0 = any(t_norm == 0 & !is.na(correct)),
    .groups = "drop"
  )

d_aoi_age <- make_age_bins(d_aoi)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


acc_trial_cdi_summarize <- function(d, flags, t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac) {
  if (look_at_start == "yes") {
    flags <- filter(flags, has_correct_at_0)
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
    mutate(valid_n = n()) |>
    left_join(trial_totals, by = c("dataset_name", "administration_id")) |>
    mutate(frac_valid = valid_n / max_trials) |>
    filter(!is.na(max_trials), frac_valid >= min_frac) |>
    select(-valid_n, -max_trials, -frac_valid) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      count = n(), .groups = "drop"
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data, by = c("dataset_name", "administration_id"))
}

acc_params <- acc_params_trial

accs_cdi_summarized <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_trial_cdi_summarize(d_aoi, trial_flags, t_s, t_e, e, las, m, mf)
  ))

accs_cdi_summarized_age <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_trial_cdi_summarize(d_aoi_age, trial_flags, t_s, t_e, e, las, m, mf)
  ))

kid_accs_cdi_summarized <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_trial_cdi_summarize(d_aoi, trial_flags, t_s, t_e, e, las, m, mf)
  ))

kid_accs_cdi_summarized_age <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_trial_cdi_summarize(d_aoi_age, trial_flags, t_s, t_e, e, las, m, mf)
  ))

rm(d_aoi, d_aoi_age, trial_flags)
gc()

accs_cdi <- accs_cdi_summarized |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi, "../cached_intermediates/4_acc_trial_cdi.rds")

accs_cdi_byage <- accs_cdi_summarized_age |>
  mutate(cdi = map(summary_data, \(d) calc_cdi(d, by_age = TRUE))) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi_byage, "../cached_intermediates/4_acc_trial_cdi_byage.rds")


kid_accs_cdi <- kid_accs_cdi_summarized |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_cdi, "../cached_intermediates/4_acc_kid_cdi.rds")

kid_accs_cdi_byage <- kid_accs_cdi_summarized_age |>
  mutate(cdi = map(summary_data, \(d) calc_cdi(d, by_age = TRUE))) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_cdi_byage, "../cached_intermediates/4_acc_kid_cdi_byage.rds")
