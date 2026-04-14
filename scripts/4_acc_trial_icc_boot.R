source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

trial_totals <- get_total_trials(d_aoi) |>
  distinct(dataset_name, administration_id, max_trials)

# Trial-level flags for semi_join (t_norm == 0 is onset)
trial_flags <- d_aoi |>
  group_by(dataset_name, trial_id, dataset_id, subject_id, administration_id, target_label) |>
  summarise(
    has_correct_at_0 = any(t_norm == 0 & !is.na(correct)),
    .groups = "drop"
  )

d_aoi_age <- make_age_bins(d_aoi)

# Summarize to trial-level accuracy after applying exclusion filters.
# min_frac: numerator = trials passing window + prop + look_at_start; denominator = max_trials from get_total_trials.
summarize_trial_exclusion <- function(d, flags, t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac) {
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
    mutate(
      # frac_valid: trials retained under this param row / max trials that administration could contribute
      frac_valid = valid_n / max_trials
    ) |>
    filter(!is.na(max_trials), frac_valid >= min_frac) |>
    select(-valid_n, -max_trials, -frac_valid) |>
    mutate(count = n()) |>
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

acc_params <- acc_params_trial

accs_summarized <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) summarize_trial_exclusion(d_aoi, trial_flags, t_s, t_e, e, las, m, mf)
  ))

accs_summarized_age <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) summarize_trial_exclusion(d_aoi_age, trial_flags, t_s, t_e, e, las, m, mf)
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

kid_accs_summarized <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) summarize_trial_exclusion(d_aoi, trial_flags, t_s, t_e, e, las, m, mf)
  ))

kid_accs_summarized_age <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) summarize_trial_exclusion(d_aoi_age, trial_flags, t_s, t_e, e, las, m, mf)
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
