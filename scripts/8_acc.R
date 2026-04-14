source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

pairs_aoi_data <- make_test_retest_pairs(d_aoi) |>
  left_join(d_aoi)

pairs_sim <- pairs_aoi_data |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label, pair_number, session_num
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE)
  ) |>
  left_join(pairs_aoi_data)
trial_flags <- d_aoi |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE),
    .groups = "drop"
  )

prep_accuracy <- function(d, flags, t_start, t_end, exclude_less_than, look_both, min_trial) {
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
    ))), across(any_of(c("age_bin", "subject_id", "pair_number", "session_num")))) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, across(any_of(c("age_bin", "subject_id", "pair_number", "session_num")))) |>
    mutate(
      count = n(),
    ) |>
    filter(count >= min_trial) |>
    select(-count)
}

prep_bc <- function(d, flags, t_start, t_end, b_start, b_end, exclude_less_than, look_both, min_trial) {
  if (look_both == "ever") {
    flags <- filter(flags, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    flags <- filter(flags, pre_looking > 0, pre_looking < 1)
  }

  join_cols <- c("dataset_name", "trial_id", "dataset_id", "administration_id", "target_label")

  d |>
    semi_join(flags, by = join_cols) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
    ))), across(any_of(c("age_bin", "subject_id", "pair_number", "session_num")))) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end], na.rm = TRUE),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end], na.rm = TRUE),
      bc_accuracy = window_accuracy - baseline_accuracy,
      prop_data = mean(!is.na(correct[t_norm >= t_start & t_norm <= t_end])),
      .groups = "drop"
    ) |>
    filter(!is.na(bc_accuracy)) |>
    filter(prop_data >= exclude_less_than) |>
    group_by(administration_id, dataset_name, across(any_of(c("age_bin", "subject_id", "pair_number", "session_num")))) |>
    mutate(
      count = n(),
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    rename(accuracy = bc_accuracy)
}

for_icc <- function(d) {
  d |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

run_icc <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(est = map_dbl(data, \(x) get_icc(x, "accuracy"))) |>
    select(-data)
}


cdi_summarize <- function(d) {
  d |>
    group_by(administration_id, dataset_name, across(any_of("age_bin"))) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      count = n(), .groups = "drop"
    ) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data)
}

do_test_retest <- function(d) {
  d |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      .groups = "drop",
    ) |>
    group_by(dataset_name, subject_id, pair_number) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    ungroup() |>
    calc_test_retest()
}

acc_params <- acc_params_summary_recommended
bc_params <- bc_params_summary_alternative


accs_icc <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) prep_accuracy(d_aoi, trial_flags, t_s, t_e, e, l, m) |> for_icc()
  )) |>
  mutate(icc = map(summary_data, run_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_icc, "../cached_intermediates/8_acc_icc.rds")


accs_cdi_summarized <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) prep_accuracy(d_aoi, trial_flags, t_s, t_e, e, l, m) |> cdi_summarize()
  )) |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi_summarized, "../cached_intermediates/8_acc_cdi.rds")

accs_boot_test_retest <- acc_params |>
  mutate(corr = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) prep_accuracy(pairs_sim, trial_flags, t_s, t_e, e, l, m) |> do_test_retest()
  )) |>
  unnest(corr)

saveRDS(accs_boot_test_retest, "../cached_intermediates/8_acc_test_retest.rds")


bc_icc_summarized <- bc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, b_start, b_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, b_s, b_e, e, l, m) prep_bc(d_aoi, trial_flags, t_s, t_e, b_s, b_e, e, l, m) |> for_icc()
  )) |>
  mutate(icc = map(summary_data, run_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(bc_icc_summarized, "../cached_intermediates/8_bc_icc.rds")


bc_cdi_summarized <- bc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, b_start, b_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, b_s, b_e, e, l, m) prep_bc(d_aoi, trial_flags, t_s, t_e, b_s, b_e, e, l, m) |> cdi_summarize()
  )) |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(bc_cdi_summarized, "../cached_intermediates/8_bc_cdi.rds")

bc_boot_test_retest <- bc_params |>
  mutate(corr = pmap(
    list(t_start, t_end, b_start, b_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, b_s, b_e, e, l, m) prep_bc(pairs_sim, trial_flags, t_s, t_e, b_s, b_e, e, l, m) |> do_test_retest()
  )) |>
  unnest(corr)

saveRDS(bc_boot_test_retest, "../cached_intermediates/8_bc_test_retest.rds")
