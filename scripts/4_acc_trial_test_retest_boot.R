source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

trial_totals <- get_total_trials(d_aoi) |>
  distinct(dataset_name, administration_id, max_trials)

pairs_long <- make_test_retest_pairs(d_aoi)
pairs_aoi_data <- pairs_long |> left_join(d_aoi)

pairs_sim <- pairs_aoi_data |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label, pair_number, session_num
  ) |>
  summarise(
    has_correct_at_0 = any(t_norm == 0 & !is.na(correct)),
    .groups = "drop"
  ) |>
  left_join(pairs_aoi_data)


acc_test_retest <- function(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac) {
  df_temp <- pairs_sim
  if (look_at_start == "yes") {
    df_temp <- filter(df_temp, has_correct_at_0)
  }

  df_temp |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name) |>
    mutate(valid_n = n_distinct(trial_id)) |>
    left_join(trial_totals, by = c("dataset_name", "administration_id")) |>
    mutate(frac_valid = valid_n / max_trials) |>
    filter(!is.na(max_trials), frac_valid >= min_frac) |>
    select(-valid_n, -max_trials, -frac_valid) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      count = n(),
      .groups = "drop",
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    group_by(dataset_name, subject_id, pair_number) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    ungroup() |>
    calc_test_retest()
}


acc_params <- acc_params_trial

accs_test_retest <- acc_params |>
  mutate(corr = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_test_retest(t_s, t_e, e, las, m, mf)
  )) |>
  unnest(corr)

saveRDS(accs_test_retest, "../cached_intermediates/4_acc_trial_test_retest.rds")

kid_accs_test_retest <- acc_params_kid |>
  mutate(corr = pmap(
    list(t_start, t_end, exclude_less_than, look_at_start, min_trial, min_frac),
    \(t_s, t_e, e, las, m, mf) acc_test_retest(t_s, t_e, e, las, m, mf)
  )) |>
  unnest(corr)

saveRDS(kid_accs_test_retest, "../cached_intermediates/4_acc_kid_test_retest.rds")
