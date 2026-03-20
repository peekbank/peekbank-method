source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

pairs_long <- make_test_retest_pairs(d_aoi)
# pairs_aoi_data <- pairs_long |> left_join(d_aoi)

pairs_aoi_age <- pairs_long |>
  mutate(age_bin = case_when(
    mean_age < 18 ~ "<18",
    mean_age < 24 ~ "18-24",
    mean_age < 36 ~ "24-36",
    mean_age >= 36 ~ ">=36"
  )) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= min_per_bin) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= 5) |>
  ungroup() |>
  filter(age_bin %in% c("18-24", "24-36")) |>
  left_join(d_aoi)

# pairs_sim <- pairs_aoi_data |>
#   group_by(
#     dataset_name, trial_id, dataset_id, subject_id, administration_id,
#     target_label, pair_number, session_num
#   ) |>
#   summarise(
#     total_target_prop = mean(correct, na.rm = TRUE),
#     pre_looking = mean(correct[t_norm < 400], na.rm = TRUE)
#   ) |>
#   left_join(pairs_aoi_data)

pairs_sim_age <- pairs_aoi_data_age |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label, pair_number, session_num, age_bin
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE)
  ) |>
  left_join(pairs_aoi_data_age)

# rm(d_aoi, pairs_aoi_data, pairs_long)
# gc()

# acc_test_retest <- function(t_start, t_end, exclude_less_than, look_both, min_trial) {
#   df_temp <- pairs_sim
#   if (look_both == "ever") {
#     df_temp <- filter(df_temp, total_target_prop > 0, total_target_prop < 1)
#   } else if (look_both == "before") {
#     df_temp <- filter(df_temp, pre_looking > 0, pre_looking < 1)
#   }
#
#   df_temp |>
#     filter(t_norm > t_start, t_norm < t_end) |>
#     group_by(administration_id, dataset_name, subject_id, pair_number, session_num, target_label, trial_id, ) |>
#     summarise(
#       accuracy = mean(correct, na.rm = TRUE),
#       prop_data = mean(!is.na(correct)),
#       .groups = "drop"
#     ) |>
#     filter(prop_data >= exclude_less_than) |>
#     filter(!is.na(accuracy)) |>
#     group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
#     summarize(
#       mean_var = mean(accuracy, na.rm = T),
#       count = n(),
#       .groups = "drop",
#     ) |>
#     filter(count >= min_trial) |>
#     select(-count) |>
#     group_by(dataset_name, subject_id, pair_number) |>
#     mutate(count = n()) |>
#     filter(count == 2) |>
#     ungroup() |>
#     boot_test_retest()
# }


acc_test_retest_age <- function(t_start, t_end, exclude_less_than, look_both, min_trial) {
  df_temp <- pairs_sim_age
  if (look_both == "ever") {
    df_temp <- filter(df_temp, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    df_temp <- filter(df_temp, pre_looking > 0, pre_looking < 1)
  }

  df_temp |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, target_label, trial_id, age_bin) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, age_bin) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      count = n(),
      .groups = "drop",
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    group_by(dataset_name, subject_id, pair_number, age_bin) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    ungroup() |>
    mutate(dataset_name = str_c(age_bin, dataset_name)) |>
    boot_test_retest()
}
cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c(
    "safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest",
    "pairs_sim", "acc_test_retest",
    "pairs_sim_age", "acc_test_retest_age"
  )
)


acc_params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
  look_both = c("before", "ever", "no_need"),
  min_trial = 1
)

acc_params_kid <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .2, .5),
  look_both = c("no_need"),
  min_trial = c(1, 2, 3, 4, 5, 8, 10, 15)
)
#
#  accs_boot_test_retest <- acc_params |>
#    partition(cluster) |>
#    # head(1) |>
#    mutate(corr = pmap(list(t_start, t_end, exclude_less_than, look_both, min_trial), \(t_s, t_e, e, l, m) acc_test_retest(t_s, t_e, e, l, m))) |>
#    collect() |>
#    unnest(corr)
#
#  saveRDS(accs_boot_test_retest, "../cached_intermediates/4_acc_trial_test_retest_boot.rds")
#
# kid_accs_boot_test_retest <- acc_params_kid |>
#   partition(cluster) |>
#   mutate(corr = pmap(list(t_start, t_end, exclude_less_than, look_both, min_trial), \(t_s, t_e, e, l, m) acc_test_retest(t_s, t_e, e, l, m))) |>
#   collect() |>
#   unnest(corr)
#
# saveRDS(kid_accs_boot_test_retest, "../cached_intermediates/4_acc_kid_test_retest_boot.rds")


accs_boot_test_retest_age <- acc_params |>
  partition(cluster) |>
  # head(1) |>
  mutate(corr = pmap(list(t_start, t_end, exclude_less_than, look_both, min_trial), \(t_s, t_e, e, l, m) acc_test_retest_age(t_s, t_e, e, l, m))) |>
  collect() |>
  unnest(corr)

saveRDS(accs_boot_test_retest_age, "../cached_intermediates/4_acc_trial_test_retest_boot_age.rds")

kid_accs_boot_test_retest_age <- acc_params_kid |>
  partition(cluster) |>
  mutate(corr = pmap(list(t_start, t_end, exclude_less_than, look_both, min_trial), \(t_s, t_e, e, l, m) acc_test_retest_age(t_s, t_e, e, l, m))) |>
  collect() |>
  unnest(corr)

saveRDS(kid_accs_boot_test_retest_age, "../cached_intermediates/4_acc_kid_test_retest_boot_age.rds")
