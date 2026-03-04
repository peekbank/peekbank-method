source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

admins <- d_aoi |>
  select(dataset_name, subject_id, administration_id, age) |>
  distinct()

repeated <- admins |>
  group_by(dataset_name, subject_id) |>
  tally() |>
  filter(n > 1)

repeated_subjects <- admins |> inner_join(repeated)

pairs <- repeated_subjects |>
  group_by(dataset_name, subject_id) |>
  mutate(
    forward_age = lead(age),
    forward_diff = forward_age - age,
    test_num = case_when(
      forward_diff < 1.5 ~ 1,
    ),
    mean_age = case_when(
      test_num == 1 ~ (age + forward_age) / 2,
    ),
    second_admin = case_when(
      test_num == 1 ~ lead(administration_id)
    )
  ) |>
  filter(!is.na(test_num)) |>
  rename(first_admin = administration_id) |>
  select(-n, -age) |>
  left_join(repeated_subjects |> select(-age, -n), by = c("dataset_name", "subject_id", "second_admin" = "administration_id")) |>
  ungroup() |>
  mutate(pair_number = row_number()) |>
  select(-forward_age, -forward_diff, -test_num)

pairs_long <- pairs |> pivot_longer(c("first_admin", "second_admin"), names_to = "session_num", values_to = "administration_id")

pairs_aoi_data <- pairs_long |> left_join(d_aoi)

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

acc_test_retest <- function(t_start, t_end, exclude_less_than, look_both) {
  df_temp <- pairs_sim
  if (look_both == "ever") {
    df_temp <- filter(df_temp, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    df_temp <- filter(df_temp, pre_looking > 0, pre_looking < 1)
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
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(mean_var = mean(accuracy, na.rm = T), .groups = "drop") |>
    boot_test_retest()
}


library(multidplyr)
library(boot)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_library(cluster, "boot")

cluster_copy(cluster, "test_retest_corr")
cluster_copy(cluster, "boot_test_retest")
cluster_copy(cluster, "pairs_aoi_data")
cluster_copy(cluster, "acc_test_retest")


acc_params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
  look_both = c("before", "ever", "no_need")
)
accs_boot_test_retest <- acc_params |>
  partition(cluster) |>
  # head(1) |>
  mutate(corr = pmap(list(t_start, t_end, exclude_less_than, look_both), \(t_s, t_e, e, l) acc_test_retest(t_s, t_e, e, l))) |>
  collect() |>
  unnest(corr)

saveRDS(accs_boot_test_retest, "../cached_intermediates/4_acc_trial_test_retest_boot.rds")
