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

baseline_lengths <- d_aoi |>
  group_by(dataset_name, trial_id) |>
  summarise(t_min = min(t_norm))

d_aoi_bc <- d_aoi |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

pairs_aoi_data_bc <- pairs_long |> left_join(d_aoi_bc)



bc_test_retest <- function(b_start, b_end, t_start = -500, t_end = 4000) {

  pairs_aoi_data_bc |>
    group_by(dataset_name, dataset_id, administration_id, subject_id, pair_number, session_num, trial_id) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end],
                             na.rm = TRUE
      ),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end],
                               na.rm = TRUE
      ),
      bc_accuracy = window_accuracy - baseline_accuracy) |> 
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T), .groups = "drop") |>
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
cluster_copy(cluster, "pairs_aoi_data_bc")
cluster_copy(cluster, "bc_test_retest")


bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

bc_boot_test_retest <- bc_acc_params |>
  partition(cluster) |>
  mutate(corr = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_test_retest(b_s, b_e, t_s, t_e))) |>
  collect() |>
  unnest(corr)

saveRDS(bc_boot_test_retest, "../cached_intermediates/2_bc_test_retest_boot.rds")
