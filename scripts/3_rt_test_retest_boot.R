source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds")

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

d_rt_dt_long <- rts |>
  filter(shift_type == "D-T") |>
  mutate(
    land_rt = rt,
    first_launch_rt = shift_start_rt,
  ) |>
  mutate(
    across(c("land_rt", "first_launch_rt"), log, .names = "log_{.col}"),
    across(
      c(
        "land_rt", "first_launch_rt",
        "log_land_rt", "log_first_launch_rt",
      ),
      ~ ifelse(shift_length <= 600, .x, NA),
      .names = "trim_first_{.col}"
    ),
    across(
      c(
        "land_rt",
        "log_land_rt",
      ),
      ~ ifelse(last_shift_length <= 600, .x, NA),
      .names = "trim_last_{.col}"
    )
  ) |>
  select(-rt, -shift_start_rt, -last_shift_rt) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt")

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ))


rt_bootstrap_test_retest <- rt_pairs |>
  group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  group_by(measure, window, time_0, time_end, during, frac) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_test_retest)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_bootstrap_test_retest, "../cached_intermediates/3_rt_test_retest_boot.rds")
