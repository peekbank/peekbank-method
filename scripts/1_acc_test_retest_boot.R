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


test_retest_corr <- function(data, indices) {
  summ <- data |>
    slice(indices) |>
    summarise(cor_test_retest = ifelse(sum(!is.na(first_admin) & !is.na(second_admin)) > 2, cor.test(first_admin, second_admin)$estimate, NA))
  return(summ$cor_test_retest[1])
}

boot_test_retest <- function(data) {
  data |>
    select(-administration_id) |>
    pivot_wider(names_from = session_num, values_from = mean_var) |>
    group_by(dataset_name) |>
    nest() |>
    # head(1) |>
    mutate(corr = map(data, \(d) {
      b <- boot::boot(d, test_retest_corr, 2000)
      ci <- boot::boot.ci(b, type = "basic")
      print(ci)
      tibble(est = b$t0, lower = ci$basic[4], upper = ci$basic[5])
    })) |>
    select(-data) |>
    unnest(corr)
}

acc_test_retest <- function(t_start = -500, t_end = 4000) {
  print(paste(t_start, t_end))

  pairs_aoi_data |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
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
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)
accs_boot_test_retest <- acc_params |>
  partition(cluster) |>
  # head(1) |>
  mutate(corr = pmap(list(t_start, t_end), \(t_s, t_e) acc_test_retest(t_s, t_e))) |>
  collect() |>
  unnest(corr)

saveRDS(accs_boot_test_retest, "../cached_intermediates/1_acc_test_retest_boot.rds")
