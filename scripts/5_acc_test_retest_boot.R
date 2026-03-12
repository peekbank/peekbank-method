source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

pairs_long <- make_test_retest_pairs(d_aoi)
pairs_aoi_data <- pairs_long |> left_join(d_aoi)

rm(d_aoi, pairs_long)
gc()

acc_downsample_test_retest <- function(t_start, t_end, start_point, sample_down, iter) {
  d <- pairs_aoi_data |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, target_label, trial_id) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    mutate(count = n()) |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iter)) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, iteration) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, iteration) |>
    summarize(mean_var = mean(accuracy, na.rm = T), .groups = "drop")

  wide_data <- d |>
    select(-administration_id) |>
    filter(!is.na(mean_var)) |>
    group_by(dataset_name, subject_id, pair_number) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    select(-count)

  wide_data |>
    pivot_wider(names_from = session_num, values_from = mean_var) |>
    group_by(dataset_name) |>
    nest() |>
    test_retest_corr( 1:nrow(d))
}


cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest", "pairs_aoi_data", "acc_downsample_test_retest")
)


params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)

acc_downsample <- params |>
  mutate(iteration = 100) |>
  partition(cluster) |>
  # head(1) |>
  mutate(corr = pmap(list(t_start, t_end, start_point, sample_down, iteration), \(t_s, t_e, s_p, s_d, iteration) acc_downsample_test_retest(t_s, t_e, s_p, s_d, iteration))) |>
  collect() |>
  unnest(corr)

saveRDS(acc_downsample, "../cached_intermediates/5_acc_downsample_test_retest_boot.rds")
