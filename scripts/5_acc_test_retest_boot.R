source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

pairs_long <- make_test_retest_pairs(d_aoi)
pairs_aoi_data <- pairs_long |> left_join(d_aoi)

rm(d_aoi, pairs_long)
gc()

acc_downsample_test_retest <- function(t_start, t_end, start_point, sample_down) {
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
    mutate(count = n()) |>
    filter(count >= start_point) |>
    select(-count) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(mean_var = mean(accuracy, na.rm = T), .groups = "drop") |>
    boot_test_retest()
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
  partition(cluster) |>
  # head(1) |>
  mutate(corr = pmap(list(t_start, t_end, start_point, sample_down), \(t_s, t_e, s_p, s_d) acc_downsample_test_retest(t_s, t_e, s_p, s_d))) |>
  collect() |>
  unnest(corr)

saveRDS(acc_downsample, "../cached_intermediates/5_acc_downsample_test_retest_boot.rds")
