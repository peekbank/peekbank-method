source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

pairs_long <- make_test_retest_pairs(d_aoi)
pairs_aoi_data <- pairs_long |> left_join(d_aoi)

rm(d_aoi, pairs_long)
gc()

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


cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest", "pairs_aoi_data", "acc_test_retest")
)


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
