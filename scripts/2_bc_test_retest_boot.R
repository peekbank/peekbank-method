source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

pairs_long <- make_test_retest_pairs(d_aoi)

d_aoi_bc <- make_baseline_corrected(d_aoi)
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
      bc_accuracy = window_accuracy - baseline_accuracy
    ) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    summarize(mean_var = mean(bc_accuracy, na.rm = T), .groups = "drop") |>
    calc_test_retest()
}




bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

bc_test_retest_result <- bc_acc_params |>
  mutate(corr = pmap(list(b_start, b_end, t_start, t_end), \(b_s, b_e, t_s, t_e) bc_test_retest(b_s, b_e, t_s, t_e))) |>
  unnest(corr)

saveRDS(bc_test_retest_result, "../cached_intermediates/2_bc_test_retest.rds")
