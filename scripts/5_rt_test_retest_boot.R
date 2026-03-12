source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest", "empirical_ci")
)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")


params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_bootstrap_test_retest <- params |>
  pmap_dfr(\(start_point, sample_down) {
    rt_pairs |>
      group_by(administration_id, dataset_name, subject_id, pair_number, session_num, time_0, time_end, measure, window, during, frac) |>
      mutate(count = sum(!is.na(rt))) |>
      filter(count >= start_point) |>
      select(-count) |>
      cross_join(tibble(iteration = 1:1000)) |>
      group_by(administration_id, dataset_name, subject_id, pair_number, session_num, time_0, time_end, measure, window, during, frac, iteration) |>
      slice_sample(n = sample_down) |>
      ungroup() |>
      mutate(start_point = start_point, sample_down = sample_down)
  }) |>
  group_by(administration_id, dataset_name, subject_id, pair_number, session_num, time_0, time_end, measure, window, during, frac, start_point, sample_down, iteration) |>
  summarize(mean_var = mean(rt, na.rm = T), .groups = "drop") |>
  group_by(measure, window, time_0, time_end, during, frac, start_point, sample_down) |>
  nest() |>
  partition(cluster) |>
  mutate(corr = map(data, \(d) {
    wide_data <- d |>
      select(-administration_id) |>
      filter(!is.na(mean_var)) |>
      group_by(dataset_name, subject_id, pair_number, iteration) |>
      mutate(count = n()) |>
      filter(count == 2) |>
      select(-count)

    wide_data |>
      pivot_wider(names_from = session_num, values_from = mean_var) |>
      group_by(dataset_name, iteration) |>
      nest() |>
      mutate(corr = map(data, \(d) test_retest_corr(d, 1:nrow(d)))) |>
      select(-data) |>
      unnest(corr) |>
      ungroup() |>
      empirical_ci()
  })) |>
  collect() |>
  select(-data) |>
  unnest(corr)

saveRDS(rt_bootstrap_test_retest, "../cached_intermediates/5_rt_test_retest_boot.rds")
