source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest")
)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")

rm(d_aoi, rts, d_rt_dt_long, pairs_long)
gc()


params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_bootstrap_test_retest <- rt_pairs |>
  group_by(administration_id, dataset_name, subject_id, pair_number, session_num, time_0, time_end, measure, window, during, frac) |>
  mutate(count = sum(!is.na(rt))) |>
  cross_join(params) |>
  filter(count >= start_point) |>
  select(-count) |>
  slice_sample(n = sample_down) |>
  group_by(administration_id, dataset_name, subject_id, pair_number, session_num, time_0, time_end, measure, window, during, frac, start_point, sample_down) |>
  summarize(mean_var = mean(rt, na.rm = T), .groups = "drop") |>
  group_by(measure, window, time_0, time_end, during, frac, start_point, sample_down) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_test_retest)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_bootstrap_test_retest, "../cached_intermediates/5_rt_test_retest_boot.rds")
