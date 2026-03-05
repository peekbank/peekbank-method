source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds")

pairs_long <- make_test_retest_pairs(d_aoi)

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "test_retest_corr", "boot_test_retest")
)

d_rt_dt_long <- preprocess_rt_dt(rts) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt")

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ))

rm(d_aoi, rts, d_rt_dt_long, pairs_long)
gc()

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
