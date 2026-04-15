source("../helper/common.R")
source("../helper/rt_helper.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> left_join(rt_params)

rts_weird <- readRDS("../cached_intermediates/3_weird_rts.rds") |> left_join(rt_params_weird)

pairs_long <- make_test_retest_pairs(d_aoi)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, min_rt, max_rt, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")


rt_pairs_weird <- pairs_long |>
  left_join(preprocess_rt_dt(rts_weird) |> filter(measure %in% c("log_land_rt", "land_rt")) |> select(
    administration_id, rt, measure, min_rt, max_rt, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")

rt_test_retest <- rt_pairs |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_test_retest, "../cached_intermediates/3_rt_test_retest.rds")

rt_test_retest_weird <- rt_pairs_weird |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_test_retest_weird, "../cached_intermediates/3_rt_test_retest_weird.rds")
