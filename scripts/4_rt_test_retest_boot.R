source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |>
  filter(time_0, time_end, frac == 1, window == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")


params <- expand_grid(min_trial = c(1, 2, 3, 4, 5, 6))

rt_test_retest <- rt_pairs |>
  group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(
    mean_var = mean(rt, na.rm = T),
    count = sum(!is.na(rt))
  ) |>
  filter(!is.na(mean_var)) |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(measure, window, time_0, time_end, during, frac, min_trial) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_test_retest, "../cached_intermediates/4_rt_test_retest.rds")
