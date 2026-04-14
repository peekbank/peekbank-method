source("../helper/common.R")
source("../helper/rt_helper.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |>
  filter(time_0, time_end, frac == 1, min_rt == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, min_rt, max_rt, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")


params <- rt_params_min_trial

rt_test_retest <- rt_pairs |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(
    mean_var = mean(rt, na.rm = T),
    count = sum(!is.na(rt))
  ) |>
  filter(!is.na(mean_var)) |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(measure, min_rt, max_rt, time_0, time_end, during, frac, min_trial) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_test_retest, "../cached_intermediates/4_rt_test_retest.rds")
