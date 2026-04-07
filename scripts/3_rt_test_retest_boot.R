source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |>
  filter(time_0, time_end, frac == 1)

rts_weird <- readRDS("../cached_intermediates/3_rts.rds") |> filter(window == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

d_rt_dt_long <- preprocess_rt_dt(rts)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")



# rt_pairs_weird <- pairs_long |>
#   left_join(preprocess_rt_dt(rts_weird) |> filter(measure %in% c("log_land_rt", "land_rt")) |> select(
#     administration_id, rt, measure, window, time_0, time_end, during,
#     frac, dataset_name
#   ), relationship = "many-to-many")

rt_test_retest <- rt_pairs |>
  group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  group_by(measure, window, time_0, time_end, during, frac) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_test_retest, "../cached_intermediates/3_rt_test_retest.rds")

# rt_test_retest_weird <- rt_pairs_weird |>
#   group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
#   summarize(mean_var = mean(rt, na.rm = T)) |>
#   group_by(measure, window, time_0, time_end, during, frac) |>
#   nest() |>
#   mutate(cdi = map(data, calc_test_retest)) |>
#   select(-data) |>
#   unnest(cdi)
#
# saveRDS(rt_test_retest_weird, "../cached_intermediates/3_rt_test_retest_boot_weird.rds")
