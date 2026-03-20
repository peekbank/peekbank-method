source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |>
  filter(time_0, time_end, frac == 1)

rts_weird <- readRDS("../cached_intermediates/3_rts.rds") |> filter(window == 400)

pairs_long <- make_test_retest_pairs(d_aoi)

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "test_retest_corr", "boot_test_retest")
)

d_rt_dt_long <- preprocess_rt_dt(rts)

# rt_pairs <- pairs_long |>
#   left_join(d_rt_dt_long |> select(
#     administration_id, rt, measure, window, time_0, time_end, during,
#     frac, dataset_name
#   ), relationship = "many-to-many")


rt_pairs_age <- pairs_long |>
  mutate(age_bin = case_when(
    mean_age < 18 ~ "<18",
    mean_age < 24 ~ "18-24",
    mean_age < 36 ~ "24-36",
    mean_age >= 36 ~ ">=36"
  )) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= 5) |>
  ungroup() |>
  filter(age_bin %in% c("18-24", "24-36")) |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, measure, window, time_0, time_end, during,
    frac, dataset_name
  ), relationship = "many-to-many")

# rt_pairs_weird <- pairs_long |>
#   left_join(preprocess_rt_dt(rts_weird) |> filter(measure %in% c("log_land_rt", "land_rt")) |> select(
#     administration_id, rt, measure, window, time_0, time_end, during,
#     frac, dataset_name
#   ), relationship = "many-to-many")

# rt_bootstrap_test_retest <- rt_pairs |>
#   group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
#   summarize(mean_var = mean(rt, na.rm = T)) |>
#   group_by(measure, window, time_0, time_end, during, frac) |>
#   nest() |>
#   partition(cluster) |>
#   mutate(cdi = map(data, boot_test_retest)) |>
#   collect() |>
#   select(-data) |>
#   unnest(cdi)
#
# saveRDS(rt_bootstrap_test_retest, "../cached_intermediates/3_rt_test_retest_boot.rds")

# rt_bootstrap_test_retest_weird <- rt_pairs_weird |>
#   group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num) |>
#   summarize(mean_var = mean(rt, na.rm = T)) |>
#   group_by(measure, window, time_0, time_end, during, frac) |>
#   nest() |>
#   partition(cluster) |>
#   mutate(cdi = map(data, boot_test_retest)) |>
#   collect() |>
#   select(-data) |>
#   unnest(cdi)
#
# saveRDS(rt_bootstrap_test_retest_weird, "../cached_intermediates/3_rt_test_retest_boot_weird.rds")


rt_bootstrap_test_retest_age <- rt_pairs_age |>
  group_by(measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num, age_bin) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  mutate(dataset_name = str_c(age_bin, dataset_name)) |>
  group_by(measure, window, time_0, time_end, during, frac) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_test_retest)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_bootstrap_test_retest_age, "../cached_intermediates/3_rt_test_retest_boot_age.rds")
