source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")




rt_params <- expand_grid(
  option = c("recommended"),
  window = c(400),
  measure = c("log_land_rt"),
  min_trial = c(1)
) |> bind_rows(
  expand_grid(
    option = c("alternative"),
    window = c(200),
    measure = c("first_launch_rt"),
    min_trial = c(5)
  )
)

rt_icc <- rt_params |>
  left_join(preprocess_rt_dt(rts) |>
    group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, measure) |>
    mutate(repetition = row_number())) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure, option) |>
  mutate(count = sum(!is.na(rt))) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, measure, min_trial, option) |>
  nest() |>
  mutate(icc_admin = map(data, \(d) {
    get_icc(d, column = "rt")
  })) |>
  select(-data) |>
  unnest(icc_admin)

saveRDS(rt_icc, "../cached_intermediates/8_rt_icc.rds")

rt_cdi <- rt_params |>
  left_join(preprocess_rt_dt(rts) |>
    group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure) |>
    summarize(
      mean_var = mean(rt, na.rm = T),
      count = sum(!is.na(rt))
    ) |>
    filter(!is.na(mean_var))) |>
  filter(count >= min_trial) |>
  select(-count) |>
  left_join(cdi_data) |>
  group_by(time_0, window, time_end, during, frac, measure, min_trial, option) |>
  nest() |>
  mutate(cdi = map(data, calc_cdi)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_cdi, "../cached_intermediates/8_rt_cdi.rds")


rt_trt <- rt_params |>
  left_join(make_test_retest_pairs(d_aoi) |>
    left_join(preprocess_rt_dt(rts) |> select(
      administration_id, rt, measure, window, time_0, time_end, during,
      frac, dataset_name
    ), relationship = "many-to-many")) |>
  group_by(min_trial, measure, window, time_0, time_end, during, frac, dataset_name, administration_id, subject_id, pair_number, session_num, option) |>
  summarize(
    mean_var = mean(rt, na.rm = T),
    count = sum(!is.na(rt))
  ) |>
  filter(!is.na(mean_var)) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(measure, window, time_0, time_end, during, frac, min_trial, option) |>
  nest() |>
  mutate(cdi = map(data, calc_test_retest)) |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_trt, "../cached_intermediates/8_rt_test_retest.rds")
