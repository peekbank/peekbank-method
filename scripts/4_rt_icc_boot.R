source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

d_rt_dt <- preprocess_rt_dt(rts) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, measure) |>
  mutate(repetition = row_number())

d_rt_dt_byage <- preprocess_rt_dt(rts) |>
  inner_join(age_bin_cutoff) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, age_bin, measure) |>
  mutate(repetition = row_number())

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "agreement"),
  copy_names = c("bootstrap_icc")
)


params <- expand_grid(min_trial = c(1, 2, 3, 4, 5, 6))

rt_iccs <- d_rt_dt |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure) |>
  mutate(count = sum(!is.na(rt))) |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, measure, min_trial) |>
  nest() |>
  partition(cluster) |>
  mutate(icc_admin = map(data, \(d) {
    bootstrap_icc(d, column = "rt", bootstrap = 2000)
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin)

saveRDS(rt_iccs, "../cached_intermediates/4_rt_icc_boot.rds")

rt_iccs_age <- d_rt_dt_byage |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, age_bin, measure) |>
  mutate(count = sum(!is.na(rt))) |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, age_bin, measure, min_trial) |>
  nest() |>
  partition(cluster) |>
  mutate(icc_admin = map(data, \(d) {
    bootstrap_icc(d, column = "rt", bootstrap = 2000)
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin)
saveRDS(rt_iccs_age, "../cached_intermediates/4_rt_icc_boot_byage.rds")
