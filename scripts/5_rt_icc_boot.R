source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

d_rt_dt <- preprocess_rt_dt(rts)

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "agreement"),
  copy_names = c("bootstrap_icc")
)


params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_iccs <- d_rt_dt |>
  cross_join(params) |>
  group_by(dataset_name, administration_id, time_0, window, time_end, during, frac, measure, start_point, sample_down) |>
  mutate(count = sum(!is.na(rt))) |>
  filter(count >= start_point) |>
  select(-count) |>
  slice_sample(n = sample_down) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, measure) |>
  mutate(repetition = row_number()) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, measure) |>
  nest() |>
  partition(cluster) |>
  mutate(icc_admin = map(data, \(d) {
    bootstrap_icc(d, column = "rt", bootstrap = 2000)
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin)

saveRDS(rt_iccs, "../cached_intermediates/5_rt_icc_boot.rds")
