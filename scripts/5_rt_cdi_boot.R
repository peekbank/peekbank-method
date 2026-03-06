source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")



cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "do_cdi", "cdi_data", "boot_cdi")
)

params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)

rt_boot_cdi <- preprocess_rt_dt(rts) |>
  cross_join(params) |>
  group_by(dataset_name, administration_id, time_0, window, time_end, during, frac, measure, start_point, sample_down) |>
  mutate(count = sum(!is.na(rt))) |>
  filter(count >= start_point) |>
  select(-count) |>
  slice_sample(n = sample_down) |>
  group_by(dataset_name, administration_id, time_0, window, time_end, during, frac, measure, start_point, sample_down) |>
  summarize(mean_var = mean(rt, na.rm = T), .groups = "drop") |>
  filter(!is.na(mean_var)) |>
  group_by(time_0, window, time_end, during, frac, measure, start_point, sample_down) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/5_rt_cdi_boot.rds")
