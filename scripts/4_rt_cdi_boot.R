source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

d_rt_dt <- preprocess_rt_dt(rts) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure) |>
  summarize(
    mean_var = mean(rt, na.rm = T),
    count = sum(!is.na(rt))
  ) |>
  filter(!is.na(mean_var))

d_rt_dt_byage <- preprocess_rt_dt(rts) |>
  inner_join(age_bin_cutoff) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, age_bin, measure) |>
  summarize(
    mean_var = mean(rt, na.rm = T),
    count = sum(!is.na(rt))
  ) |>
  filter(!is.na(mean_var))

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "do_cdi", "boot_cdi")
)

params <- expand_grid(min_trial = c(1, 2, 3, 4, 5, 6))

rt_boot_cdi <- d_rt_dt |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  left_join(cdi_data) |> 
  group_by(time_0, window, time_end, during, frac, measure, min_trial) |>
  nest() |> 
  #partition(cluster) |>
  mutate(cdi = map(data, boot_cdi)) |>
  #collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/4_rt_cdi_boot.rds")

rt_boot_cdi_byage <- d_rt_dt_byage |>
  cross_join(params) |>
  filter(count >= min_trial) |>
  select(-count) |>
  left_join(cdi_data) |> 
  group_by(age_bin, time_0, window, time_end, during, frac, measure, min_trial) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi_byage, "../cached_intermediates/4_rt_cdi_boot_byage.rds")
