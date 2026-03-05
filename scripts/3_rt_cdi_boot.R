source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds")

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

d_rt_dt <- preprocess_rt_dt(rts) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt") |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var))

d_rt_dt_byage <- preprocess_rt_dt(rts) |>
  inner_join(age_bin_cutoff) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt") |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, age_bin, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var))

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "do_cdi", "cdi_data", "boot_cdi")
)


rt_boot_cdi <- d_rt_dt |>
  group_by(time_0, window, time_end, during, frac) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/3_rt_cdi_boot.rds")

rt_boot_cdi_byage <- d_rt_dt_byage |>
  group_by(age_bin, time_0, window, time_end, during, frac) |>
  nest() |>
  partition(cluster) |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi_byage, "../cached_intermediates/3_rt_cdi_boot_byage.rds")
