source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- d_aoi |>
  filter(!is.na(correct)) |>
  distinct(administration_id, age, dataset_name) |>
  mutate(age_bin = case_when(
    age < 18 ~ "<18",
    age < 24 ~ "18-24",
    age < 36 ~ "24-36",
    age >= 36 ~ ">=36"
  )) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= 5) |>
  ungroup()

rts <- readRDS("../cached_intermediates/3_rts.rds")

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")

d_rt_dt <- rts |>
  filter(shift_type == "D-T") |>
  mutate(
    land_rt = rt,
    first_launch_rt = shift_start_rt,
  ) |>
  mutate(
    across(c("land_rt", "first_launch_rt"), log, .names = "log_{.col}"),
    across(
      c(
        "land_rt", "first_launch_rt",
        "log_land_rt", "log_first_launch_rt",
      ),
      ~ ifelse(shift_length <= 600, .x, NA),
      .names = "trim_first_{.col}"
    ),
    across(
      c(
        "land_rt",
        "log_land_rt",
      ),
      ~ ifelse(last_shift_length <= 600, .x, NA),
      .names = "trim_last_{.col}"
    )
  ) |>
  select(-rt, -shift_start_rt, -last_shift_rt) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt") |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var))

d_rt_dt_byage <- rts |>
  filter(shift_type == "D-T") |>
  mutate(
    land_rt = rt,
    first_launch_rt = shift_start_rt,
  ) |>
  mutate(
    across(c("land_rt", "first_launch_rt"), log, .names = "log_{.col}"),
    across(
      c(
        "land_rt", "first_launch_rt",
        "log_land_rt", "log_first_launch_rt",
      ),
      ~ ifelse(shift_length <= 600, .x, NA),
      .names = "trim_first_{.col}"
    ),
    across(
      c(
        "land_rt",
        "log_land_rt",
      ),
      ~ ifelse(last_shift_length <= 600, .x, NA),
      .names = "trim_last_{.col}"
    )
  ) |>
  select(-rt, -shift_start_rt, -last_shift_rt) |>
  inner_join(age_bin_cutoff) |>
  pivot_longer(ends_with("rt"), names_to = "measure", values_to = "rt") |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, age_bin, measure) |>
  summarize(mean_var = mean(rt, na.rm = T)) |>
  filter(!is.na(mean_var))

library(boot)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "stats")
cluster_library(cluster, "tibble")
cluster_library(cluster, "boot")
cluster_copy(cluster, "do_cdi")
cluster_copy(cluster, "cdi_data")
cluster_copy(cluster, "boot_cdi")
cluster_copy(cluster, "boot_cdi_age")


rt_boot_cdi <- d_rt_dt |>
  group_by(time_0, window, time_end, during, frac) |>
  nest() |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/3_rt_cdi_boot.rds")

rt_boot_cdi_byage <- d_rt_dt_byage |>
  group_by(age_bin, time_0, window, time_end, during, frac) |>
  nest() |>
  mutate(cdi = map(data, boot_cdi)) |>
  collect() |>
  select(-data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi_byage, "../cached_intermediates/3_rt_cdi_boot_byage.rds")

