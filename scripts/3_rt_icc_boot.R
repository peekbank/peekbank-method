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
  select(-rt, -shift_start_rt, -last_shift_rt) | |> 
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label) |>
  mutate(repetition = row_number())

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
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, age_bin) |>
  mutate(repetition = row_number())

cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "bootstrap_icc")

rt_iccs <- d_rt_dt |>
  group_by(dataset_name, time_0, window, time_end, during, frac) |>
  nest() |>
  partition(cluster) |>
  mutate(icc_admin = map(data, \(d) {
    rt_names <- colnames(d)[str_ends(colnames(d), "rt")]
    n_total <- map_dbl(rt_names, \(rt_col) {
      d |>
        filter(!is.na(d[rt_col])) |>
        nrow()
    })
    n_admin <- map_dbl(rt_names, \(rt_col) {
      d |>
        filter(!is.na(d[rt_col])) |>
        select(administration_id) |>
        unique() |>
        nrow()
    })
    icc_values <- map(rt_names, \(rt_col) {
      bootstrap_icc(d |> filter(!is.na(d[rt_col])), column = rt_col, bootstrap = 2000)
    })
    tibble(
      measure = rt_names, foo = icc_values,
      n_total = n_total, n_admin = n_admin
    )
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin) |>
  unnest(foo)

saveRDS(rt_iccs, "../cached_intermediates/3_rt_icc_boot.rds")

rt_iccs_age <- d_rt_dt_byage |>
  group_by(dataset_name, time_0, window, time_end, during, frac, age_bin) |>
  nest() |>
  partition(cluster) |>
  mutate(icc_admin = map(data, \(d) {
    rt_names <- colnames(d)[str_ends(colnames(d), "rt")]
    n_total <- map_dbl(rt_names, \(rt_col) {
      d |>
        filter(!is.na(d[rt_col])) |>
        nrow()
    })
    n_admin <- map_dbl(rt_names, \(rt_col) {
      d |>
        filter(!is.na(d[rt_col])) |>
        select(administration_id) |>
        unique() |>
        nrow()
    })
    icc_values <- map(rt_names, \(rt_col) {
      bootstrap_icc(d |> filter(!is.na(d[rt_col])), column = rt_col, bootstrap = 2000)
    })
    tibble(
      measure = rt_names, foo = icc_values,
      n_total = n_total, n_admin = n_admin
    )
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin) |>
  unnest(foo)

saveRDS(rt_iccs_age, "../cached_intermediates/3_rt_icc_boot_byage.rds")
