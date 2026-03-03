source("../helper/common.R")

rts <- readRDS("../cached_intermediates/3_rts.rds")

cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "stringr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "bootstrap_icc")

rts <- readRDS("rts.rds")
d_aoi <- readRDS("d_aoi.Rds")
d_rt_dt <- rts |>
  left_join(d_aoi |> select(age, administration_id) |> distinct()) |>
  filter(window %in% c(200, 350, 375, 400, 425, 450, 500)) |>
  filter(frac %in% c(0, 1)) |>
  filter(time_0 == T, time_end == T) |>
  filter(shift_type == "D-T") |>
  mutate(
    land_rt = rt,
    first_launch_rt = shift_start_rt,
    # last_launch_rt = last_shift_rt
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
  ) |>
  select(-rt, -shift_start_rt, -last_shift_rt) |>
  mutate(younger = age < 24) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, younger) |>
  mutate(repetition = row_number())


rt_iccs_age <- d_rt_dt |>
  group_by(dataset_name, younger, time_0, window, time_end, during, frac) |>
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
    # icc_values
  })) |>
  collect() |>
  select(-data) |>
  unnest(icc_admin) |>
  unnest(foo)

saveRDS(rt_iccs_age, "4_rt_iccs_age_boot.rds")
