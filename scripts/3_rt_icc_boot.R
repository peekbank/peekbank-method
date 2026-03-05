source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- get_age_bin_cutoff(d_aoi)

rts <- readRDS("../cached_intermediates/3_rts.rds")

d_rt_dt <- preprocess_rt_dt(rts) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label) |>
  mutate(repetition = row_number())

d_rt_dt_byage <- preprocess_rt_dt(rts) |>
  inner_join(age_bin_cutoff) |>
  group_by(dataset_name, time_0, window, time_end, during, frac, administration_id, target_label, age_bin) |>
  mutate(repetition = row_number())

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "agreement"),
  copy_names = c("bootstrap_icc")
)

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
