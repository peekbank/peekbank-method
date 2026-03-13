source("../helper/common.R")
source("../helper/rt_helper.R")

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)
d_rt_dt <- preprocess_rt_dt(rts) |> filter(measure == "log_land_rt")

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


downsample_rt_cdi <- function(start_point, sample_down, iters) {
  d_rt_dt |>
    group_by(dataset_name, administration_id) |>
    mutate(count = sum(!is.na(rt))) |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iters)) |>
    group_by(administration_id, dataset_name, iteration) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name, iteration) |>
    summarize(mean_var = mean(rt, na.rm = T)) |>
    filter(!is.na(mean_var)) |>
    ungroup() |>
    left_join(cdi_data, by = c("administration_id", "dataset_name"))
}

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_cor", "do_cdi", "empirical_ci_cdi")
)

params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_boot_cdi <- params |>
  mutate(iters = 1000) |>
  mutate(summary_data = pmap(list(start_point, sample_down, iters), \(s_p, s_d, iters) downsample_rt_cdi(s_p, s_d, iters))) |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, \(d) {
    suppressWarnings(d |>
      group_by(dataset_name, iteration) |>
      nest() |>
      mutate(corr = map(data, \(d) do_cdi(d, 1:nrow(d)))) |>
      select(-data) |>
      unnest_wider(corr) |>
      ungroup() |>
      empirical_ci_cdi())
  })) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/5_rt_cdi_boot.rds")
