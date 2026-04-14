source("../helper/common.R")
source("../helper/rt_helper.R")
source("../helper/params.R")

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, min_rt == 400)
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

params <- rt_downsample_params


rt_boot_cdi <- params |>
  mutate(iters = 1000) |>
  mutate(
    cdi = pmap(
      list(start_point, sample_down, iters),
      \(s_p, s_d, iters) {
        suppressWarnings(downsample_rt_cdi(s_p, s_d, iters)) |>
          group_by(dataset_name, iteration) |>
          nest() |>
          mutate(corr = map2(data, dataset_name, \(d, dn) {
            calc_cdi(mutate(d, dataset_name = dn)) |> select(-any_of("dataset_name"))
          })) |>
          select(-data) |>
          unnest_wider(corr) |>
          ungroup() |>
          empirical_ci_cdi()
      }
    )
  ) |>
  unnest(cdi)

saveRDS(rt_boot_cdi, "../cached_intermediates/5_rt_cdi_boot.rds")
