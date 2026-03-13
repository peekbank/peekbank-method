source("../helper/common.R")
source("../helper/rt_helper.R")

# d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)


d_rt_dt <- preprocess_rt_dt(rts) |> filter(measure == "log_land_rt")

downsample_summarize_rt <- function(start_point, sample_down, iterations) {
  d_rt_dt |>
    group_by(dataset_name, administration_id) |>
    mutate(count = sum(!is.na(rt))) |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iterations)) |>
    group_by(dataset_name, administration_id, iteration) |>
    slice_sample(n = sample_down) |>
    group_by(dataset_name, administration_id, target_label, iteration) |>
    mutate(repetition = row_number()) |>
    ungroup() |>
    run_icc()
}



run_icc <- function(d) {
  d |>
    group_by(dataset_name, iteration) |>
    nest() |>
    mutate(corr = map(data, \(x) get_icc(x, "rt"))) |>
    select(-data) |>
    unnest(corr) |>
    ungroup() |>
    empirical_ci()
}

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "agreement"),
  copy_names = c("d_rt_dt", "run_icc", "get_icc", "empirical_ci", "downsample_summarize_rt")
)


params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_iccs <- params |>
  mutate(iters = 1000) |>
  partition(cluster) |>
  mutate(corr = pmap(list(start_point, sample_down, iters), \(s_p, s_d, iters) downsample_summarize_rt(s_p, s_d, iters))) |>
  collect() |>
  unnest(corr)

saveRDS(rt_iccs, "../cached_intermediates/5_rt_icc_boot.rds")
