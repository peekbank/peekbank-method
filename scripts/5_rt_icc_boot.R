source("../helper/common.R")
source("../helper/rt_helper.R")
source("../helper/params.R")

# d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, min_rt == 400)


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
    summarize_icc_resamples("rt")
}

params <- rt_downsample_params |> filter(sample_down > 1)


rt_iccs <- params |>
  mutate(iters = 1000) |>
  mutate(corr = pmap(list(start_point, sample_down, iters), \(s_p, s_d, iters) downsample_summarize_rt(s_p, s_d, iters))) |>
  unnest(corr)

saveRDS(rt_iccs, "../cached_intermediates/5_rt_icc_boot.rds")
