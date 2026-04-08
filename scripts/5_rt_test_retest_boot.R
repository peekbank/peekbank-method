source("../helper/common.R")
source("../helper/rt_helper.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")
rts <- readRDS("../cached_intermediates/3_rts.rds") |> filter(time_0, time_end, frac == 1, window == 400)

d_rt_dt_long <- preprocess_rt_dt(rts) |> filter(measure == "log_land_rt")

pairs_long <- make_test_retest_pairs(d_aoi)

rt_pairs <- pairs_long |>
  left_join(d_rt_dt_long |> select(
    administration_id, rt, dataset_name
  ), relationship = "many-to-many")

rm(d_aoi, rts, d_rt_dt_long, pairs_long)
gc()

downsample_rt_test_retest <- function(start_point, sample_down, iters) {
  d <- rt_pairs |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num) |>
    mutate(count = sum(!is.na(rt))) |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iters)) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, iteration) |>
    slice_sample(n = sample_down) |>
    group_by(administration_id, dataset_name, subject_id, pair_number, session_num, iteration) |>
    summarize(mean_var = mean(rt, na.rm = T), .groups = "drop")

  wide_data <- d |>
    select(-administration_id) |>
    filter(!is.na(mean_var)) |>
    group_by(dataset_name, subject_id, pair_number, iteration) |>
    mutate(count = n()) |>
    filter(count == 2) |>
    select(-count)

  wide_data |>
    pivot_wider(names_from = session_num, values_from = mean_var) |>
    group_by(dataset_name, iteration) |>
    nest() |>
    mutate(corr = map_dbl(data, cor_test_retest_wide)) |>
    select(-data) |>
    empirical_ci()
}

params <- expand_grid(
  start_point = c(3, 5, 7, 10, 15),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
) |> filter(sample_down <= start_point)


rt_bootstrap_test_retest <- params |>
  mutate(iters = 1000) |>
  mutate(corr = pmap(list(start_point, sample_down, iters), \(s_p, s_d, iters) downsample_rt_test_retest(s_p, s_d, iters))) |>
  unnest(corr)

saveRDS(rt_bootstrap_test_retest, "../cached_intermediates/5_rt_test_retest_boot.rds")
