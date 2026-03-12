source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

# Summarize timepoint-level data to trial-level accuracy.
# Includes age_bin in grouping if present in the data.
downsample_summarize_accuracy <- function(d, t_start, t_end, start_point, sample_down, iter) {
  d |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
    ))), across(any_of("age_bin"))) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(!is.na(accuracy)) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id"
    ))), across(any_of("age_bin"))) |>
    mutate(count = n()) |>
    filter(count >= start_point) |>
    select(-count) |>
    cross_join(tibble(iteration = 1:iter)) |>
    group_by(iteration, across(all_of(c(
      "dataset_name", "dataset_id", "administration_id"
    ))), across(any_of("age_bin"))) |>
    slice_sample(n = sample_down) |>
    group_by(iteration, across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

# Bootstrap ICCs on pre-computed trial-level summaries.
run_icc <- function(d) {
  d |>
    group_by(dataset_name, iteration) |>
    nest() |>
    mutate(corr = map(data, \(x) get_icc(x, "accuracy"))) |>
    select(-data) |>
    unnest(corr) |>
    ungroup() |>
    empirical_ci()
}

params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)


# Pre-compute trial-level summaries on main process
accs_summarized <- params |>
  mutate(iters = 1000) |>
  mutate(summary_data = pmap(list(t_start, t_end, start_point, sample_down, iters), \(t_s, t_e, s_p, s_d, iters) downsample_summarize_accuracy(d_aoi, t_s, t_e, s_p, s_d, iters)))


cluster <- setup_cluster(
  libs = c("dplyr", "tidyr", "purrr", "agreement"),
  copy_names = c("bootstrap_icc", "run_icc_bootstrap", "get_icc", "empirical_ci")
)

accs_boot <- accs_summarized |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_icc)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_boot, "../cached_intermediates/5_acc_icc_boot.rds")
