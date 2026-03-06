source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

# Summarize timepoint-level data to trial-level accuracy.
# Includes age_bin in grouping if present in the data.
downsample_summarize_accuracy <- function(d, t_start, t_end, start_point, sample_down) {
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
    slice_sample(n = sample_down) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

# Bootstrap ICCs on pre-computed trial-level summaries.
run_icc_bootstrap <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(icc = map(data, \(x) bootstrap_icc(x, "accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}

params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)


# Pre-compute trial-level summaries on main process
accs_summarized <- params |>
  mutate(summary_data = pmap(list(t_start, t_end, start_point, sample_down), \(t_s, t_e, s_p, s_d) downsample_summarize_accuracy(d_aoi, t_s, t_e, s_p, s_d)))


cluster <- setup_cluster(
  libs = c("dplyr", "tidyr", "purrr", "agreement"),
  copy_names = c("bootstrap_icc", "run_icc_bootstrap")
)

accs_boot <- accs_summarized |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_boot, "../cached_intermediates/5_acc_icc_boot.rds")
