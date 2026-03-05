source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

# Summarize timepoint-level data to trial-level accuracy.
# Includes age_bin in grouping if present in the data.
summarize_accuracy <- function(d, t_start, t_end) {
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

acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)

# Pre-compute trial-level summaries on main process
accs_summarized <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) summarize_accuracy(d_aoi, t_s, t_e)))

accs_summarized_age <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) summarize_accuracy(d_aoi_age, t_s, t_e)))

rm(d_aoi, d_aoi_age)
gc()

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

saveRDS(accs_boot, "../cached_intermediates/1_acc_icc_boot.rds")

accs_boot_age <- accs_summarized_age |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_boot_age, "../cached_intermediates/1_acc_icc_boot_byage.rds")
