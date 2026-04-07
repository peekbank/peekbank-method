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

run_icc <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(est = map_dbl(data, \(x) get_icc(x, "accuracy"))) |>
    select(-data)
}

acc_params <- expand_grid(
  t_start = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000),
  t_end = c(2000, 3000, 4000),
)

accs_summarized <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) summarize_accuracy(d_aoi, t_s, t_e)))

accs_summarized_age <- acc_params |>
  mutate(summary_data = pmap(list(t_start, t_end), \(t_s, t_e) summarize_accuracy(d_aoi_age, t_s, t_e)))

accs_icc <- accs_summarized |>
  mutate(icc = map(summary_data, run_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_icc, "../cached_intermediates/1_acc_icc.rds")

accs_icc_age <- accs_summarized_age |>
  mutate(icc = map(summary_data, run_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(accs_icc_age, "../cached_intermediates/1_acc_icc_byage.rds")
