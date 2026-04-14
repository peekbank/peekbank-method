source("../helper/common.R")
source("../helper/params.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

d_aoi_age <- make_age_bins(d_aoi)

d_aoi_bc <- make_baseline_corrected(d_aoi)
d_aoi_bc_age <- make_baseline_corrected(d_aoi_age)

# Summarize to trial-level baseline-corrected accuracy.
# Includes age_bin in grouping if present in the data.
summarize_bc_accuracy <- function(d, b_start, b_end, t_start, t_end) {
  d |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
    ))), across(any_of("age_bin"))) |>
    summarise(
      window_accuracy = mean(correct[t_norm >= t_start & t_norm <= t_end], na.rm = TRUE),
      baseline_accuracy = mean(correct[t_norm >= b_start & t_norm <= b_end], na.rm = TRUE),
      bc_accuracy = window_accuracy - baseline_accuracy,
      .groups = "drop"
    ) |>
    filter(!is.na(bc_accuracy)) |>
    group_by(across(all_of(c(
      "dataset_id", "dataset_name", "administration_id", "target_label"
    ))), across(any_of("age_bin"))) |>
    mutate(repetition = row_number()) |>
    ungroup()
}

run_bc_icc <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(est = map_dbl(data, \(x) get_icc(x, "bc_accuracy"))) |>
    select(-data)
}

bc_acc_params <- bc_params

# Pre-compute trial-level summaries on main process
bc_summarized <- bc_acc_params |>
  mutate(summary_data = pmap(
    list(b_start, b_end, t_start, t_end),
    \(b_s, b_e, t_s, t_e) summarize_bc_accuracy(d_aoi_bc, b_s, b_e, t_s, t_e)
  ))

bc_summarized_age <- bc_acc_params |>
  mutate(summary_data = pmap(
    list(b_start, b_end, t_start, t_end),
    \(b_s, b_e, t_s, t_e) summarize_bc_accuracy(d_aoi_bc_age, b_s, b_e, t_s, t_e)
  ))

rm(d_aoi, d_aoi_age, d_aoi_bc, d_aoi_bc_age)
gc()

bc_accs <- bc_summarized |>
  mutate(icc = map(summary_data, run_bc_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(bc_accs, "../cached_intermediates/2_bc_icc.rds")

bc_accs_age <- bc_summarized_age |>
  mutate(icc = map(summary_data, run_bc_icc)) |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(bc_accs_age, "../cached_intermediates/2_bc_icc_byage.rds")
