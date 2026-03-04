source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

age_bin_cutoff <- d_aoi |>
  filter(!is.na(correct)) |>
  distinct(administration_id, age, dataset_name) |>
  mutate(age_bin = case_when(
    age < 18 ~ "<18",
    age < 24 ~ "18-24",
    age < 36 ~ "24-36",
    age >= 36 ~ ">=36"
  )) |>
  group_by(dataset_name, age_bin) |>
  mutate(count = n()) |>
  filter(count >= 5) |>
  ungroup()

d_aoi_age <- d_aoi |> inner_join(age_bin_cutoff)

baseline_lengths <- d_aoi |>
  group_by(dataset_name, trial_id) |>
  summarise(t_min = min(t_norm))

d_aoi_bc <- d_aoi |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

d_aoi_bc_age <- d_aoi_age |>
  left_join(baseline_lengths) |>
  filter(t_min < 0)

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

# Bootstrap ICCs on pre-computed trial-level summaries.
run_bc_icc_bootstrap <- function(d) {
  d |>
    group_by(dataset_name, across(any_of("age_bin"))) |>
    nest() |>
    mutate(icc = map(data, \(x) bootstrap_icc(x, "bc_accuracy", 2000))) |>
    select(-data) |>
    unnest(icc)
}

bc_acc_params <- expand_grid(
  t_start = 400,
  t_end = c(2000, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(-500, 0),
)

# Pre-compute trial-level summaries on main process
bc_summarized <- bc_acc_params |>
  mutate(summary_data = pmap(list(b_start, b_end, t_start, t_end),
    \(b_s, b_e, t_s, t_e) summarize_bc_accuracy(d_aoi_bc, b_s, b_e, t_s, t_e)))

bc_summarized_age <- bc_acc_params |>
  mutate(summary_data = pmap(list(b_start, b_end, t_start, t_end),
    \(b_s, b_e, t_s, t_e) summarize_bc_accuracy(d_aoi_bc_age, b_s, b_e, t_s, t_e)))

rm(d_aoi, d_aoi_age, d_aoi_bc, d_aoi_bc_age, age_bin_cutoff, baseline_lengths)
gc()

# Only bootstrap_icc and run_bc_icc_bootstrap go to workers (NOT d_aoi_bc)
cluster <- new_cluster(16)
cluster_library(cluster, "dplyr")
cluster_library(cluster, "tidyr")
cluster_library(cluster, "purrr")
cluster_library(cluster, "agreement")
cluster_copy(cluster, "bootstrap_icc")
cluster_copy(cluster, "run_bc_icc_bootstrap")

bc_accs <- bc_summarized |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_bc_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(bc_accs, "../cached_intermediates/2_bc_icc_boot.rds")

bc_accs_age <- bc_summarized_age |>
  partition(cluster) |>
  mutate(icc = map(summary_data, run_bc_icc_bootstrap)) |>
  collect() |>
  select(-summary_data) |>
  unnest(col = icc)

saveRDS(bc_accs_age, "../cached_intermediates/2_bc_icc_boot_byage.rds")
