source("../helper/common.R")

d_aoi <- readRDS("../cached_intermediates/0_d_aoi.rds")

trial_flags <- d_aoi |>
  group_by(
    dataset_name, trial_id, dataset_id, subject_id, administration_id,
    target_label
  ) |>
  summarise(
    total_target_prop = mean(correct, na.rm = TRUE),
    pre_looking = mean(correct[t_norm < 400], na.rm = TRUE),
    .groups = "drop"
  )

d_aoi_age <- make_age_bins(d_aoi)

cdi_data <- readRDS("../cached_intermediates/0_cdi_subjects.rds")


acc_trial_cdi_summarize <- function(d, flags, t_start, t_end, exclude_less_than, look_both, min_trial) {
  if (look_both == "ever") {
    flags <- filter(flags, total_target_prop > 0, total_target_prop < 1)
  } else if (look_both == "before") {
    flags <- filter(flags, pre_looking > 0, pre_looking < 1)
  }


  join_cols <- c("dataset_name", "trial_id", "dataset_id", "administration_id", "target_label")

  d |>
    semi_join(flags, by = join_cols) |>
    filter(t_norm > t_start, t_norm < t_end) |>
    group_by(across(all_of(c(
      "dataset_name", "dataset_id", "administration_id", "target_label", "trial_id"
    ))), across(any_of("age_bin"))) |>
    summarise(
      accuracy = mean(correct, na.rm = TRUE),
      prop_data = mean(!is.na(correct)),
      .groups = "drop"
    ) |>
    filter(prop_data >= exclude_less_than) |>
    filter(!is.na(accuracy)) |>
    group_by(administration_id, dataset_name, across(any_of("age_bin"))) |>
    summarize(
      mean_var = mean(accuracy, na.rm = T),
      count = n(), .groups = "drop"
    ) |>
    filter(count >= min_trial) |>
    select(-count) |>
    filter(!is.na(mean_var)) |>
    left_join(cdi_data)
}

acc_params <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
  look_both = c("before", "ever", "no_need"),
  min_trial = c(1) # no restrictions
  # side_bias =c(1.5, 1, .95, .9), #1.5 is equal to none
  # min_trials=c(1,2,3,4,5,8,10, 15)
)


acc_params_kid <- expand_grid(
  t_start = c(400),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .2, .5),
  look_both = c("no_need"),
  min_trial = c(1, 2, 3, 4, 5, 8, 10, 15)
)

accs_cdi_summarized <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) acc_trial_cdi_summarize(d_aoi, trial_flags, t_s, t_e, e, l, m)
  ))

accs_cdi_summarized_age <- acc_params |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) acc_trial_cdi_summarize(d_aoi_age, trial_flags, t_s, t_e, e, l, m)
  ))

kid_accs_cdi_summarized <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) acc_trial_cdi_summarize(d_aoi, trial_flags, t_s, t_e, e, l, m)
  ))

kid_accs_cdi_summarized_age <- acc_params_kid |>
  mutate(summary_data = pmap(
    list(t_start, t_end, exclude_less_than, look_both, min_trial),
    \(t_s, t_e, e, l, m) acc_trial_cdi_summarize(d_aoi_age, trial_flags, t_s, t_e, e, l, m)
  ))

rm(d_aoi, d_aoi_age, trial_flags)
gc()

cluster <- setup_cluster(
  libs = c("dplyr", "stringr", "purrr", "tidyr", "stats", "tibble", "boot"),
  copy_names = c("safe_boot_ci", "safe_cor", "do_cdi", "boot_cdi")
)

accs_boot_cdi <- accs_cdi_summarized |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, boot_cdi)) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_boot_cdi, "../cached_intermediates/4_acc_trial_cdi_boot.rds")

accs_boot_cdi_byage <- accs_cdi_summarized_age |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, \(d) boot_cdi(d, by_age = TRUE))) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_boot_cdi_byage, "../cached_intermediates/4_acc_trial_cdi_boot_byage.rds")


kid_accs_boot_cdi <- kid_accs_cdi_summarized |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, boot_cdi)) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_boot_cdi, "../cached_intermediates/4_acc_kid_cdi_boot.rds")

kid_accs_boot_cdi_byage <- kid_accs_cdi_summarized_age |>
  partition(cluster) |>
  mutate(cdi = map(summary_data, \(d) boot_cdi(d, by_age = TRUE))) |>
  collect() |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_boot_cdi_byage, "../cached_intermediates/4_acc_kid_cdi_boot_byage.rds")
