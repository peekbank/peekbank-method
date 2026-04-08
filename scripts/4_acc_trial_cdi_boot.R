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
    left_join(cdi_data, by = c("dataset_name", "administration_id"))
}

acc_params <- expand_grid(
  t_start = c(200, 400, 600),
  t_end = c(2000, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1),
  look_both = c("before", "ever", "no_need"),
  min_trial = c(1) # no restrictions
)


acc_params_kid <- expand_grid(
  t_start = c(400, 600, 200),
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

accs_cdi <- accs_cdi_summarized |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi, "../cached_intermediates/4_acc_trial_cdi.rds")

accs_cdi_byage <- accs_cdi_summarized_age |>
  mutate(cdi = map(summary_data, \(d) calc_cdi(d, by_age = TRUE))) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(accs_cdi_byage, "../cached_intermediates/4_acc_trial_cdi_byage.rds")


kid_accs_cdi <- kid_accs_cdi_summarized |>
  mutate(cdi = map(summary_data, calc_cdi)) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_cdi, "../cached_intermediates/4_acc_kid_cdi.rds")

kid_accs_cdi_byage <- kid_accs_cdi_summarized_age |>
  mutate(cdi = map(summary_data, \(d) calc_cdi(d, by_age = TRUE))) |>
  select(-summary_data) |>
  unnest(cdi)

saveRDS(kid_accs_cdi_byage, "../cached_intermediates/4_acc_kid_cdi_byage.rds")
