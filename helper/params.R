# Shared parameter definitions for scripts in ../scripts.
# This file is sourced by scripts after common helpers are loaded.

acc_params <- expand_grid(
  t_start = seq(-500, 1000, 100),
  t_end = seq(1500, 4000, 100)
)


bc_params <- expand_grid(
  t_start = c(300, 400, 600),
  t_end = c(1800, 3000, 4000),
  b_start = seq(-4000, -1000, 1000),
  b_end = c(0)
)


rt_params <- expand_grid(
  min_rt = seq(0, 1000, 100),
  max_rt = seq(1600, 4000, 200),
  time_0 = c(TRUE),
  time_end = c(TRUE),
  during = TRUE,
  frac = c(1)
)


rt_params_weird <- expand_grid(
  min_rt = c(400),
  max_rt = c(4000),
  time_0 = c(FALSE, TRUE),
  time_end = c(FALSE, TRUE),
  during = TRUE,
  frac = c(0, .5, .75, 1)
)

acc_params_trial <- expand_grid(
  t_start = c(300, 400, 600),
  t_end = c(1800, 3000, 4000),
  exclude_less_than = c(0, .1, .2, .3, .4, .5, .6, .7, .8),
  look_at_start = c("yes", "no"),
  min_trial = c(1),
  min_frac = c(0),
)

acc_params_kid <- expand_grid(
  t_start = c(300, 400, 600),
  t_end = c(1800, 3000, 4000),
  exclude_less_than = c(0, .2, .5, .7),
  look_at_start = c("no"),
  min_trial = c(1, 2, 3, 4, 6, 8, 10, 12, 14, 16),
  min_frac = c(0, .2, .4, .5, .6, .8)
) |> filter(!(min_trial > 1 & min_frac > 0)) # only use one at a time



rt_params_min_trial <- expand_grid(min_trial = c(1, 2, 3, 4, 5, 6))


# for simulated downsampling
acc_downsample_params <- expand_grid(
  t_start = c(600),
  t_end = c(4000),
  start_point = c(5, 10, 15, 20),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
) |> filter(sample_down <= start_point)

rt_downsample_params <- expand_grid(
  start_point = c(3, 5, 7, 10),
  sample_down = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
) |> filter(sample_down <= start_point)


# recommended comparisons
acc_params_summary_recommended <- expand_grid(
  option = "recommended",
  t_start = c(600),
  t_end = c(4000),
  exclude_less_than = c(0),
  look_at_start = c("no"),
  min_trial = c(1),
  min_frac = c(0)
)


bc_params_summary_alternative <- expand_grid(
  option = "alternative",
  t_start = c(300),
  t_end = c(1800),
  b_start = c(-2000),
  b_end = c(0),
  exclude_less_than = c(.6),
  look_at_start = c("yes"),
  min_trial = c(1),
  min_frac = c(.5)
)

rt_params_summary_options <- expand_grid(
  option = c("recommended"),
  min_rt = c(400),
  max_rt = c(2000),
  measure = c("log_land_rt"),
  min_trial = c(1)
) |>
  bind_rows(expand_grid(
    option = c("alternative"),
    min_rt = c(300),
    max_rt = c(1800),
    measure = c("first_launch_rt"),
    min_trial = c(2)
  ))
