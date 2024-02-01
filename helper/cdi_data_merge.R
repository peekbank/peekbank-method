# Importing CDI data for two datasets from accuracies analysis. 
# First get Garrison & Bergelson / Swingley & Aslin data. 
library(peekbankr)
subjects <- get_subjects()
sa_administrations <- get_administrations(dataset_name = "swingley_aslin_2002")
sa_trial_types <- get_trial_types(dataset_name = "swingley_aslin_2002")
sa_trials <- get_trials(dataset_name = "swingley_aslin_2002")
sa_aoi_timepoints <- get_aoi_timepoints(dataset_name = "swingley_aslin_2002")

sa_data <- sa_aoi_timepoints |>
  left_join(sa_administrations) |>
  left_join(sa_trials) |>
  left_join(sa_trial_types) |>
  left_join(subjects) |>
  filter(condition != "filler") |>
  mutate(condition = if_else(condition == "cp", "Correct", "Mispronounced"))

gb_administrations <- get_administrations(dataset_name = "garrison_bergelson_2020")
gb_trial_types <- get_trial_types(dataset_name = "garrison_bergelson_2020")
gb_trials <- get_trials(dataset_name = "garrison_bergelson_2020")
gb_aoi_timepoints <- get_aoi_timepoints(dataset_name = "garrison_bergelson_2020")

gb_data <- gb_aoi_timepoints |>
  left_join(gb_administrations) |>
  left_join(gb_trials) |>
  left_join(gb_trial_types) |>
  left_join(subjects)

cdis <- read_csv(here::here("aux_data","cdi_data.csv"))

vanilla_cdi_datasets <- bind_rows(filter(sa_data, condition == "Correct"), 
                                  gb_data) |>
  mutate(correct = ifelse(aoi == "target", 1, 
                          ifelse(aoi == "distractor", 0, NA)))
