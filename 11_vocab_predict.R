## Looking at whether initial LWL accuracy predicts later CDI (first vs second administrations)
# vocab_size dataload
source(here::here("helper/common.R"))
load(here::here("cached_intermediates", "1_cdi_subjects.Rds"))
theme_set(theme_bw())
en_prod_pars <- read_csv(
  here("aux_data", "EN_production_2PL_params_slopeint.csv")
) |> mutate(d = -d)
en_comp_pars <- read_csv(
  here("aux_data", "EN_comprehension_2PL_params_slopeint.csv")
) |> mutate(d = -d)
sp_prod_pars <- read_csv(
  here("aux_data", "SP_production_2PL_params_slopeint.csv")
) |> mutate(d = -d)
sp_comp_pars <- read_csv(
  here("aux_data", "SP_comprehension_2PL_params_slopeint.csv")
) |> mutate(d = -d)

d_trial <- d_trial |> 
  mutate(target_cdi_label = case_when(
    target_label == "fish" ~ "fish (animal)",
    target_label == "chicken" ~ "chicken (animal)",
    target_label == "bike" ~ "bicycle",
    target_label == "carrot" ~ "carrots",
    target_label == "blocks" ~ "block",
    target_label == "slide" ~ "slide (object)",
    target_label == "glove" ~ "gloves",
    target_label == "dress" ~ "dress (object)",
    target_label == "can" ~ "can (object)",
    target_label == "boot" ~ "boots",
    target_label == "toy" ~ "toy (object)",
    .default = target_label
  ))

acc_start <- 300
acc_end <- 4000

trial_acc <- d_trial |>
  group_by(
    dataset_name, subject_id, administration_id, age,
    trial_id, trial_order, dataset_id, target_cdi_label, distractor_label
  ) |>
  filter(t_norm >= acc_start) |>
  filter(t_norm <= acc_end) |>
  filter(!is.na(correct)) |>
  summarize(acc = mean(correct)) |>
  left_join(
    en_prod_pars |>
      select(definition, d), by = c("target_cdi_label" = "definition")
  ) |>
  rename(target_label_difficulty = d) |>
  left_join(
    en_prod_pars |>
      select(definition, d), by = c("distractor_label"="definition")
  ) |>
  rename(distractor_label_difficulty = d) |>
  mutate(
    target_weighted_acc = target_label_difficulty * acc,
    weighted_acc = (target_label_difficulty+distractor_label_difficulty) * acc
  )

child_acc <- trial_acc |>
  group_by(dataset_name, subject_id, administration_id, age) |>
  summarize(
    LWL_accuracy = mean(acc),
    LWL_weighted_acc = mean(weighted_acc),
    LWL_target_weighted_acc = mean(target_weighted_acc),
    mean_target_difficulty = mean(target_label_difficulty),
    sd = sd(acc),
    n = n()
  )

# Testing growth
colnames(cdi_data)
colnames(child_acc)
demographics_nmisc <- cdi_data %>%
  select(
    dataset_name, subject_id, sex,
    native_language, language, instrument_type
  ) %>%
  distinct
cdi_wide <- cdi_data %>%
  #filter(lab_subject_id!="L50") %>%
  group_by(subject_id, measure) %>%
  arrange(age) %>%
  distinct(subject_id, measure, age, CDI_percent) %>%
  mutate(rccdi = seq_len(n())) %>%
  pivot_wider(
    names_from = rccdi,
    values_from = c(CDI_percent, age)
  )
lwl_wide <- child_acc %>%
  group_by(subject_id) %>%
  arrange(age) %>%
  distinct(subject_id, LWL_accuracy, age) %>%
  mutate(rclwl = seq_len(n())) %>%
  pivot_wider(
    names_from = rclwl,
    values_from = c(LWL_accuracy, age)
  )
longitudinal <- left_join(
  merge(lwl_wide, cdi_wide, by = "subject_id", all = TRUE),
  demographics_nmisc,
  by = "subject_id",
  multiple = "any"
)
