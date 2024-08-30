## Looking at whether initial LWL accuracy and processing speed predict later CDI scores (first vs second administrations)
# vocab_size dataload
source(here::here("helper/common.R"))
# load(here::here("cached_intermediates","1_d_trial.Rds")) # ToDo: use Mike's new by-trial df..but where is it?
load(here::here("cached_intermediates","1_cdi_subjects.Rds"))
d_trial <- readRDS(here("..", "peekbank-development", "cached_intermediates", "1_d_trial.Rds"))
d_sub <- readRDS(here("..", "peekbank-development", "cached_intermediates", "1_d_sub.Rds"))

# get IRT item difficulties
en_prod_pars <- read_csv(here("aux_data","EN_production_2PL_params_slopeint.csv")) |> mutate(d = -d)
en_comp_pars <- read_csv(here("aux_data","EN_comprehension_2PL_params_slopeint.csv")) |> mutate(d = -d)
sp_prod_pars <- read_csv(here("aux_data","SP_production_2PL_params_slopeint.csv")) |> mutate(d = -d)
sp_comp_pars <- read_csv(here("aux_data","SP_comprehension_2PL_params_slopeint.csv")) |> mutate(d = -d)
theme_set(theme_bw())

# mismatched items (ToDo: correct distractor_label, too, if desired)
d_trial <- d_trial |> 
  mutate(target_cdi_label = case_when(target_label=="fish" ~ "fish (animal)",
                                      target_label=="chicken" ~ "chicken (animal)",
                                      target_label=="bike" ~ "bicycle",
                                      target_label=="carrot" ~ "carrots",
                                      target_label=="blocks" ~ "block",
                                      target_label=="slide" ~ "slide (object)",
                                      target_label=="glove" ~ "gloves",
                                      target_label=="dress" ~ "dress (object)",
                                      target_label=="can" ~ "can (object)",
                                      target_label=="boot" ~ "boots",
                                      target_label=="toy" ~ "toy (object)",
                                      .default = target_label)) |>
  mutate(distractor_cdi_label = case_when(distractor_label=="fish" ~ "fish (animal)",
                                      distractor_label=="chicken" ~ "chicken (animal)",
                                      distractor_label=="bike" ~ "bicycle",
                                      distractor_label=="carrot" ~ "carrots",
                                      distractor_label=="blocks" ~ "block",
                                      distractor_label=="slide" ~ "slide (object)",
                                      distractor_label=="glove" ~ "gloves",
                                      distractor_label=="dress" ~ "dress (object)",
                                      distractor_label=="can" ~ "can (object)",
                                      distractor_label=="boot" ~ "boots",
                                      distractor_label=="toy" ~ "toy (object)",
                                      .default = distractor_label))

trial_acc <- d_trial |>
  left_join(en_prod_pars |> select(definition, d), by=c("target_cdi_label"="definition")) |> # FixMe: should use Spanish and comprehension params in relevant datasets...
  rename(target_label_difficulty = d) |>
  left_join(en_prod_pars |> select(definition, d), by=c("distractor_cdi_label"="definition")) |>
  rename(distractor_label_difficulty = d) |>
  mutate(target_weighted_acc = target_label_difficulty * long_window_accuracy,
         weighted_acc = (target_label_difficulty+distractor_label_difficulty) * long_window_accuracy)

child_acc <- trial_acc |>
  group_by(dataset_name, subject_id, administration_id, age) |>
  summarize(LWL_accuracy = mean(long_window_accuracy, na.rm=T),
            LWL_weighted_acc = mean(weighted_acc, na.rm=T),
            LWL_target_weighted_acc = mean(target_weighted_acc, na.rm=T),
            mean_target_difficulty = mean(target_label_difficulty, na.rm=T),
            sd = sd(long_window_accuracy, na.rm=T),
            n = n())

by_subj <- trial_acc |>
  # left_join(cdi_data) |> # many-to-many
  group_by(administration_id, age) |>
  summarise(target_label_difficulty = mean(target_label_difficulty),
            distractor_label_difficulty = mean(distractor_label_difficulty),
            LWL_accuracy = mean(short_window_accuracy, na.rm=T), 
            #CDI_percent = first(CDI_percent),
            RT = mean(log(rt), na.rm=T),
            n = n()) |>
  ungroup()


cdi_acc <- cdi_data |>
  left_join(child_acc)

cdi_wide <- cdi_acc %>%
  group_by(subject_id, measure) %>%
  arrange(age) %>%
  distinct(dataset_name, subject_id, measure, age, CDI_percent) %>%
  mutate(rccdi = seq_len(n())) %>%
  pivot_wider(
    names_from = rccdi,
    values_from = c(CDI_percent, age)
  )

lwl_wide <- by_subj %>%
  left_join(
    child_acc %>%
      select(subject_id, administration_id),
    by = "administration_id"
  ) %>%
  group_by(subject_id) %>%
  arrange(age) %>%
  distinct(subject_id, LWL_accuracy, age, RT) %>%
  mutate(rclwl = seq_len(n())) %>%
  pivot_wider(
    names_from = rclwl,
    values_from = c(LWL_accuracy, age, RT)
  )
longitudinal <- merge(lwl_wide, cdi_wide, by = "subject_id", all = TRUE)


df <- longitudinal %>%
  filter(measure == "prod") %>%
  mutate(
    age_diff = age_2.y - age_1.x,
    lwl_diff = LWL_accuracy_2 - LWL_accuracy_1,
    rt_diff = RT_2 - RT_1,
    cdi_diff = CDI_percent_2 - CDI_percent_1
  ) %>%
  select(
    dataset_name,
    LWL_accuracy_1, RT_1,
    CDI_percent_1, CDI_percent_2, age_2.y
  )
df <- df[complete.cases(df), ]
table(df$dataset_name)
m0 <- df %>%
  lm(data = ., CDI_percent_2 ~ LWL_accuracy_1 + RT_1 + CDI_percent_1)
summary(m0)
gvlma::gvlma(m0)
m1 <- df %>%
  lm(
    data = .,
    CDI_percent_2 ~ RT_1 + CDI_percent_1 + age_2.y
  )
summary(m1)
gvlma::gvlma(m1)
anova(m0, m1)

#slide plot
c1 <- cor.test(df$LWL_accuracy_1, df$CDI_percent_2)$p.value
c2 <- cor.test(df$RT_1, df$CDI_percent_2)$p.value
p1 <- ggplot(data = df, aes(LWL_accuracy_1, CDI_percent_2)) +
  geom_point() +
  geom_smooth(method = "loess") + 
  annotate("text", x=0.4, y= 1.2, label= paste0("p = ", round(c1, 4))) +
  labs(x = "mean(short_window_accuracy_t1)" )
p2 <- ggplot(data = df, aes(RT_1, CDI_percent_2)) +
  geom_point() +
  geom_smooth(method = "loess") + 
  annotate("text", x=5.8, y=1.2, label= paste0("p = ", round(c2, 4))) +
  labs(x = "mean(log(RT_t1))" )
gridExtra::grid.arrange(p1, p2, ncol=2)
