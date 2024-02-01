# source me for all analyses 

library(here)
library(tidyverse)
library(multidplyr)
# library(ggpmisc)
library(ggrepel)
# library(ggthemes)
library(viridis)
# library(cowplot)
# remotes::install_github("jmgirard/agreement")
library(agreement)
library(tictoc)

# Seed for random number generation
set.seed(42)
knitr::opts_chunk$set(cache.extra = knitr::rand_seed, cache = TRUE, 
                      message=FALSE, warning=FALSE, error=FALSE)
options(dplyr.summarise.inform = FALSE)

# load data
load(file = here("cached_intermediates","1_d_trial.Rds"))

# key ICC function
get_icc <- function (x, column = "accuracy", object = "stimulus") {
  if (object == "stimulus") {
    iccs <- dim_icc(x, 
                    model = "2A", 
                    type = "consistency", 
                    unit = "average",
                    object = target_label, 
                    rater = administration_id,
                    trial = trial_id, 
                    score = {{column}}, 
                    bootstrap = 0)
  } else {
    iccs <- dim_icc(x, 
                    model = "2A", 
                    type = "consistency", 
                    unit = "average",
                    object = administration_id, 
                    rater = target_label,
                    trial = trial_id, 
                    score = {{column}}, 
                    bootstrap = 0)
  }
  
  return(iccs$Inter_ICC)
}