# source me for all analyses 

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(ggpmisc))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
# remotes::install_github("jmgirard/agreement")
library(agreement)

# Seed for random number generation
set.seed(42)
knitr::opts_chunk$set(cache.extra = knitr::rand_seed, cache = TRUE, 
                      message=FALSE, warning=FALSE, error=FALSE)
options(dplyr.summarise.inform = FALSE)

# load data
load(file = here("cached_intermediates","1_d_trial.Rds"))
