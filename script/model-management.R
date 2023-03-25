# The model-management.R file is intended to be a scratchpad for doing things
# like defining, submitting, tagging, etc. your models. There is no need to keep
# a "record in code" of these activities because they can all be reconstructed
# later via functions like `run_log()`, as demonstrated in `model-summary.Rmd`
#
# The `Model Management Demo` (rendered from the `model-management-demo.Rmd`
# file) shows code for a range of these activities at different stages in the
# modeling process. It exists purely for reference; the intent is _not_ for you
# to replicate the full narrative.
# https://ghe.metrumrg.com/pages/example-projects/bbr-nonmem-poppk-foce/model-management-demo
#
# This script assumes you have already installed and set up bbi. For details
# on getting set up with bbr, see:
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#setup


library(bbr)
library(tidyverse)
library(here)
library(glue)

source(here("script/functions-model.R"))


# define model dir and load tags
MODEL_DIR <- here("model/nonmem")
# TAGS <- yaml::read_yaml("tags.yaml")


############################################
# CHECK FOR BBI INSTALLATION AND CONFIG
############################################

# The code below checks that you have everything configured to begin modeling
# with bbr. This code can be deleted once everything is working correctly.

# To check if you have bbi installed and configured, run `bbi_version()`
bbi_version()
# If this doesn't return a version number, see
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#installing-bbi
# for installation details.


# The first time you are modeling in a new directory, you will need to "intialize" bbi.
# The bbi_init() function will create a bbi.yaml, with the default settings, in the
# specified directory.
#
# To check you have initialized bbi, try to read in the `bbi.yaml` file.
file.path(MODEL_DIR, "bbi.yaml") %>%
  yaml::read_yaml() %>%
  names()
# If this errors, run `bbi_init()`:

bbi_init(
  .dir = MODEL_DIR, # the directory to create the bbi.yaml in
  .nonmem_dir = "/opt/NONMEM", # location of NONMEM installation
  .nonmem_version = "nm75"
) # default NONMEM version to use

# Note this only needs to be done _once for each folder_ you are modeling in. Once the bbi.yaml exists,
# you will not need to run `bbi_init()` again unless you want to create another one; for example if you
# move to modeling in a different directory.
#
# For more details on the `bbi.yaml` file and its usage, see:
# https://metrumresearchgroup.github.io/bbr/articles/getting-started.html#bbi-yaml-configuration-file

# Create new model object for FOCE model (101)
mod101 <- new_model(file.path(MODEL_DIR, 101))

# Create new Bayes model (200) from FOCE model (101)
mod <- copy_model_from(
  .parent_mod = mod101,
  .new_model = 200
)

run_number <- 200
mod <- read_model(file.path(MODEL_DIR, run_number))
bbr::model_diff(mod, .viewer = TRUE)

submit_model(
  mod,
  .bbi_args = list(overwrite = TRUE),
  .mode = "local"
)

read.table(file.path(MODEL_DIR, glue("{run_number}/{run_number}.chn")), header = TRUE) %>%
  select(ITERATION:`OMEGA.3.3.`) %>%
  pivot_longer(cols = -ITERATION) %>%
  pivot_wider(names_from = "ITERATION") %>%
  print(n = Inf)

run_chains(MODEL_DIR, run_number,
  .bbi_args = list(
    overwrite = TRUE, parallel = TRUE, threads = 4
    # overwrite = TRUE, parallel = FALSE
  ),
  # .mode = "local"
  .mode = "sge"
)

source("script/functions-diagnostics-npde.R")

theme_set(theme_bw())
ext <- get_chain_files(MODEL_DIR, run_number, .n_chains = 4, .ext = "ext")
p <- ext %>% 
  mutate(chain = factor(chain)) %>% 
  filter(ITERATION > -1e9) %>% 
  filter(ITERATION > 0) %>%
  # thinning: take every nth sample
  filter(row_number() %% 1 == 0) %>% 
  pivot_longer(-c(ITERATION, chain)) %>% 
  group_by(name) %>% 
  filter(max(value) > min(value)) %>% 
  ungroup() %>% 
  trace_plot()
p[[1]]

p[[2]]
p[[3]]
