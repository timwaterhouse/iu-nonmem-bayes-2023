library(tidyverse)
library(mrgsolve)
library(here)
library(glue)
library(bbr)
theme_set(theme_bw())

modelName <- 200

mod_sim <- mread(here(glue("model/mrgsolve/{modelName}.mod")))
param(mod_sim)

# read bbr model for first chain
mod_bbr <- read_model(here(glue("model/nonmem/{modelName}/{modelName}_1")))

data <- nm_join(mod_bbr)

out <- mrgsim(
  mod_sim, 
  data      = data, 
  etasrc    = "data.all", 
  obsonly   = TRUE, 
  Req       = "MRG = IPRED",
  carry_out = "NM = IPRED, PRED"
) %>% 
  as_tibble() %>% 
  mutate(
    diff = MRG - NM,
    pct_diff = diff / NM * 100
  )

out %>% select(diff, pct_diff) %>% summary()

ggplot(out, aes(NM, MRG)) +
  geom_point() +
  geom_abline()
