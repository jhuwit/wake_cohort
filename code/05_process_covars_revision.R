
library(readr)
library(dplyr)
library(tidyverse)
library(lubridate)
library(purrr)
library(stringr)
library(here)
library(slider)
library(furrr)
library(future)
library(patchwork)
# library(fuzzyjoin)
options(digits.secs = 2)
force = TRUE
# source(here::here("code", "utilities.R"))

covars = read_rds(here::here("data", "processed", "covars_proc.rds"))

icu = read_csv(here::here("data", "raw", "icu.csv"),
               col_types = cols(id = col_character())) %>% 
  janitor::clean_names()
cr = read_csv(here::here("data", "raw", "cr.csv"),
              col_types = cols(ID = col_character())) %>% 
  janitor::clean_names()
los = read_csv(here::here("data", "raw", "los.csv"),
               col_types = cols(ID = col_character())) %>% 
  janitor::clean_names() %>% 
  rename(bin_hospmort = in_hospmor) %>% 
  mutate(bin_hospmort = if_else(bin_hospmort =="N", 0, 1))

covars_updated = 
  covars %>% 
  left_join(icu, by = "id")  %>% 
  left_join(cr, by = "id") %>% 
  left_join(los, by = "id") %>% 
  mutate(increase = max_of_cr / val_creatlst,
         change = max_of_cr - val_creatlst,
         bin_aki7d = if_else(increase >= 1.5 | bin_aki == 1, 1, 0),
         aki_stage = case_when(bin_aki == 0 & bin_aki7d == 0 ~ 0,
                           increase >= 3 | max_of_cr >= 4 ~ 3,
                           increase >= 2 & increase < 3 ~ 2,
                           .default = 1))

write_rds(covars_updated, here::here("data", "processed", "covars_proc_extra.rds"))

covars = read_rds(here::here("data", "processed", "covars_proc_extra.rds"))
# 
# cr = read_csv(here::here("data", "raw", "cr.csv"),
#               col_types = cols(ID = col_character())) %>% janitor::clean_names()
# 
# 
# covars_aki = 
#   covars %>% 
#   select(id, bin_aki, val_creatlst) %>% 
#   left_join(cr, by = "id") %>% 
#   mutate(increase = max_of_cr / val_creatlst,
#          change = max_of_cr - val_creatlst,
#          bin_aki7d = if_else(increase >= 1.5, 1, 0),
#          stage = case_when(bin_aki == 0 & bin_aki7d == 0 ~ 0,
#                            increase >= 3 | max_of_cr >= 4 & change >= 0.5 ~ 3,
#                            increase >= 2 & increase < 3 ~ 2,
#                            .default = 1))
# table(covars_aki$stage, covars_aki$bin_aki)
# sum(covars$bin_aki)
# 
# 
# ## calculate additional creatinine 