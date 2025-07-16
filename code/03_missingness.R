check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Example usage:
required_packages <- c("readr", "dplyr", "tidyverse",
                       "lubridate", "purrr", "stringr", "here",
                       "slider", "furrr", "future", "patchwork", "zoo")
check_and_install_packages(required_packages)
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

hemo_interp = readRDS(here::here("data", "processed", "hemo_timeseries_interp.rds"))


covars = readxl::read_excel(
  here::here("data", "raw", "WAKE.flatfile.IDsfixed.11.27.24.xlsx")) %>% 
  mutate(ID = as.character(ID)) %>% 
  janitor::clean_names()

mult_proc = 
  covars %>% 
  # filter(cabg %in% c("CABG ON PUMP", "CABG ON  PUM")) %>% 
  group_by(id) %>% 
  mutate(n_proc = n()) %>% 
  filter(n_proc > 1) %>% 
  ungroup() %>% 
  pull(id)

ecmo = covars %>% 
  filter(!is.na(ecmo)) %>% 
  pull(id) 

crt4 = 
  covars %>% 
  filter(base_line_cr > 4) %>% 
  pull(id) 

# include_ids = include$id[include$id %in% unique(hemo_interp$id)]

# hemo_interp = 
#   hemo_interp %>% 
#   filter(id %in% include_ids)

hemo_interp %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  mutate(start = min(time[!is.na(val_map)]),
         end = max(time[!is.na(val_map)])) %>% 
  filter(between(time, start, end)) %>% 
  group_by(id, cat_cpb) %>% 
  summarize(n = n(),
            map_values = sum(is.na(val_map)),
            ci_values = sum(is.na(val_ci)),
            cvp_values  = sum(is.na(val_cvp))) %>% 
  mutate(minutes = n  / 12,
         map_miss = map_values / 12,
         ci_miss = ci_values / 12,
         cvp_miss = cvp_values / 12,
         pct_map = (map_miss / minutes) * 100,
         pct_ci = (ci_miss / minutes) * 100,
         pct_cvp = (cvp_miss / minutes) * 100) %>% 
  group_by(cat_cpb) %>% 
  summarize(across(minutes:pct_cvp, list(mean = ~mean(.x), sd = ~sd(.x))))



### missingness 

missing_map_all =
  hemo_interp %>% 
  group_by(id) %>% 
  summarize(non_na_map = sum(!is.na(val_map))) %>% 
  filter(non_na_map == 0) %>% 
  pull(id)

missing_map_all

missing_cvp_all = 
  hemo_interp %>% 
  group_by(id) %>% 
  summarize(non_na_cvp = sum(!is.na(val_cvp))) %>% 
  filter(non_na_cvp == 0) %>% 
  pull(id)

missing_cvp_all


missing_ci_all = 
  hemo_interp %>% 
  group_by(id) %>% 
  summarize(non_na_ci = sum(!is.na(val_ci))) %>% 
  filter(non_na_ci == 0) %>% 
  pull(id)

missing_ci_all
# none 

# missing map or ci during pre or post bypass 
get_missing_period = function(variable, period, hemodf){
  hemodf %>% 
    rename(var = all_of(variable)) %>% 
    filter(cat_cpb == period) %>% 
    group_by(id) %>% 
    summarize(non_na = sum(!is.na(var))) %>% 
    filter(non_na == 0) %>% 
    pull(id) 
}

missing_map_precpb <- get_missing_period("val_map", "pre", hemo_interp); missing_map_precpb
missing_map_postcpb <- get_missing_period("val_map", "post", hemo_interp); missing_map_postcpb

missing_ci_precpb <- get_missing_period("val_ci", "pre", hemo_interp); missing_ci_precpb
missing_ci_postcpb <- get_missing_period("val_ci", "post", hemo_interp); missing_ci_postcpb

missing_cvp_precpb <- get_missing_period("val_cvp", "pre", hemo_interp); missing_cvp_precpb
missing_cvp_postcpb <- get_missing_period("val_cvp", "post", hemo_interp); missing_cvp_postcpb


# last one: both >50% missing and <20 minutes of map or ci in one period
# denominator is minutes between first and last arterial BP measurement 
# get df that's only between first and last arterial BP measurement 


# summarize the percent missing and number of minutes present for each id and period
missing_summary = 
  hemo_interp %>%
  filter(cat_cpb!= "intra") %>% 
  group_by(id, cat_cpb) %>%
  summarize(
    elapsed = as.numeric(difftime(max(time), min(time), units = "mins"))+1,
    na_map = sum(is.na(val_map))/12,
    na_ci = sum(is.na(val_ci))/12,
    na_cvp = sum(is.na(val_cvp))/12,
    present_map = sum(!is.na(val_map))/12,
    present_ci = sum(!is.na(val_ci))/12,
    present_cvp = sum(!is.na(val_cvp))/12
  ) %>%
  ungroup() %>% 
  mutate(across(starts_with("na"), ~.x/elapsed*100, .names = "{.col}_pct")) %>% 
  rowwise() %>% 
  mutate(missing_map_50_20 = ifelse(na_map_pct > 50 & present_map < 20, 1, 0),
         missing_ci_50_20 = ifelse(na_ci_pct > 50 & present_ci < 20, 1, 0),
         missing_cvp_50_20 = ifelse(na_cvp_pct > 50 & present_cvp < 20, 1, 0))


missing_map_50_20_precpb <- 
  missing_summary %>% 
  filter(cat_cpb == "pre" & missing_map_50_20 == 1) %>% 
  pull(id)

missing_map_50_20_postcpb <- 
  missing_summary %>% 
  filter(cat_cpb == "post" & missing_map_50_20 == 1) %>% 
  pull(id)

missing_cvp_50_20_precpb <- 
  missing_summary %>% 
  filter(cat_cpb == "pre" & missing_cvp_50_20 == 1) %>% 
  pull(id)

missing_cvp_50_20_postcpb <- 
  missing_summary %>% 
  filter(cat_cpb == "post" & missing_cvp_50_20 == 1) %>% 
  pull(id)

missing_ci_50_20_precpb <- 
  missing_summary %>% 
  filter(cat_cpb == "pre" & missing_ci_50_20 == 1) %>% 
  pull(id)


missing_ci_50_20_postcpb <- 
  missing_summary %>% 
  filter(cat_cpb == "post" & missing_ci_50_20 == 1) %>% 
  pull(id)



## TO DO - MULTIPLE AXC? 
exclusion_summary =
  tibble(
    id = unique(hemo_interp$id)
  ) %>% 
  mutate(
    # has_multiple_axc = id %in% mult_axc_id, 
    is_missing_map_all = id %in% missing_map_all,
    is_missing_cvp_all = id %in% missing_cvp_all,
    is_missing_ci_all = id %in% missing_ci_all,
    is_missing_map_precpb = id %in% missing_map_precpb,
    is_missing_map_postcpb = id %in% missing_map_postcpb,
    is_missing_cvp_precpb = id %in% missing_cvp_precpb,
    is_missing_cvp_postcpb = id %in% missing_cvp_postcpb,
    is_missing_cvp_precpb = id %in% missing_ci_precpb,
    is_missing_cvp_postcpb = id %in% missing_ci_postcpb,
    is_missing_map_50_20_precpb = id %in% missing_map_50_20_precpb,
    is_missing_map_50_20_postcpb = id %in% missing_map_50_20_postcpb,
    is_missing_cvp_50_20_precpb = id %in% missing_cvp_50_20_precpb,
    is_missing_cvp_50_20_postcpb = id %in% missing_cvp_50_20_postcpb,
    is_missing_ci_50_20_precpb = id %in% missing_ci_50_20_precpb,
    is_missing_ci_50_20_postcpb = id %in% missing_ci_50_20_postcpb
    # manual_reivew = id %in% manual_review_ids
  ) 

exclusion_summary = exclusion_summary %>% 
  mutate(wrong_procedure = !(id %in% covars$id),
         excluded_creat = id %in% crt4,
         mult_proc = id %in% mult_proc) %>% 
  rowwise() %>% 
  mutate(exclude_hemo = if_any(contains("is_"), ~ .x == TRUE),
           excluded = case_when(if_any(contains("is_"), ~ .x == TRUE) ~ TRUE,
                              if_any(contains("has"), ~ .x == TRUE) ~ TRUE,
                              .default = FALSE)) %>% 
  ungroup() 

sum(exclusion_summary$excluded)

exclusion_summary %>% 
  filter(!wrong_procedure) %>% 
  filter(excluded_creat) %>% 
  nrow()

exclusion_summary %>% 
  filter(!wrong_procedure) %>% 
  filter(!excluded_creat) %>% 
  filter(mult_proc) %>% 
  nrow()

exclusion_summary %>% 
  filter(!wrong_procedure) %>% 
  filter(!excluded_creat) %>% 
  filter(!mult_proc) %>% 
  filter(exclude_hemo) %>% 
  nrow()


exclusion_summary %>% 
  filter(!wrong_procedure) %>% 
  filter(!excluded_creat) %>% 
  filter(!mult_proc) %>% 
  filter(!exclude_hemo) %>% 
  nrow()

write_rds(exclusion_summary, here::here("data", "processed", "exclusion_summary.rds"))
sum(exclusion_summary$excluded)/nrow(exclusion_summary)

ex_ids = exclusion_summary %>% 
  filter(excluded) %>% 
  pull(id)

hemo_analytic = 
  hemo_interp %>% 
  filter(!(id %in% ex_ids))

write_rds(hemo_analytic, here::here("data", "processed", "hemo_analytic.rds"),
          compress = "xz")



