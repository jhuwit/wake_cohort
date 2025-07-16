# make filenames df 

library(tidyverse)


covars = readxl::read_excel(
  here::here("data", "raw", "WAKE.flatfile.IDsfixed.11.27.24.xlsx")) %>% 
  mutate(ID = as.character(ID)) %>% 
  janitor::clean_names()


or_times = read_csv(here::here("data", "raw", "WAKE_OR_Times.IDfixed.11.29.24.csv")) %>% 
  janitor::clean_names()

individual_files = list.files(here::here("data", "raw", "individual"),
                              recursive = TRUE,
                              full.names = TRUE)


ts_ids = sub(".*individual\\/(.+).csv.*", "\\1", individual_files)
covar_ids = unique(covars$id)
or_ids = unique(or_times$id)

all_ids = c(ts_ids, covar_ids, or_ids) %>% unique() 

mutate(across(contains("file"), ~file.exists(.x), .names = "{.col}_exists"))


filenames_df = 
  tibble(id = all_ids) %>% 
  mutate(individualfile = here::here("data", "raw", "individual", paste0(id, ".csv")),
         processedfile = here::here("data", "processed", "interp_individual", paste0(id, ".rds")),
         covars_exists = id %in% covar_ids,
         ortimes_exists = id %in% or_ids) %>% 
  mutate(across(contains("file"), ~file.exists(.x), .names = "{.col}_exists"))

           
