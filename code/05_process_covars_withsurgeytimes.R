# install.packages("fuzzyjoin")
# library(fuzzyjoin)
rm(list = ls()) 
# TAKE INDIVIDUAL FILES AND CREATE ONE DF WITH ALL TIME SERIES DATA 
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
                       "slider", "furrr", "future", "patchwork", "zoo", "fuzzyjoin")
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
#library(fuzzyjoin)
options(digits.secs = 2)
force = TRUE
source(here::here("code", "utilities.R"))
ncores = parallelly::availableCores() - 1
options(dplyr.summarise.inform=FALSE)

hemo_data = 
  read_rds(here::here("data", "processed", "hemo_analytic.rds"))

cpb_time = 
  hemo_data %>% 
  group_by(id) %>% 
  summarize(start = max(time[cat_cpb == "pre" & !is.na(val_ci)]),
            end = min(time[cat_cpb == "post" & !is.na(val_ci)]),
            .groups = "drop") %>% 
  mutate(val_cpbtime = difftime(end, start, unit = "mins") %>% as.numeric)

# Read cleaned anesthesia and surgery event times
or_event_times <- readr::read_rds(
  here::here(
    "data",
    "processed",
    "covariates_with_or_times.rds"
  )
) %>%
  transmute(
    id = as.character(id),
    anesthesia_ready = as.POSIXct(anes_ready, tz = "UTC"),
    anesthesia_start = as.POSIXct(anes_start, tz = "UTC"),
    surgery_start = as.POSIXct(surgery_start, tz = "UTC"),
    surgery_start_source = surgery_start_source,
    cpb_start = as.POSIXct(cpb_start, tz = "UTC"),
    cpb_end = as.POSIXct(cpb_end, tz = "UTC"),
    anesthesia_end = as.POSIXct(anes_end, tz = "UTC")
  ) %>%
  distinct()

duplicate_or_event_ids <- or_event_times %>%
  count(id, name = "n") %>%
  filter(n > 1L)

if (nrow(duplicate_or_event_ids) > 0L) {
  print(duplicate_or_event_ids)
  
  stop(
    paste0(
      "Some IDs have multiple OR-event records. ",
      "Do not join by ID alone until the correct operation is identified."
    )
  )
}

covars = readxl::read_excel(
  here::here("data", "raw", "WAKE.flatfile.IDsfixed.11.27.24.xlsx")) %>% 
  mutate(ID = as.character(ID)) %>% 
  janitor::clean_names()


covars_lee = readxl::read_xls(here::here("data", "raw", "duplicated_covars_Goeddel_ver2.xls")) %>% 
  mutate(id = as.character(id)) %>% 
  janitor::clean_names()

covars_keep = 
  covars_lee %>% 
  filter(is.na(exclude)) %>% 
  mutate(across(c(iabp, impella, no_sternotomy, redo,
                  pat_height, pat_weight), as.numeric),
         across(contains("date"), lubridate::as_datetime))

covars = 
  covars %>% 
  filter(!(id %in% covars_lee$id)) %>% 
  bind_rows(covars_keep %>% select(-contains("exclude")))


# duplicated_cov = 
#   covars %>% group_by(id) %>% mutate(n = n()) %>% ungroup() %>% filter(n>1) %>% arrange(id) %>% distinct()
# write_csv(duplicated_cov, here::here("data_qc", "duplicated_covars.csv"))

comorbidities = read_csv(
  here::here("data",
             "raw", "table2.comorbidities.wake.2.12.25.csv"),
  col_types = cols(ID = col_character())
) %>% janitor::clean_names() %>% 
  select(-x2) %>% 
  group_by(id) %>% 
  slice(1) %>% 
  ungroup() 




# duplicated_comob = 
#   comorbidities %>% group_by(id) %>% mutate(n = n()) %>% ungroup() %>% filter(n>1) %>% arrange(id) %>% distinct()
# write_csv(duplicated_comob, here::here("data_qc", "duplicated_comorbidities.csv"))

comorbidities_proc = 
  comorbidities %>% 
  mutate(across(hypertension:acei_arb, ~if_else(.x == "Y", 1, 0))) %>% 
  rename(bin_htn = hypertension,
         bin_diabetes = diabates,
         bin_mi = myocardial_infarction,
         bin_chf = congestive_heart_failure,
         bin_stroke = stroke,
         bin_pvd = peripheral_vascular_disease,
         bin_cld = chronic_lung_disease,
         bin_emergent = emergency,
         bin_betablocker = beta_blocker,
         bin_statin = statin,
         bin_acearb = acei_arb,
         val_crystalloid = crystalliods_in_ml,
         cat_ef = ef) %>% 
  mutate(ef_30_50 = if_else(cat_ef %in% c("30-35%", "35-40%", "40-45%", "45-50%"), 1, 0),
         ef_l30 = if_else(cat_ef %in% c("10-15%", "15-20%", "20-25%", "25-30%"), 1, 0),
         cat_rbc = factor(case_when(rbc == 0 ~ "0",
                                    rbc > 350 ~ ">1", 
                                    rbc <= 350 ~ "1",
                                    .default = NA_character_)))


covars_proc <- covars %>%
  mutate(
    val_proctime = as.numeric(
      difftime(
        an_stop_datetime,
        an_start_datetime,
        units = "mins"
      )
    )
  ) %>%
  select(
    id,
    val_age = patient_age,
    val_bmi = bmi,
    cat_ethnicity = patient_ethnicity,
    val_creatlst = base_line_cr,
    cat_aki = aki,
    bin_iabp = iabp,
    gender = gender_1_female_2_male,
    val_proctime,
    surgery_date,
    redo
  ) %>%
  mutate(
    id = as.character(id),
    
    cat_gender = if_else(
      gender == 1,
      "Female",
      "Male"
    ),
    
    bin_aki = case_when(
      cat_aki == "No"  ~ 0,
      cat_aki == "Aki" ~ 1,
      TRUE             ~ NA_real_
    ),
    
    bin_redo = case_when(
      redo == 1 ~ 1,
      TRUE      ~ 0
    ),
    
    bin_iabp = case_when(
      is.na(bin_iabp) ~ 0,
      TRUE            ~ 1
    )
  ) %>%
  left_join(
    or_event_times,
    by = "id"
  ) %>%
  mutate(
    # Time from anesthesia start to incision/procedure start
    val_anes_to_surgery_min = as.numeric(
      difftime(
        surgery_start,
        anesthesia_start,
        units = "mins"
      )
    ),
    
    # Recalculate anesthesia duration from the cleaned timestamps
    val_anesthesia_duration_min = as.numeric(
      difftime(
        anesthesia_end,
        anesthesia_start,
        units = "mins"
      )
    ),
    
    val_cpb_event_duration_min = as.numeric(
      difftime(
        cpb_end,
        cpb_start,
        units = "mins"
      )
    )
  )


covars_all = 
  covars_proc %>% 
  full_join(comorbidities_proc, by = "id")

covars_all = 
  covars_all %>% 
  mutate(cr_g2_26 = if_else(val_creatlst > 2.26, 1, 0))

final = calc_euroscore(covars_all)

final = final %>% 
  left_join(cpb_time %>% select(id, val_cpbtime), by = "id") 



## calculate egfr
final = final %>% 
  mutate(
    val_egfr = case_when(
      cat_gender == "Female" & val_creatlst <= 0.7 ~ (142 * (val_creatlst/0.7)^(-0.241) * (0.9938^val_age) * 1.012),
      cat_gender == "Female" & val_creatlst > 0.7 ~ (142 * (val_creatlst/0.7)^(-1.2) * (0.9938^val_age) * 1.012),
      cat_gender == "Male" & val_creatlst <= 0.9 ~ (142 * (val_creatlst/0.9)^(-0.302) * (0.9938^val_age)),
      cat_gender == "Male" & val_creatlst > 0.9 ~ (142 * (val_creatlst/0.9)^(-1.2) * (0.9938^val_age)),
      .default = NA_real_
    )
  )


### add in extra covariates from final wake extract 
extra_covars = readxl::read_excel(here::here("data", "raw", "jhu.9.18.25.preophemoglobin.xlsx"),
                                  col_types = c("text", "text", "numeric", "numeric")) %>% 
  janitor::clean_names() %>% 
  mutate(anemia_before_surgery = if_else(anemia_before_surgery == "N", 0, 1)) %>% 
  rename(val_hematocrit = hematocrit, 
         val_hemoglobin = hemoglobin,
         bin_preopanemia = anemia_before_surgery)

final = 
  final %>% 
  left_join(extra_covars, by = "id")

write_rds(final, here::here("data", "processed", "covars_proc_with_surgerytimes.rds"))



### extra 
ex_euro = read_csv(here::here("data", "ex_euroscore.csv"))
ex_euro = 
  ex_euro %>% janitor::clean_names()
summary(ex_euro$x3)
hist(ex_euro$x3)
summary(final$euro_predmort)
hist(final$euro_predmort * 100)
nrow(final)

sample = final %>% 
  slice_sample(n = 10)
g_20 = 
  final %>% 
  filter(euro_predmort > 20)

write_csv(g_20, here::here("data", "euroscore_above_20.csv"))
write_csv(sample, here::here("data", "euroscore_sample10pts.csv"))