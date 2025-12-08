library(here)
# library(ggpubr)
library(tidyverse)
library(tidymodels)
library(patchwork)
library(gt) 
library(gtsummary)
library(paletteer)
library(survival)
library(survminer)
theme_set(theme_light())
force = FALSE
options(dplyr.summarise.inform=F)
col1 = "#F06719FF"; col2 = "#1BA3C6FF"
col1a = "#7873C0FF"; col2a = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
source(here::here("code", "utilities.R"))


##### ----- loading data --------
source(here::here("code", "prep_data.R"))


df = 
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars) %>% 
  left_join(map_only, by = "id") %>% 
  left_join(map_auc, by = "id") 


### --- proc survival data --- ####

los = read_csv(here::here("data", "raw", "los.csv"),
               col_types = cols(ID = col_character())) %>% janitor::clean_names()
icu = read_csv(here::here("data", "raw", "icu.csv"),
               col_types = cols(id = col_character())) %>% janitor::clean_names()

los %>% glimpse

icu %>% glimpse

df_surv = df %>% 
  left_join(los, by = "id") %>% 
  left_join(icu, by = "id") %>% 
  mutate(event = if_else(in_hospmor == "Y", 0, 1)) # event 1: disharged alive, 0 = not discharged (died)

nrow(df_surv)

df_surv_nd = 
  df_surv %>% 
  filter(in_hospmor == "N")
# summary(df_surv$icu_los)
# summary(df_surv$los)


result_uni = 
  map_dfr(.x = df_surv %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("Surv(los, event)", "~", x))
            model = coxph(formula, data = df_surv %>% mutate(across(contains("q"), ~.x / 5)))
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(1)
          })

result_uni %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 

result_uni_nd = 
  map_dfr(.x = df_surv_nd %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("Surv(los, event)", "~", x))
            model = coxph(formula, data = df_surv_nd %>% mutate(across(contains("q"), ~.x / 5)))
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(1)
          })

result_uni_nd %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 

los_model = coxph(Surv(los, event) ~ hypo_q1 + hypo_q2 + hypo_q3 + 
                    hypo_q4 + normo_q1 + normo_q2 + normo_q3 + normo_q4 +
                    cat_gender + val_creatlst + val_age + val_bmi +
                    bin_htn + bin_diabetes+ bin_stroke + bin_ef40 + 
                    bin_mi + bin_chf + bin_pvd + bin_cld + bin_betablocker +
                    bin_acearb + bin_statin + bin_redo + bin_emergent + 
                    val_crystalloid + val_cpbtime + cat_rbc + val_hematocrit,
                    data = df_surv %>% mutate(across(contains("q"), ~.x / 5)))

los_model %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 

los_model_nd = coxph(Surv(los, event) ~ hypo_q1 + hypo_q2 + hypo_q3 + 
                    hypo_q4 + normo_q1 + normo_q2 + normo_q3 + normo_q4 +
                    cat_gender + val_creatlst + val_age + val_bmi +
                    bin_htn + bin_diabetes+ bin_stroke + bin_ef40 + 
                    bin_mi + bin_chf + bin_pvd + bin_cld + bin_betablocker +
                    bin_acearb + bin_statin + bin_redo + bin_emergent + 
                    val_crystalloid + val_cpbtime + cat_rbc + val_hematocrit,
                  data = df_surv_nd %>% mutate(across(contains("q"), ~.x / 5)))

los_model_nd %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 


result_uni_icu = 
  map_dfr(.x = df_surv %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("Surv(icu_los, event)", "~", x))
            model = coxph(formula, data = df_surv %>% mutate(across(contains("q"), ~.x / 5)))
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(1)
          })

result_uni_icu %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 


result_uni_icu_nd = 
  map_dfr(.x = df_surv_nd %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("Surv(icu_los, event)", "~", x))
            model = coxph(formula, data = df_surv_nd %>% mutate(across(contains("q"), ~.x / 5)))
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(1)
          })

result_uni_icu_nd %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 




icu_model = coxph(Surv(icu_los, event) ~ hypo_q1 + hypo_q2 + hypo_q3 + 
                    hypo_q4 + normo_q1 + normo_q2 + normo_q3 + normo_q4 +
                    cat_gender + val_creatlst + val_age + val_bmi +
                    bin_htn + bin_diabetes+ bin_stroke + bin_ef40 + 
                    bin_mi + bin_chf + bin_pvd + bin_cld + bin_betablocker +
                    bin_acearb + bin_statin + bin_redo + bin_emergent + 
                    val_crystalloid + val_cpbtime + cat_rbc + val_hematocrit,
                  data = df_surv %>% mutate(across(contains("q"), ~.x / 5)))

icu_model %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 

icu_model_nd = coxph(Surv(icu_los, event) ~ hypo_q1 + hypo_q2 + hypo_q3 + 
                    hypo_q4 + normo_q1 + normo_q2 + normo_q3 + normo_q4 +
                    cat_gender + val_creatlst + val_age + val_bmi +
                    bin_htn + bin_diabetes+ bin_stroke + bin_ef40 + 
                    bin_mi + bin_chf + bin_pvd + bin_cld + bin_betablocker +
                    bin_acearb + bin_statin + bin_redo + bin_emergent + 
                    val_crystalloid + val_cpbtime + cat_rbc + val_hematocrit,
                  data = df_surv_nd %>% mutate(across(contains("q"), ~.x / 5)))

icu_model_nd  %>%
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(p.value = format.pval(p.value, digits = 3)) 

