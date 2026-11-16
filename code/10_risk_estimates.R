library(here)
# library(ggpubr)
# library(tidyverse)
library(dplyr)
library(tidymodels)
library(patchwork)
library(gt) 
library(gtsummary)
library(paletteer)
theme_set(theme_light())
force = FALSE
options(dplyr.summarise.inform=F)
col1 = "#F06719FF"; col2 = "#1BA3C6FF"
col1a = "#7873C0FF"; col2a = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
source(here::here("code", "utilities.R"))

covars = read_rds(here::here("data", "processed", "covars_proc.rds"))
map_ci = read_rds(here::here("data", "analysis", "mapci_ranges.rds"))


wake_covars = 
  covars %>% 
  filter(id %in% map_ci$id) %>% 
  rename(bin_aki48h = bin_aki) %>% 
  select(id, val_age, cat_gender, val_bmi, val_egfr,
         bin_htn, bin_diabetes, bin_mi,
         bin_chf, cat_ef,
         bin_iabp, bin_stroke, bin_pvd,
         bin_cld, bin_betablocker, bin_acearb,
         bin_statin, val_creatlst, bin_redo, bin_emergent,
         bin_aki48h, euro_predmort, val_proctime, val_crystalloid,
         cat_rbc, val_cpbtime, val_hematocrit) %>% 
  mutate(cohort = "WAKE") %>% 
  mutate(val_ef = 
           case_when(cat_ef == "20-25%" ~ 22.5,
                     cat_ef == "40-45%" ~ 42.5,
                     cat_ef == "45-50%" ~ 47.5,
                     cat_ef == "65-70%" ~ 67.5,
                     cat_ef == "15-20%" ~ 17.5,
                     cat_ef == "30-35%" ~ 32.5,
                     cat_ef == "35-40%" ~ 37.5,
                     cat_ef == "25-30%" ~ 27.5,
                     cat_ef == "70-75%" ~ 77.5,
                     cat_ef == "50-55%" ~ 52.5,
                     cat_ef == "60-65%" ~ 62.5,
                     cat_ef == "55-60%" ~ 57.5)) %>% 
  mutate(bin_ef40 = case_when(
    cat_ef %in% c("15-20%", "20-25%", "25-30%", "30-35%", "35-40%") ~ "<=40%",
    cat_ef %in% c("40-45%" ,"45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%") ~ ">40%",
    .default = NA_character_)) 

wake_covars =
  wake_covars %>% 
  drop_na()

map_ci = map_ci %>% 
  filter(id %in% wake_covars$id) 

hemo = read_rds(here::here("data", "processed", "hemo_analytic.rds"))  %>% 
  filter(!is.na(cat_cpb)) %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "intra", "post"))) %>% 
  filter(id %in% wake_covars$id) 


map_only = 
  hemo %>% 
  group_by(id) %>% 
  summarize(map_65 = sum(val_map < 65, na.rm = TRUE),
            .groups = "drop")

four_cats = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id) %>% 
  summarize(low_map_low_ci = sum(val_map < 65 & val_ci <= 2, na.rm = TRUE),
            low_map_hi_ci = sum(val_map < 65 & val_ci > 2, na.rm = TRUE),
            hi_map_low_ci = sum(val_map >= 65 & val_ci <= 2, na.rm = TRUE),
            hi_map_hi_ci = sum(val_map >= 65 & val_ci > 2, na.rm = TRUE),
            .groups = "drop")


four_cats_cpb = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id, cat_cpb) %>% 
  summarize(low_map_low_ci = sum(val_map < 65 & val_ci <= 2, na.rm = TRUE),
            low_map_hi_ci = sum(val_map < 65 & val_ci > 2, na.rm = TRUE),
            hi_map_low_ci = sum(val_map >= 65 & val_ci <= 2, na.rm = TRUE),
            hi_map_hi_ci = sum(val_map >= 65 & val_ci >2, na.rm = TRUE),
            .groups = "drop")

four_cats_cpb %>% 
  group_by(cat_cpb) %>% 
  summarize(med = median(low_map_low_ci),
            q1 = quantile(low_map_low_ci, 0.25),
            q3 =  quantile(low_map_low_ci, 0.75))

ci_ranges = read_rds(here::here("data", "analysis", "ci_ranges.rds")) %>% 
  filter(id %in% wake_covars$id)
ci_data = read_rds(here::here("data", "analysis", "ci_data.rds")) %>% 
  filter(id %in% wake_covars$id)


map_ci_post = read_rds(here::here("data", "analysis", "mapci_ranges_post.rds")) %>% 
  filter(id %in% wake_covars$id)

map_ci_pre = read_rds(here::here("data", "analysis", "mapci_ranges_pre.rds")) %>% 
  filter(id %in% wake_covars$id)


map_ci_cpb = 
  map_ci_post %>% mutate(cat_cpb = "post") %>% 
  bind_rows(map_ci_pre %>% mutate(cat_cpb = "pre")) 

ci_ranges4 = read_rds(here::here("data", "analysis", "ci_ranges4.rds")) %>% 
  filter(id %in% wake_covars$id)

ci_ranges7 = read_rds(here::here("data", "analysis", "ci_ranges7.rds")) %>% 
  filter(id %in% wake_covars$id)


map_data = read_rds(here::here("data", "analysis", "map_ranges.rds")) %>% 
  filter(id %in% wake_covars$id)

mapci_tert = read_rds(here::here("data", "analysis", "mapci_ranges_tertile.rds")) %>% 
  filter(id %in% wake_covars$id)

mapci_quint = read_rds(here::here("data", "analysis", "mapci_ranges_quintile.rds")) %>% 
  filter(id %in% wake_covars$id)

df = 
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars) 
df = 
  df %>% 
  left_join(map_only, by = "id")

df2 = 
  df %>% 
  left_join(four_cats, by = "id") %>% 
  select(-map_65) %>% 
  mutate(across(contains("map"), ~.x / 5))


# original model 
m1 = glm(
  bin_aki48h ~ low_map_low_ci+  low_map_hi_ci + 
    hi_map_low_ci + hi_map_hi_ci + 
    cat_gender+
    val_creatlst+
    val_age+
    val_bmi+
    bin_htn+
    bin_diabetes+
    bin_stroke+
    bin_ef40+
    bin_mi+
    bin_chf+
    bin_pvd+
    bin_cld+
    bin_betablocker+
    bin_acearb+
    bin_statin+
    bin_redo+
    bin_emergent + 
    val_crystalloid + 
    val_cpbtime + 
    cat_rbc+
    val_hematocrit,
  data = df2,
  family = binomial()
)

### how many pts hvae <5 min in hypotension/CI >2

df2 %>% 
  filter(low_map_hi_ci < 1) %>% 
  nrow()

### how many pts hvae <5 min in normotension/CI<=2

df2 %>% 
  filter(hi_map_low_ci < 1) %>% 
  nrow()



### increase MAP 
### Sub 5 min MAP < 65, CI < 2 to 5 minutes MAP >= 65, CI < 2
aki_probs0 = predict(m1, newdata = df2, type = "response")

data_shifted = 
  df2 %>% 
  mutate(low_map_low_ci = low_map_low_ci + 1,
         hi_map_low_ci = hi_map_low_ci - 1)

aki_probs_shifted0 = predict(m1, newdata = data_shifted, type = "response")
# 
# mean(aki_probs) - mean(aki_probs_shifted)
# mean(aki_probs) / mean(aki_probs_shifted)
# mean(aki_probs)
# mean(aki_probs_shifted)
boot_res = map_dfr(.x = 1:1000,
                   .f = function(s){
                     set.seed(s)
                     sampled_data = df2 %>% slice_sample(n = nrow(df2), replace = TRUE)
                     model = update(m1, data = sampled_data)
                     aki_probs = predict(model, newdata = sampled_data, type = "response")
                     data_shifted = 
                       sampled_data %>% 
                       mutate(low_map_low_ci = low_map_low_ci + 1,
                              hi_map_low_ci = hi_map_low_ci - 1) 
                     
                     aki_probs_shifted = predict(model, newdata = data_shifted, type = "response")
                     
                     tibble(risk_diff =  mean(aki_probs_shifted) - mean(aki_probs),
                            risk_ratio = mean(aki_probs_shifted) / mean(aki_probs))
                   })



res = tibble(
  risk_orig = mean(aki_probs0),
  risk_shift = mean(aki_probs_shifted0),
  risk_diff = risk_shift - risk_orig,
  lower_bound_diff = quantile(boot_res$risk_diff, .025),
  upper_bound_diff = quantile(boot_res$risk_diff, .975),
  risk_ratio = risk_shift / risk_orig,
  lower_bound_ratio = quantile(boot_res$risk_ratio, .025),
  upper_bound_ratio = quantile(boot_res$risk_ratio, .975))

res_1 =
  res %>% mutate(across(everything(), ~round(.x, 3))) %>%
  mutate(risk_diff = glue::glue("{risk_diff} ({lower_bound_diff}, {upper_bound_diff})"),
         risk_ratio = glue::glue("{risk_ratio} ({lower_bound_ratio}, {upper_bound_ratio})")) %>% 
  select(risk = risk_orig, 
         risk_shifted = risk_shift, 
         risk_diff, 
         risk_ratio) %>% 
  mutate(shift = "Increase MAP")

res_1 

### average person

mode_cat = function(x) names(which.max(table(x)))
average_person =
  df2 %>% 
  summarize(across(c(val_creatlst, val_age, val_bmi, val_crystalloid, val_cpbtime, val_hematocrit, contains("map")), mean),
            across(c(contains("cat"), contains("bin")), mode_cat)) %>% 
  mutate(across(c(contains("cat")), as.factor),
         across(c(bin_redo, bin_acearb, bin_betablocker, bin_emergent, bin_iabp, bin_chf, bin_htn,
                  bin_diabetes, bin_mi, bin_stroke, bin_pvd, bin_cld, bin_statin, bin_acearb),
                as.numeric))

aki_prob = predict(m1, newdata = average_person, type = "response"); aki_prob

average_person_mod = 
  average_person %>% 
  mutate(low_map_low_ci = low_map_low_ci + 1,
         hi_map_low_ci = hi_map_low_ci - 1)

aki_prob_mod = predict(m1, newdata = average_person_mod, type = "response"); aki_prob_mod

aki_prob_mod - aki_prob 
aki_prob_mod / aki_prob 

### Sub 5 min MAP < 65, CI < 2 to 5 minutes MAP < 65, CI >= 2
##### increase CI 

data_shifted = 
  df2 %>% 
  mutate(low_map_low_ci = low_map_low_ci + 1,
         low_map_hi_ci = low_map_hi_ci - 1)

aki_probs_shifted_ci = predict(m1, newdata = data_shifted, type = "response")
# 
# mean(aki_probs) - mean(aki_probs_shifted)
# mean(aki_probs) / mean(aki_probs_shifted)
# mean(aki_probs)
# mean(aki_probs_shifted)
boot_res = map_dfr(.x = 1:1000,
                   .f = function(s){
                     set.seed(s)
                     sampled_data = df2 %>% slice_sample(n = nrow(df2), replace = TRUE)
                     model = update(m1, data = sampled_data)
                     aki_probs = predict(model, newdata = sampled_data, type = "response")
                     data_shifted = 
                       sampled_data %>% 
                       mutate(low_map_low_ci = low_map_low_ci + 1,
                              low_map_hi_ci = low_map_hi_ci - 1)
                     
                     
                     aki_probs_shifted = predict(model, newdata = data_shifted, type = "response")
                     
                     tibble(risk_diff =  mean(aki_probs_shifted) - mean(aki_probs),
                            risk_ratio = mean(aki_probs_shifted) / mean(aki_probs))
                   })



res = tibble(
  risk_orig = mean(aki_probs0),
  risk_shift = mean(aki_probs_shifted_ci),
  risk_diff = risk_shift - risk_orig,
  lower_bound_diff = quantile(boot_res$risk_diff, .025),
  upper_bound_diff = quantile(boot_res$risk_diff, .975),
  risk_ratio = risk_shift / risk_orig,
  lower_bound_ratio = quantile(boot_res$risk_ratio, .025),
  upper_bound_ratio = quantile(boot_res$risk_ratio, .975))

res_2 =
  res %>% mutate(across(everything(), ~round(.x, 3))) %>%
  mutate(risk_diff = glue::glue("{risk_diff} ({lower_bound_diff}, {upper_bound_diff})"),
         risk_ratio = glue::glue("{risk_ratio} ({lower_bound_ratio}, {upper_bound_ratio})")) %>% 
  select(risk = risk_orig, 
         risk_shifted = risk_shift, 
         risk_diff, 
         risk_ratio) %>% 
  mutate(shift = "Increase CI")
res_2 



res_2 %>% 
  bind_rows(res_1)

### average person
average_person_mod = 
  average_person %>% 
  mutate(low_map_low_ci = low_map_low_ci + 1,
         low_map_hi_ci = low_map_hi_ci - 1)

aki_prob_mod = predict(m1, newdata = average_person_mod, type = "response"); aki_prob_mod

aki_prob_mod - aki_prob 
aki_prob_mod / aki_prob 


average_person %>% t()

average_person %>% mutate(across(contains("map"), ~.x * 5)) %>% select(contains("map"))
