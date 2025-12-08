library(here)
library(ggpubr)
library(tidyverse)
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


## Substitution 

# Sub 5 minutes MAP < 65, CI in quartile 1 (<2) to 5 minutes MAP < 65, CI in quartile 2 (2-2.4)


m1 = 
  glm(
    bin_aki48h ~ .,
    data = df %>% select(
      contains("q"),
      bin_aki48h,
      cat_gender,
      val_creatlst,
      val_age,
      val_bmi,
      bin_htn,
      bin_diabetes,
      bin_stroke,
      bin_ef40,
      bin_mi,
      bin_chf,
      bin_pvd,
      bin_cld,
      bin_betablocker,
      bin_acearb,
      bin_statin,
      bin_redo,
      bin_emergent,
      val_crystalloid,
      val_cpbtime,
      cat_rbc,
      val_hematocrit
    ) %>% mutate(across(c(contains("q"), contains("map")), ~ .x / 5)),
    family = binomial()
  )

k  = matrix(c(0, -1, 1, rep(0, length(coef(m1)) - 3)), 1)
contr =  multcomp::glht(m1, linfct = k)
summary(contr)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 1 == 0 -0.17617    0.08846  -1.992   0.0464 *
## 95 % CI 
pt_est = -0.17617
se = 0.08846
exp(pt_est)
# [1] 0.8384754
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1] 0.7050041 0.9972155


# Sub 5 minutes MAP < 65, CI in quartile 1 (<2) to 5 minutes MAP >= 65, CI in quartile 1 

k  = matrix(c(0, -1, 0, 0, 0, 1 , rep(0, length(coef(m1)) - 6)), 1)
contr =  multcomp::glht(m1, linfct = k)
summary(contr)



# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 1 == 0 -0.13051    0.05307  -2.459   0.0139 *
  
## 95 % CI 
pt_est = -0.13051
se = 0.05307
exp(pt_est)
# [1] 0.8776477
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1] 0.7909447 0.9738551



# Sub 5 minutes MAP < 65, CI <= 2 to 5 minutes MAP < 65, CI > 2
df2 = 
  df %>% 
  left_join(four_cats, by = "id") %>% 
  select(-map_65)
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
  data = df2 %>%  mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
k  <- matrix(c(0, -1, 1, rep(0, length(coef(m1)) - 3)), 1)
contr <- multcomp::glht(m1, linfct = k)
summary(contr)

# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 1 == 0 -0.13616    0.05819   -2.34   0.0193 *
  
## 95 % CI 
pt_est =-0.13616  
se = 0.05819 
exp(pt_est)
# [1] 0.872703
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1] 0.7786354 0.9781350


# Sub 5 minutes MAP < 65, CI <= 2 with 5 minutes MAP >= 65, CI <= 2


k  <- matrix(c(0, -1, 0, 1, rep(0, length(coef(m1)) - 4)), 1)
contr <- multcomp::glht(m1, linfct = k)
summary(contr)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 1 == 0 -0.11773    0.05052   -2.33   0.0198 *
pt_est = -0.11773
se = 0.05052
exp(pt_est)
# [1] 0.888936
exp(c(pt_est - (1.96 * se), pt_est + (1.96 * se)))
# [1]  0.8051319 0.9814632


df2


library(codaredistlm)

# Use the 'fat_data' dataset included in the package as an example
# y = "fat" (the outcome)
# comps = c("sl", "sb", "lpa", "mvpa") (sleep, sedentary, light, mod-vig activity)
# covars = c("sibs", "parents", "ed") (confounders)
# deltas = time to reallocate (in units of a day, so 60 min = 60 / (24*60))

# 1. Run the analysis
# This function does it all: fits the model and calculates substitutions
# 'comparisons = "one-v-one"' calculates all pairwise substitutions
sub_analysis <- predict_delta_comps(
  dataf = fat_data,
  y = "fat",
  comps = c("sl", "sb", "lpa", "mvpa"),
  covars = c("sibs", "parents", "ed"),
  deltas = seq(-60, 60, by = 10) / (24 * 60), # Check -60 to +60 min in 10-min steps
  comparisons = "one-v-one", # Or "prop-realloc" for one-vs-all
  alpha = 0.05 # For 95% CIs
)

# 2. Plot the results
# This will generate a plot showing the predicted change in the "fat" outcome
# when reallocating time between any two behaviors.
plot_delta_comp(sub_analysis)

sub_analysis = 
  predict_delta_comps(
    y = "bin_aki48h",
    comps = c("low_map_low_ci", "low_map_hi_ci",
              "hi_map_low_ci", "hi_map_hi_ci"),
    covars = c("val_age", "val_creatlst"),
    comparisons = "one-v-one",
    deltas = seq(0, 1, .2),
    dataf = df2 %>% 
      mutate(across(contains("map"), ~.x * 5 + 1))
  )

plot_delta_comp(sub_analysis)
