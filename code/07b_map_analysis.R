# map analysis 
library(here)
# library(ggpubr)
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


df = df %>% 
  mutate(map_auc = map_auc / 10)
#### --- map regressions ---- ### 

df %>% colnames

df %>% 
  ggplot(aes(x= factor(bin_aki48h), y = map_65)) + 
  geom_boxplot() +
  labs(x = "AKI", y = "Minutes with MAP < 65") +
  scale_x_discrete(labels = c("No", "Yes"))

df %>% 
  ggplot(aes(x= factor(bin_aki48h), y = map_60)) + 
  geom_boxplot() +
  labs(x = "AKI", y = "Minutes with MAP < 60") +
  scale_x_discrete(labels = c("No", "Yes"))

df %>% 
  ggplot(aes(x= factor(bin_aki48h), y = map_55)) + 
  geom_boxplot() +
  labs(x = "AKI", y = "Minutes with MAP < 55") +
  scale_x_discrete(labels = c("No", "Yes"))

df %>% 
  ggplot(aes(x= factor(bin_aki48h), y = map_auc)) + 
  geom_boxplot() +
  labs(x = "AKI", y = "MAP AUC 65") +
  scale_x_discrete(labels = c("No", "Yes"))

# uniariate regressions 
result_uni = 
  map_dfr(.x = df %>% select(contains("map")) %>% colnames,
          .f = function(x){
            formula = as.formula(paste("bin_aki48h", "~", x))
            model = glm(formula, data = df %>% mutate(across(c(contains("map") & !contains("auc")), ~.x / 5)), family = binomial)
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(2)
          })
result_uni

## fully adjusted model with MAP < 65
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
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
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map65.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map65.csv"))
} 


## fully adjusted model with MAP < 60
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_60,
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
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)

summary(df$map_60)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   14.00   22.00   26.55   35.00  179.00 
if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map60.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map60.csv"))
} 

# fully adjusted model with MAP 55
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_55,
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
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
# 
# summary(df$map_55)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    5.00    9.00   11.65   15.00   95.00 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map55.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map55.csv"))
} 

# fully adjusted model with MAP AUC

# 
# summary(df$map_auc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0   204.8   319.0   378.0   472.2  2791.0 


model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_auc,
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
  ),
  family = binomial()
)
# 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_mapauc.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_mapauc.csv"))
} 


### ---- remove CPB time ---- ### 

## fully adjusted model with MAP < 65
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_65,
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
    # val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map65_nocpb.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map65_nocpb.csv"))
} 


## fully adjusted model with MAP < 60
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_60,
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
    # val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)


# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   14.00   22.00   26.55   35.00  179.00 
if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map60_nocpb.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map60_nocpb.csv"))
} 

# fully adjusted model with MAP 55
model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_55,
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
    # val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("map"), ~ .x / 5)),
  family = binomial()
)
# 
# summary(df$map_55)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    5.00    9.00   11.65   15.00   95.00 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_map55_nocpb.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_map55_nocpb.csv"))
} 

# fully adjusted model with MAP AUC

# 
# summary(df$map_auc)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.0   204.8   319.0   378.0   472.2  2791.0 


model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_auc,
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
    # val_cpbtime,
    cat_rbc,
    val_hematocrit
  ),
  family = binomial()
)
# 

model = glm(
  bin_aki48h ~ .,
  data = df %>% select(
    map_auc,
    bin_aki48h,
    val_creatlst
  ),
  family = binomial()
)
# 

if(!file.exists(here::here("manuscript", "revision", "adjusted_regression_mapauc_nocpb.csv")) ||
   force) {
  model %>%
    broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(p.value = format.pval(p.value, digits = 3)) %>%
    write_csv(., here::here("manuscript", "revision", "adjusted_regression_mapauc_nocpb.csv"))
} 
