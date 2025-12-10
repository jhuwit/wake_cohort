library(here)
library(ggpubr)
library(tidyverse)
# library(tidymodels)
library(patchwork)
library(gt) 
library(gtsummary)
library(paletteer)
theme_set(theme_light())
force = FALSE 
options(dplyr.summarise.inform=F)
col1 = "#F06719FF"; col2 = "#1BA3C6FF"
# "#009E73FF"
# "#D55E00FF"

col1a = "#7873C0FF"; col2a = "#21B087FF"; col3 = "#F06719FF";col4 = "#1BA3C6FF"; col5 = "#F64971FF"; col6= "#F8B620FF"
source(here::here("code", "utilities.R"))
paletteer::paletteer_d("ggthemes::colorblind")
covars = read_rds(here::here("data", "processed", "covars_proc.rds"))
map_ci = read_rds(here::here("data", "analysis", "mapci_ranges.rds"))


if(!dir.exists(here::here("manuscript", "figures"))) dir.create(here::here("manuscript", "figures"))
if(!dir.exists(here::here("manuscript", "final_figures"))) dir.create(here::here("manuscript", "final_figures"))

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
            hi_map_hi_ci = sum(val_map >= 65 & val_ci > 2, na.rm = TRUE),
            .groups = "drop")
## map ci 


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

p = map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  summarize(across(contains("map"), 
                   list(mean = mean, 
                        sd = sd, 
                        n = ~ sum(.x >= 5, na.rm = TRUE), 
                        pctn = ~sum(.x >= 5, na.rm = TRUE) / n()))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(map = sub(".*map_([^_]*)_.*", "\\1", name),
         ci = sub(".*ci_([^_]*)_.*", "\\1", name),
         metric = sub(".*_(.*)", "\\1", name)) %>% 
  select(value, map, ci, metric) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(map = factor(map, levels = c("[0,64]", "(64,Inf]"), labels = c("<65", ">=65")),
         ci = factor(ci, levels = c("[0,2]", "(2,2.4]", "(2.4,2.8]", "(2.8,Inf]"),
                     labels =  c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8"))) %>% 
  mutate(fcol = factor(1:8)) %>% 
  ggplot(aes(x = map, y = ci, fill = fcol)) + 
  geom_tile(col = "black", alpha = .5) + 
  geom_text(aes(label = paste0(n, " (", round(pctn*100, 0), "%)")), size = 8) + 
  scale_fill_manual(values = c("#e51931", "#ca001b","#f37076", "#f37076",  
                               "#55caff", "#0f5fa5" ,"#137fc7", "#55caff"))  + 
  theme_classic(base_size = 20) +
  labs(x = "Mean Arterial Pressure (mmHg)", 
       y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8"))

png(here::here("manuscript", "figures", "n_pct_range.png"),
               width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

paletteer::paletteer_d("colorBlindness::Green2Magenta16Steps")
p =map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  summarize(across(contains("map"), 
                   list(mean = mean, 
                        sd = sd, 
                        n = ~ sum(.x >= 5, na.rm = TRUE), 
                        pctn = ~sum(.x >= 5, na.rm = TRUE) / n()))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(map = sub(".*map_([^_]*)_.*", "\\1", name),
         ci = sub(".*ci_([^_]*)_.*", "\\1", name),
         metric = sub(".*_(.*)", "\\1", name)) %>% 
  select(value, map, ci, metric) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(map = factor(map, levels = c("[0,64]", "(64,Inf]"), labels = c("<65", ">=65")),
         ci = factor(ci, levels = c("[0,2]", "(2,2.4]", "(2.4,2.8]", "(2.8,Inf]"),
                     labels =  c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8"))) %>% 
  mutate(fcol = factor(1:8)) %>% 
  ggplot(aes(x = map, y = ci, fill = fcol)) + 
  geom_tile(col = "black", alpha = .5) + 
  geom_text(aes(label = paste0(n, " (", round(pctn*100, 0), "%)")), size = 8) + 
  scale_fill_manual(values = c("#860086FF", "#500050FF","#FF50FFFF", "#FF50FFFF",  
                               "#86FF86FF", "#005000FF" ,"#00BB00FF", "#86FF86FF"))  + 
  theme_classic(base_size = 20) +
  labs(x = "Mean Arterial Pressure (mmHg)", 
       y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8"))

paletteer_c("ggthemes::Orange-Blue Diverging", 30)

p =map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  summarize(across(contains("map"), 
                   list(mean = mean, 
                        sd = sd, 
                        n = ~ sum(.x >= 5, na.rm = TRUE), 
                        pctn = ~sum(.x >= 5, na.rm = TRUE) / n()))) %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(map = sub(".*map_([^_]*)_.*", "\\1", name),
         ci = sub(".*ci_([^_]*)_.*", "\\1", name),
         metric = sub(".*_(.*)", "\\1", name)) %>% 
  select(value, map, ci, metric) %>% 
  pivot_wider(names_from = metric, values_from = value) %>% 
  mutate(map = factor(map, levels = c("[0,64]", "(64,Inf]"), labels = c("<65", ">=65")),
         ci = factor(ci, levels = c("[0,2]", "(2,2.4]", "(2.4,2.8]", "(2.8,Inf]"),
                     labels =  c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8"))) %>% 
  mutate(fcol = factor(1:8)) %>% 
  ggplot(aes(x = map, y = ci, fill = fcol)) + 
  scale_fill_paletteer_d("ggsci::grey_material") +
  geom_tile(col = "black", alpha = .5) + 
  geom_text(aes(label = paste0(n, " (", round(pctn*100, 0), "%)")), size = 8) + 
  # scale_fill_manual(values = c("#A94322FF", "#9E3D22FF","#BF4F22FF", "#BF4F22FF",  
  #                              "#5082B0FF", "#2B5C8AFF" ,"#3A6B99FF", "#5082B0FF"))  + 
  theme_classic(base_size = 20) +
  labs(x = "Mean Arterial Pressure (mmHg)", 
       y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8"))


png(here::here("manuscript", "figures", "n_pct_range_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

png(here::here("manuscript", "final_figures", "n_pct_range_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

p = 
  ci_ranges4 %>% 
  filter(id %in% wake_covars$id) %>% 
  pivot_longer(cols = contains("ci")) %>% 
  mutate(name = sub(".*ci\\_", "", name),
         name = factor(name, levels = c("[0,2]", "(2,2.4]", "(2.4,2.8]", "(2.8,Inf]"),
                       labels = c("<=2", "(2, 2.4]", "(2.4, 2.8]", ">2.8"))) %>% 
  ggplot(aes(x = name, y = value)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha = .3, size = .4, width = .1, color = "#0072B2FF") +
  labs(y = "Total Minutes in Range", x = expression("Quartile of Cardiac Index (L/min/" * m^2 * ")")) + 
  scale_y_continuous(breaks=seq(0, 360, 60), limits = c(0, 180)) + 
  theme_light(base_size = 20) + 
  theme(legend.position = "none") + 
  scale_x_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8"))
  

png(here::here("manuscript", "figures", "ci_quartiles_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 


png(here::here("manuscript", "final_figures", "ci_quartiles_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

p = hemo %>% 
  filter(cat_cpb != "intra") %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "post"),
                          labels = c("Pre-CPB", "Post-CPB"))) %>% 
  ggplot(aes(x = val_ci)) + 
  geom_histogram(
    aes(y = ..density..),
    binwidth = 0.5,
    color = "black",
    alpha = 0.7,
    fill = "darkgrey"
  ) +
  geom_density(size = 1, fill = NA, color = "#0072B2FF") + # was darkblue before
  theme_light(base_size = 20) +
  labs(x = expression("Cardiac Index (L/min/" * m^2 * ")"), y = "Density") +
  scale_x_continuous(breaks = seq(0, 10, 1), limits = c(0, 10)) +
  theme(legend.position = "none")  

png(here::here("manuscript", "figures", "ci_density.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

png(here::here("manuscript", "figures", "ci_density.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()



p = map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  pivot_longer(cols = -id) %>% 
  separate_wider_delim(name, "_", names = c("map", "ci")) %>% 
  ggplot(aes(x = ci, y = value, color = map)) + 
  geom_boxplot(outlier.size = .4, outlier.alpha = .5) + 
  labs(x = expression("Quartile of Cardiac Index (L/min/" * m^2 * ")"), y = "Total Minutes in Range Per Participant") + 
  scale_y_continuous(breaks=seq(0,360,60), limits=c(0,180)) + 
  scale_x_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_color_manual(values = c("#009E73FF","#D55E00FF"),labels = c("<65", ">=65"), name = "Mean Arterial Pressure (mmHg)") +
  theme_light(base_size = 20) + 
  theme(legend.position = c(0.28, .90))



p = map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  pivot_longer(cols = -id) %>% 
  separate_wider_delim(name, "_", names = c("map", "ci")) %>% 
  mutate(map = factor(map, levels = c("hypo", "normo"), labels = c("MAP <65", expression(MAP>=65))),
         ci = factor(ci, levels = c("q1", "q2", "q3", "q4"), labels = c(
           "CI <= 2",                       # parsed as math
           "\"CI (2, 2.4]\"",              # quoted so parse() yields a string literal
           "\"CI (2.4, 2.8]\"",
           "\"CI > 2.8\""
         )),
         # map = fct_rev(map),
         ci = fct_rev(ci)) %>% 
  ggplot(aes(x = 1, y = value)) + 
  facet_grid(ci ~ map, switch = "both",labeller = label_parsed) + 
  geom_boxplot(outlier.size = .5, outlier.alpha = .5, fill ="lightgrey") + 
  labs(x = "",  y = "Total Minutes in Range Per Participant") + 
  scale_y_continuous(breaks=seq(0,360,60), limits=c(0,120), position = "right") + 
  scale_x_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  # scale_color_manual(values = c("#009E73FF","#D55E00FF"),labels = c("<65", ">=65"), name = "Mean Arterial Pressure (mmHg)") +
  theme_light(base_size = 20) +
  theme(legend.position = "none")

map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  pivot_longer(cols = -id) %>% 
  separate_wider_delim(name, "_", names = c("map", "ci")) %>% 
  mutate(map = factor(map, levels = c("hypo", "normo"), labels = c("MAP <65", expression(MAP>=65))),
         ci = factor(ci, levels = c("q1", "q2", "q3", "q4"), labels = c(
           "CI <= 2",                       # parsed as math
           "\"CI (2, 2.4]\"",              # quoted so parse() yields a string literal
           "\"CI (2.4, 2.8]\"",
           "\"CI > 2.8\""
         )),
         # map = fct_rev(map),
         ci = fct_rev(ci)) %>% 
  filter(map == "MAP <65") %>% 
  ggplot(aes(x = 1, y = value)) + 
  facet_grid(ci ~ ., switch = "both",labeller = label_parsed) +
  geom_boxplot(outlier.size = .5, outlier.alpha = .5, outlier.shape=NA, fill ="lightgrey") + 
  labs(x = "",  y = "Total Minutes in Range Per Participant") + 
  # scale_y_continuous(breaks=seq(0,360,60), limits=c(0,120), position = "right") + 
  scale_x_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  # scale_color_manual(values = c("#009E73FF","#D55E00FF"),labels = c("<65", ">=65"), name = "Mean Arterial Pressure (mmHg)") +
  theme_light(base_size = 20) +
  theme(legend.position = "none")
png(here::here("manuscript", "figures", "ci_map_min_in_range_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

png(here::here("manuscript", "final_figures", "ci_map_min_in_range_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

p = map_ci_cpb %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4", "cat_cpb")) %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "post"), labels = c("Pre-CPB", "Post-CPB"))) %>% 
  pivot_longer(cols = -c(id, cat_cpb)) %>% 
  separate_wider_delim(name, "_", names = c("map", "ci")) %>% 
  ggplot(aes(x = ci, y = value, color = map)) + 
  facet_wrap(.~cat_cpb) + 
  geom_boxplot(outlier.size = .4, outlier.alpha = .5) + 
  labs( x = expression("Quartile of Cardiac Index (L/min/" * m^2 * ")"),
        y = "Total Minutes in Range Per Participant") + 
  scale_y_continuous(breaks=seq(0,360,60), limits = c(0,120)) + 
  scale_x_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_color_manual(values = c("#009E73FF","#D55E00FF"),labels = c("<65", ">=65"), name = "Mean Arterial Pressure (mmHg)") + 
  theme_light(base_size = 20) + 
  theme(legend.position = c(0.28, .90))

p=map_ci_cpb %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4", "cat_cpb")) %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "post"), labels = c("Pre-CPB", "Post-CPB"))) %>% 
  pivot_longer(cols = -c(id, cat_cpb)) %>% 
  separate_wider_delim(name, "_", names = c("map", "ci")) %>% 
  mutate(map = factor(map, levels = c("hypo", "normo"), labels = c("MAP <65", expression(MAP>=65))),
         ci = factor(ci, levels = c("q1", "q2", "q3", "q4"), labels = c(
           "CI <= 2",                       # parsed as math
           "\"CI (2, 2.4]\"",              # quoted so parse() yields a string literal
           "\"CI (2.4, 2.8]\"",
           "\"CI > 2.8\""
         )),
         ci = fct_rev(ci)) %>% 
  ggplot(aes(x = cat_cpb, y = value)) + 
  facet_grid(ci ~ map, switch = "both", labeller = label_parsed) + 
  geom_boxplot(outlier.size = .5,fill = "lightgrey", outlier.alpha = .5, aes(outlier.color = cat_cpb)) + 
  labs( x = "",
        y = "Total Minutes in Range Per Participant") + 
  scale_y_continuous(breaks=seq(0,360,60), limits = c(0,120), position = "right") +
  scale_x_discrete(position = "top") +
  # scale_fill_manual(values = c("#009E73FF","#D55E00FF"), name = "CPB Phase") + 
  theme_light(base_size = 20) + 
  theme(legend.position = "none")

  # theme(legend.position = c(0.28, .90))

png(here::here("manuscript", "figures", "ci_map_min_in_range_cpb_v2.png"),
      width = 10, height = 8, units = "in", res = 350)
p
dev.off() 
  

png(here::here("manuscript", "final_figures", "ci_map_min_in_range_cpb_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

  

df = 
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars) 

result_uni = 
  map_dfr(.x = df %>% select(contains("q")) %>% colnames,
          .f = function(x){
            formula <- as.formula(paste("bin_aki48h", "~", x))
            model <- glm(formula, data = df %>% mutate(across(contains("q"), ~.x / 5)), family = binomial)
            broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% slice(2)
          })

p = result_uni %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  scale_x_discrete(labels = c("<65", ">=65")) + 
  scale_y_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())


p = result_uni %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())


png(here::here("manuscript", "figures", "univariate_reg_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

png(here::here("manuscript", "final_figures", "univariate_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

p_ci = 
  result_uni %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
        lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
         hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE)) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_),
         ci = paste0("(", lo, ",", hi,")")) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", est, "\n", ci))) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

png(here::here("manuscript", "figures", "revision/univariate_reg_ci.png"),
    width = 10, height = 8, units = "in", res = 350)
p_ci
dev.off() 


model = glm(
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
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


p = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  scale_x_discrete(labels = c("<65", ">=65")) + 
  scale_y_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

p = model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_text(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p)), size = 12) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

png(here::here("manuscript", "figures", "adjusted_reg_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 

png(here::here("manuscript", "final_figures", "adjusted_reg_v2.jpg"),
    width = 12, height = 8, units = "in", res = 350)
p
dev.off() 

p_ci = model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
         lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
         hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE),
         ci = paste0("(", lo, ",", hi,")")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_text(aes(label = paste0("OR = ", est, "\n", ci)), size = 12) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")")) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

p_ci

png(here::here("manuscript", "figures/revision", "adjusted_reg_v2.jpg"),
    width = 12, height = 8, units = "in", res = 350)
p_ci
dev.off() 

scale_fill_manual(values = c("#e51931", "#ca001b","#f37076", "#f37076",  
                               "#55caff", "#0f5fa5" ,"#137fc7", "#55caff"))  


## adjusted sensitivity

low_lvef = wake_covars %>% 
  filter(bin_ef40 == "<=40%") %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% low_lvef) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("LVEF <=40, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))

low_lvef_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "LVEF <= 40")

mod1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  mutate(type = "low_lvef")


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% low_lvef)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    bin_stroke,
    bin_mi,
    bin_chf,
    bin_pvd,
    bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

hi_lvef_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "LVEF > 40")


p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("LVEF >40, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


### Shock groups 

shock = wake_covars %>% 
  filter(bin_iabp == 1 | bin_emergent == 1) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% shock) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    # bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    # bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  mutate(type = "shock")

shock_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "Shock")
p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("Shock, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% shock)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
    contains("q"),
    bin_aki48h,
    cat_gender,
    val_creatlst,
    val_age,
    val_bmi,
    bin_htn,
    bin_diabetes,
    # bin_stroke,
    bin_ef40,
    bin_mi,
    bin_chf,
    bin_pvd,
    # bin_cld,
    bin_betablocker,
    bin_acearb,
    bin_statin,
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
noshock_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "No shock")

p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("No shock, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))



m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>%
  mutate(type = "no shock")



### eGFR groups 

low_egfr = wake_covars %>% 
  filter(val_egfr < 60) %>% 
  pull(id)

df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(id %in% low_egfr) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

low_egfr_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "eGFR < 60")

m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  mutate(type = "low egfr")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("eGFR < 60, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))


df_temp =  
  map_ci %>% 
  filter(id %in% wake_covars$id) %>% 
  filter(!(id %in% low_egfr)) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
    # bin_redo,
    bin_emergent,
    val_crystalloid,
    val_cpbtime,
    cat_rbc,
    val_hematocrit
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)

hi_egfr_cv = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(term == "hypo_q1") %>% 
  mutate(type = "eGFR >= 60")

m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  mutate(type = "normal egfr")
p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  mutate(name = paste0("eGFR >= 60, n = ", nrow(df_temp), ", n cases = ", nrow(df_temp %>% filter(bin_aki48h == 1))))



### Forest plot 

fplot_df = 
  hi_egfr_cv %>% 
  bind_rows(low_egfr_cv, shock_cv, noshock_cv, low_lvef_cv, hi_lvef_cv) %>% 
  mutate(type = factor(type, levels = c("Shock", "No shock", "LVEF <= 40", "LVEF > 40",
                                        "eGFR < 60", "eGFR >= 60"))) %>% 
  mutate(group = 
           case_when(grepl("hock", type) ~ "shock",
                     grepl("eGFR", type) ~ "egfr",
                     .default = "lvef")) %>% 
  mutate(n = c(1005, 267, 110, 1162, 143, 1129),
         n_cases = c(109, 75, 18, 166, 30, 154)) %>% 
  mutate(x_num = 
           case_when(type == "Shock" ~ 1,
                     type == "No shock" ~ 2,
                     type == "LVEF <= 40" ~ 4,
                     type == "LVEF > 40" ~ 5,
                     type == "eGFR < 60" ~ 7,
                     .default = 8))



fplot_df = 
  hi_egfr_cv %>% 
  bind_rows(low_egfr_cv, shock_cv, noshock_cv, low_lvef_cv, hi_lvef_cv, not_anemic_cv,
            anemic_cv, not_htn_cv, htn_cv) %>% 
  mutate(type = factor(type, levels = c("Shock", "No shock", "LVEF <= 40", "LVEF > 40",
                                        "eGFR < 60", "eGFR >= 60", "Anemic","Not Anemic",
                        "History of hypertension", "No history of hypertension"))) %>% 
  mutate(group = 
           case_when(grepl("hock", type) ~ "shock",
                     grepl("eGFR", type) ~ "egfr",
                     grepl("CHF", type) ~ "chf",
                     grepl("nemic", type) ~ "anemia",
                     grepl("hyper", type) ~ "htn",
                     .default = "lvef")) %>% 
  mutate(n = c(1005, 267, 110, 1162, 143, 1129, 597, 675, 329, 943),
         n_cases = c(109, 75, 18, 166, 30, 154, 122, 62, 56, 128)) %>% 
  mutate(x_num = 
           case_when(type == "Shock" ~ 1,
                     type == "No shock" ~ 2,
                     type == "LVEF <= 40" ~ 4,
                     type == "LVEF > 40" ~ 5,
                     type == "eGFR < 60" ~ 7,
                     type == "eGFR >= 60" ~ 8,
                     type == "Anemic" ~ 10,
                     type == "Not Anemic" ~ 11,
                     type == "History of hypertension" ~ 13,
                     type == "No history of hypertension" ~ 14))


p = fplot_df %>% 
  mutate(pval = paste0("p=", format.pval(round(p.value, 3), digits = 2))) %>% 
  ggplot(aes(x = estimate, y = x_num, color = group)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = .2) + 
  geom_vline(aes(xintercept = 1), linetype = "dashed") + 
  labs(x = "Odds Ratio of AKI\nFor 5 additional minutes with low MAP and low CI", y = "") + 
  scale_color_manual(values= c("#0072B2FF", "#D55E00FF", "#CC79A7FF", "#009E73FF", "#E69F00FF"))+
  # scale_color_paletteer_d("ggthemes::colorblind", direction = -1) + 
  theme_classic(base_size = 20) + 
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(breaks=seq(0.8, 2.8, 0.2), limits = c(0.7, 2.8)) + 
  geom_label(aes(x = 2.55, y = x_num, hjust = 1, label = paste0("n=", n, "; cases=", n_cases)), 
             border.color = NA, size = 4) +
  geom_label(aes(x = estimate, y = x_num, label = pval),
            vjust = -1, size = 3.8, color = "black",border.color = NA) + 
  scale_y_continuous(breaks =seq(1:14), 
                     labels = c("Shock", "No shock","", 
                                expression(LVEF <= 40), "LVEF > 40", "", 
                                "eGFR < 60", expression(eGFR>=60), "",
                                "Anemic", "Not Anemic", "",
                                "History of hypertension", "No history of hypertension")) 

png(here::here("manuscript", "figures", "revision/forest_plot.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 


png(here::here("manuscript", "final_figures", "forest_plot.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off() 


### Surgery phase 


df_temp =  
  map_ci_pre %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
m1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  mutate(type = "pre cpb")

p1 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_),
         est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
                lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
                hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE),
                ci = paste0("(", lo, ",", hi,")")) %>% 
  mutate(name = paste0("Pre-CPB"))


df_temp =  
  map_ci_post %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3", "hypo_q4",
                           "normo_q1", "normo_q2", "normo_q3", "normo_q4")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)
m2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE,conf.int = TRUE) %>% 
  mutate(type = "post cpb")
p2 = 
  model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_),
       est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
                lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
                hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE),
                ci = paste0("(", lo, ",", hi,")")) %>% 
  mutate(name = "Post-CPB")


p = p1 %>% 
  bind_rows(p2) %>% 
  mutate(name = factor(name, levels = c("Pre-CPB", "Post-CPB"))) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  scale_x_discrete(labels = c("<65", ">=65")) + 
  facet_grid(.~name) + 
  scale_y_discrete(labels = c("<=2", "(2,2.4]", "(2.4,2.8]", ">2.8")) +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

p = p1 %>% 
  bind_rows(p2) %>% 
  mutate(name = factor(name, levels = c("Pre-CPB", "Post-CPB"))) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  facet_grid(.~name) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank()) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") 

  
png(here::here("manuscript", "figures", "cpb_reg_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

png(here::here("manuscript", "final_figures", "cpb_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()


p = p1 %>% 
  bind_rows(p2) %>% 
  mutate(name = factor(name, levels = c("Pre-CPB", "Post-CPB"))) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", est, "\n", ci))) + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  facet_grid(.~name) + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank()) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2), "(2,2.4]", "(2.4,2.8]", ">2.8")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)") 


png(here::here("manuscript", "figures/revision", "cpb_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

### Tertiles 



df_temp =  
  mapci_tert %>% 
  filter(id %in% wake_covars$id) %>% 
  magrittr::set_colnames(c("id", "hypo_q1", "hypo_q2", "hypo_q3",
                           "normo_q1", "normo_q2", "normo_q3")) %>% 
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


p = model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  scale_x_discrete(labels = c("<65", ">=65")) + 
  scale_y_discrete(labels = c("<=2.2", "(2.2,2.7]", ">2.7")) +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

p = model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2.2), "(2.2,2.7]", ">2.7")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)")  +
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

png(here::here("manuscript", "figures", "tertile_reg_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

png(here::here("manuscript", "final_figures", "tertile_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()


p = model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_),
         est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
         lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
         hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE),
         ci = paste0("(", lo, ",", hi,")")) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", est, "\n", ci))) + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=2.2), "(2.2,2.7]", ">2.7")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)")  +
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())


png(here::here("manuscript", "figures/revision", "tertile_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

### Quintiles 


df_temp =
  mapci_quint %>%
  filter(id %in% wake_covars$id) %>%
  magrittr::set_colnames(
    c(
      "id",
      "hypo_q1",
      "hypo_q2",
      "hypo_q3",
      "hypo_q4",
      "hypo_q5",
      "normo_q1",
      "normo_q2",
      "normo_q3",
      "normo_q4",
      "normo_q5"
    )
  ) %>%
  left_join(wake_covars)

model = glm(
  bin_aki48h ~ .,
  data = df_temp %>% select(
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
  ) %>% mutate(across(contains("q"), ~ .x / 5)),
  family = binomial()
)


p = model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  scale_x_discrete(labels = c("<65", ">=65")) + 
  scale_y_discrete(labels = c("<=1.9", "(1.9,2.2]", "(2.2,2.6]", "(2.6,2.9]", ">2.9")) +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

p = model %>% 
  broom::tidy(exponentiate = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_)) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p))) + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=1.9), "(1.9,2.2]", "(2.2,2.6]", "(2.6,2.9]", ">2.9")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)")  +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

png(here::here("manuscript", "figures", "quintile_reg_v2.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

png(here::here("manuscript", "final_figures", "quintile_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()



p = model %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(grepl("q", term)) %>% 
  separate_wider_delim(term, "_", names = c("MAP", "CI")) %>% 
  mutate(p = format.pval(p.value, digits = 2),
         est_sig = if_else(p.value < .05, estimate, NA_real_),
         est = format(signif(estimate, 3), scientific = FALSE, trim = TRUE), 
         lo = format(signif(conf.low, 3), scientific = FALSE, trim = TRUE),
         hi = format(signif(conf.high, 3), scientific = FALSE, trim = TRUE),
         ci = paste0("(", lo, ",", hi,")")) %>% 
  ggplot(aes(x = MAP, y = CI, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", est, "\n", ci))) + 
  labs(x = "Mean Arterial Pressure  (mmHg)", y = expression("Cardiac Index (L/min/" * m^2 * ")"))  + 
  scale_x_discrete(labels = c("<65", expression("">=65))) + 
  scale_y_discrete(labels = c(expression(""<=1.9), "(1.9,2.2]", "(2.2,2.6]", "(2.6,2.9]", ">2.9")) + 
  scale_fill_gradient2(low = "#0f5fa5", mid = "white", high = "#ca001b", midpoint = 1, name = "Odds Ratio (OR)")  +
  theme_light(base_size = 20) + 
  theme(panel.grid = element_blank())

png(here::here("manuscript", "figures/revision", "quintile_reg_v2.jpg"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

plot_df = 
  hemo %>% 
  filter(cat_cpb != "intra") %>% 
  group_by(id, cat_cpb) %>% 
  mutate(time = as.numeric(difftime(time, min(time), units = "secs"))) %>% 
  ungroup() %>% 
  mutate(cat_map = cut(val_map, breaks=c(0, 65, Inf), right = FALSE),
         # cat_ci = cut(val_ci, breaks=c(0, 2, 2.4, 2.8, Inf), right = TRUE)) %>% 
         cat_ci = cut(val_ci, breaks = c(0, 2, Inf), right = TRUE)) %>% 
  mutate(cat = paste0(cat_map, "_", cat_ci)) %>% 
  mutate(cat = case_when(
    grepl("NA", cat) ~ NA_character_,
    cat == "[0,65)_(0,2]" ~ "Low MAP, Low CI",
    cat == "[0,65)_(2,Inf]" ~ "Low MAP, Normal CI",
    cat == "[65,Inf)_(2,Inf]" ~ "Normal MAP, Normal CI",
    cat == "[65,Inf)_(0,2]" ~ "Normal MAP, Low CI",
    .default = NA_character_
  )) %>% 
  mutate(cat_cpb = factor(cat_cpb, levels = c("pre", "post"), labels= c("Pre-CPB", "Post-CPB")))

# paletteer_d("colorBlindness::PairedColor12Steps")


n_sub = unique(plot_df$id) %>% length

p = plot_df %>%
  mutate(cat = factor(cat, levels = c("Low MAP, Low CI", "Normal MAP, Low CI","Low MAP, Normal CI", "Normal MAP, Normal CI"))) %>% 
  mutate(cat = forcats::fct_na_value_to_level(cat),
         cat = forcats::fct_rev(cat)) %>% 
  ggplot(aes(x = time / 60, fill = cat, color = cat)) +
  geom_bar() +
  facet_wrap(~cat_cpb, scales = "free_x") +
  scale_x_continuous(breaks=seq(0, 480, 30), labels = seq(0, 8, .5),
                     limits = c(0, 4*60)) +
  theme_classic(base_size = 20) +
  labs(x = "Time (hr)", y = "Proportion") +
  scale_color_manual(values = c("#19B2FFFF", "#FFBF7FFF", "#FF99BFFF","#E51932FF", "#999999"), na.translate = TRUE, 
                     labels = c("Normal MAP, Normal CI",
                                "Low MAP, Normal CI",
                                "Normal MAP, Low CI",
                                "Low MAP, Low CI",
                                "Missing"), name = "") +
  scale_fill_manual(values = c("#19B2FFFF", "#FFBF7FFF", "#FF99BFFF","#E51932FF", "#999999"), na.translate = TRUE, 
                    labels = c("Normal MAP, Normal CI",
                               "Low MAP, Normal CI",
                               "Normal MAP, Low CI",
                               "Low MAP, Low CI",
                               "Missing"), name = "") +
  theme(legend.position = c(.84, .78),
        legend.title = element_blank()) + 
  scale_y_continuous(breaks = seq(0, n_sub, n_sub / 10), labels = seq(0, 1, .1))

png(here::here("manuscript", "figures", "lasagna_plot.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

## make a four cat matrix for later 
df = 
  df %>% 
  left_join(map_only, by = "id")
df2 = 
  df %>% 
  left_join(four_cats, by = "id") %>% 
  select(-map_65)




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

m1 %>% 
  broom::tidy(exp = TRUE, conf.int = TRUE) %>% 
  filter(grepl("map", term)) %>% 
  mutate(map = c("low", "low", "high", "high"),
         ci = c("low", "high", "low", "high")) %>% 
  ggplot(aes(x = map, y = ci, fill = estimate)) + 
  geom_tile(col = "black") + 
  geom_label(aes(label = paste0("OR = ", round(estimate, 3), "\np=", p.value))) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") 

extra = tibble(estimate =c(0.8766565, 0.8697321),
               term = c("inc_map", "inc_ci"))
               
m1 %>% 
  broom::tidy(exp = TRUE, conf.int = TRUE) %>% 
  filter(grepl("map", term)) %>% 
  select(term, estimate)  %>% 
  bind_rows(extra) %>% 
  ggplot(aes(x = term, y = estimate, color = estimate)) + 
  geom_point(size = 10) +
scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1, name = "Odds Ratio (OR)") 

  
