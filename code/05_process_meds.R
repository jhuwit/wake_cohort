# Script to calculate total vasotrope-inotrope score published by De La Hoz 2022
# in WAKE cohort.
# Binary yes/no for milrinone and dobutamine.
# Exclude fields with BUPIVACAINE (skin numbing agent).
# Only one participant with DOPAMINE. Exclude from calculation.
# No VASOPRESSIN in data set. Unique medications are PHENYLEPHRINE, EPINEPHRINE,
# NOREPINEPHRINE, DOBUTAMINE, MILRINONE, DOPAMINE, BUPIVACAINE

# NOTES:
# EPINEPHRINE mcg/kg/min is a mislabel. 
# Epic records doses in mcg/min and dose values align with expectations in mcg/min.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, readr, here, stringr)

dat = read_csv(here("data", "raw", "or_meds.csv"),
               col_types= cols(ID = col_character())) |>
  mutate(med = word(`Medication Administered`, 1)) |>
  mutate(`dose value unit` = ifelse(`dose value unit` == "mcg/kg/min" & med == "EPINEPHRINE", 
                                    "mcg/min", `dose value unit`),
         mins = `Seconds fromANStart` / 60,
         `Dose Value` = ifelse(`dose value unit` == "mg", `Dose Value` * 1000, `Dose Value`),
         `dose value unit` = ifelse(`dose value unit` == "mg", "mcg", `dose value unit`)) |>
  filter(!(med %in% c("BUPIVACAINE", "DOPAMINE")),
         `Action Name` != "Canceled Entry")

subject_ids = unique(dat$ID); head(subject_ids)

# Binary coding for MILRINONE and DOBUTAMINE
milrinone_ids = dat |> filter(med == "MILRINONE", mins >= 0) |> dplyr::select(ID) |> unique()
dobutamine_ids = dat |> filter(med == "DOBUTAMINE", mins >= 0) |> dplyr::select(ID) |> unique()
mil_dobu_df = tibble(ID = subject_ids,
                     milrinone = ifelse(ID %in% milrinone_ids$ID, 1, 0),
                     dobutamine = ifelse(ID %in% dobutamine_ids$ID, 1, 0))

# Exclude MILRINONE and DOBUTAMINE 
dat = dat |> filter(!(med %in% c("DOBUTAMINE", "MILRINONE")))

# Caculate total bolus doses
# "Given" and "Bolus" Action Names have the same dosing (mcg)
dat = dat |> mutate(`Action Name` = ifelse(`Action Name` == "Given", "Bolus", `Action Name`))
bolus_temp = dat |> filter(`Action Name` == "Bolus") |> 
  filter(mins >= 0) |> # Exclude boluses outside of anesthesia period
  dplyr::select(ID, `Dose Value`, `dose value unit`, med) |>
  # Conversion according to De La Hoz 2022
  mutate(converted_value = ifelse(med == "PHENYLEPHRINE", `Dose Value` / 10, `Dose Value`)) |>
  group_by(ID) |> summarize(bolus_total = sum(converted_value)) 
no_bolus_ids = data.frame(ID = subject_ids) |> anti_join(bolus_temp, by = "ID") |>
  mutate(bolus_total = 0)
bolus_df = rbind(bolus_temp, no_bolus_ids)

# Preprocess infusion doses
infusion_temp = dat |> 
  filter(`Action Name` %in% c("New Bag", "Stopped", "Restarted", "Rate Change", "Anesthesia Infusion Stopped"),
         # only infusion units for epi/norepi/phenylepi are mcg/min and mcg
         # mcg is not a valid EPIC entry, hence exclude
         `dose value unit` == "mcg/min") |> 
  # Remove duplicate infusion entries
  # distinct(`Medication Administered`, `Med Id`, `Action Name`, `Dose Value`, `dose value unit`, ORDER_MED_ID, mins, ID, .keep_all = TRUE)
  distinct(`Medication Administered`, mins, ID, .keep_all = TRUE) |>
  # Note that for one participant, they have concurrent NOREPINEPHRINE doses, which is impossible. 
  # The same medication and doses had ORDER_MED_ID == 300338 and 300339, but at slightly different minutes
  # Lee suggests excluding ORDER_MED_ID == 300339
  filter(ORDER_MED_ID != 300339) |>
  mutate(mins = ifelse(mins < 0, 0, mins))  |> # Truncate infusions for the time after anesthesia start
  arrange(ID, `Medication Administered`, mins) 

# Check for concurrent infusions of the same medication (impossible)
start_actions <- c("New Bag", "Restarted")
stop_actions  <- c("Stopped", "Anesthesia Infusion Stopped")

concurrency_trace <- infusion_temp |> 
  group_by(ID, `Medication Administered`) |>
  arrange(mins, .by_group = TRUE) |>
  mutate(
    delta = case_when(
      `Action Name` %in% start_actions ~  1L,   # starting a line
      `Action Name` %in% stop_actions  ~ -1L,   # stopping a line
      TRUE                      ~  0L    # rate changes etc.
    ),
    open_lines = cumsum(delta) # If open_lines > 1, that means there are concurrent lines
  ) |>
  ungroup()

# Any time with >1 open line = concurrency
concurrent_events <- concurrency_trace %>%
  filter(open_lines > 1)

concurrent_pairs <- concurrent_events %>%
  distinct(ID, `Medication Administered`)

n_concurrent_pairs <- nrow(concurrent_pairs)
n_concurrent_pairs

# Drill into one example:
example <- concurrent_pairs %>% slice(47)

subject_sample = dat %>%
  filter(ID == example$ID, `Medication Administered` == example$`Medication Administered`) %>%
  arrange(mins)

# After inspection, we notice that all the inspected concurrent doses seem to 
# have miscoded "Rate Change" as a "New Bag" or "Restarted"
# We have recoded these to "Rate Change"
infusion_with_prev <- infusion_temp %>%
  group_by(ID, `Medication Administered`) %>%
  arrange(mins, .by_group = TRUE) %>%
  mutate(
    # delta based on original actions
    delta_raw = case_when(
      `Action Name` %in% start_actions ~  1L,
      `Action Name` %in% stop_actions  ~ -1L,
      TRUE                      ~  0L
    ),
    open_lines_raw = cumsum(delta_raw),
    
    # how many lines were open *before* this event?
    open_lines_prev = open_lines_raw - delta_raw,   # or lag(open_lines_raw, default = 0)
    
    # fix: if New Bag happens while ≥1 line is already open, treat it as Rate Change
    action_fixed = case_when(
      `Action Name` == "New Bag" & open_lines_prev >= 1 ~ "Rate Change",
      `Action Name` == "Restarted" & open_lines_prev >= 1 ~ "Rate Change",
      TRUE                                      ~ `Action Name`
    )
  ) %>%
  ungroup()

infusion_fixed <- infusion_with_prev %>%
  select(-"delta_raw", -"open_lines_raw", -"open_lines_prev") %>%
  mutate(`Action Name` = action_fixed)


# Calculate infusion total dose
infusion_intervals <- infusion_fixed |>
  group_by(ID, `Medication Administered`) |>
  arrange(mins, .by_group = TRUE) |>
  mutate(
    start_min    = mins,
    end_min      = lead(mins),                       # time of next event
    duration_min = end_min - start_min
  ) |>
  # keep only proper intervals with positive rate
  filter(
    !is.na(duration_min),
    duration_min > 0,
    `Dose Value` > 0
  ) |>
  mutate(dose_sum = `Dose Value` * duration_min) |>
  # Conversion according to De La Hoz 2022
  mutate(converted_dose = ifelse(med == "PHENYLPHRENINE", dose_sum / 10, dose_sum)) |> 
  ungroup()

infusion_temp = infusion_intervals |> group_by(ID) |>
  summarize(infusion_total = sum(converted_dose)) |>
  ungroup()
no_infusion_ids = data.frame(ID = subject_ids) |> 
  anti_join(infusion_intervals, by = "ID") |>
  mutate(infusion_total = 0)
infusion_df = rbind(infusion_temp, no_infusion_ids)

result = infusion_df |> left_join(bolus_df, by = "ID") |>
  left_join(mil_dobu_df, by = "ID") |>
  mutate(total = infusion_total + bolus_total) |>
  dplyr::select(ID, total, dobutamine, milrinone)

write_rds(result, here("data", "processed", "vasopressor_inotrope_dose.rds"))

# temp = dat |> filter(`Action Name` == "Anesthesia Infusion Stopped") |> 
#   group_by(`Medication Administered`, ID)

# Check that "Given" action is bolus/injection (i.e. not infusion)
# length(which(dat$`Action Name` == "Given" & !(dat$`dose value unit` %in% c("mg", "mcg"))))

# Check infusion units
# which(dat$`Action Name` %in% c("New Bag", "Stopped", "Restarted", "Rate Change", "Anesthesia Infusion Stopped") & dat$`dose value unit` != "mcg/min")
