# ============================================================
# BUILD HEMODYNAMIC TIME-SERIES DATASET WITH OR EVENT TIMES
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Packages and options
# ------------------------------------------------------------

required_packages <- c(
  "readr",
  "readxl",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "lubridate",
  "here",
  "furrr",
  "future",
  "parallelly",
  "janitor",
  "tibble"
)

check_and_install_packages <- function(pkg_list) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Uncomment if packages need to be installed
# check_and_install_packages(required_packages)

library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(lubridate)
library(here)
library(furrr)
library(future)
library(parallelly)
library(janitor)
library(tibble)

options(
  digits.secs = 2,
  dplyr.summarise.inform = FALSE
)

# TRUE rebuilds processed files even when they already exist
force <- TRUE

# Always retain at least one worker
ncores <- max(
  1L,
  as.integer(parallelly::availableCores()) - 1L
)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

# Parse several common timestamp formats into UTC
parse_datetime_utc <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = "UTC"))
  }
  
  lubridate::parse_date_time(
    as.character(x),
    orders = c(
      "mdy HM",
      "mdy HMS",
      "ymd HM",
      "ymd HMS"
    ),
    tz = "UTC",
    quiet = TRUE
  )
}

# Select an appropriate future plan
set_parallel_plan <- function(n_workers) {
  if (n_workers > 1L) {
    future::plan(
      future::multisession,
      workers = n_workers
    )
  } else {
    future::plan(future::sequential)
  }
}

# ------------------------------------------------------------
# 3. File locations
# ------------------------------------------------------------

covariate_file <- here::here(
  "data",
  "raw",
  "WAKE.flatfile.IDsfixed.11.27.24.xlsx"
)

or_times_file <- here::here(
  "data",
  "raw",
  "WAKE_OR_Times.IDfixed.11.29.24.csv"
)

individual_directory <- here::here(
  "data",
  "raw",
  "individual"
)

processed_directory <- here::here(
  "data",
  "processed"
)

raw_timeseries_file <- here::here(
  processed_directory,
  "all_timeseries_stamped_raw.rds"
)

hemo_timeseries_file <- here::here(
  processed_directory,
  "hemo_timeseries.rds"
)

covariates_with_times_file <- here::here(
  processed_directory,
  "covariates_with_or_times.rds"
)

dir.create(
  processed_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# Confirm that the required input files exist
stopifnot(file.exists(covariate_file))
stopifnot(file.exists(or_times_file))
stopifnot(dir.exists(individual_directory))

# ------------------------------------------------------------
# 4. Read covariate data
# ------------------------------------------------------------

covars <- readxl::read_excel(covariate_file) %>%
  janitor::clean_names() %>%
  mutate(id = as.character(id))

# Patient IDs occurring more than once
covar_duplicate_ids <- covars %>%
  count(id, name = "n") %>%
  filter(n > 1L)

if (nrow(covar_duplicate_ids) > 0L) {
  message(
    nrow(covar_duplicate_ids),
    " IDs occur more than once in the covariate dataset."
  )
}

# ------------------------------------------------------------
# 5. Read and process OR event times
# ------------------------------------------------------------

or_times <- readr::read_csv(
  or_times_file,
  show_col_types = FALSE
) %>%
  janitor::clean_names() %>%
  mutate(id = as.character(id))

or_times_small <- or_times %>%
  transmute(
    id = id,
    
    # Anesthesia events
    anes_ready = anesthesia_ready,
    anes_start = anesthesia_start,
    anes_end = anesthesia_stop,
    
    # Candidate surgical-start events
    surgery_start_incision = incision_time,
    surgery_start_no_incision = proc_start_no_incision,
    
    # Cardiopulmonary bypass events
    cpb_start = cv_bypass_initiated,
    cpb_end = cv_bypass_ended
  ) %>%
  mutate(
    across(
      -id,
      parse_datetime_utc
    ),
    
    # Incision is preferred. For a procedure without an incision,
    # use "Proc Start/No Incision" when available.
    surgery_start = coalesce(
      surgery_start_incision,
      surgery_start_no_incision
    ),
    
    surgery_start_source = case_when(
      !is.na(surgery_start_incision) ~ "Incision Time",
      is.na(surgery_start_incision) &
        !is.na(surgery_start_no_incision) ~
        "Proc Start/No Incision",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 6. OR-time quality checks
# ------------------------------------------------------------

or_duplicate_ids <- or_times_small %>%
  count(id, name = "n") %>%
  filter(n > 1L)

if (nrow(or_duplicate_ids) > 0L) {
  print(or_duplicate_ids)
  
  stop(
    paste0(
      "Some IDs have multiple rows in the OR-time dataset. ",
      "Joining only by ID could attach the wrong operation. ",
      "Add an operation date or encounter identifier before continuing."
    )
  )
}

missing_anesthesia_times <- or_times_small %>%
  filter(is.na(anes_start) | is.na(anes_end))

invalid_anesthesia_times <- or_times_small %>%
  filter(
    !is.na(anes_start),
    !is.na(anes_end),
    anes_end < anes_start
  )

invalid_surgery_times <- or_times_small %>%
  filter(
    !is.na(surgery_start),
    !is.na(anes_start),
    !is.na(anes_end),
    surgery_start < anes_start |
      surgery_start > anes_end
  )

invalid_cpb_times <- or_times_small %>%
  filter(
    !is.na(cpb_start),
    !is.na(cpb_end),
    cpb_end < cpb_start
  )

message(
  nrow(missing_anesthesia_times),
  " OR records have missing anesthesia start or end."
)

message(
  nrow(invalid_anesthesia_times),
  " OR records have anesthesia end before anesthesia start."
)

message(
  nrow(or_times_small %>% filter(is.na(surgery_start))),
  " OR records have no surgery-start timestamp."
)

message(
  nrow(invalid_surgery_times),
  " OR records have surgery start outside the anesthesia interval."
)

message(
  nrow(invalid_cpb_times),
  " OR records have CPB end before CPB start."
)

# ------------------------------------------------------------
# 7. Add OR times to the patient-level covariate dataset
# ------------------------------------------------------------

covars_with_times <- covars %>%
  left_join(
    or_times_small %>%
      select(
        id,
        anes_ready,
        anes_start,
        surgery_start,
        surgery_start_source,
        cpb_start,
        cpb_end,
        anes_end
      ),
    by = "id"
  ) %>%
  mutate(
    minutes_anes_start_to_surgery = as.numeric(
      difftime(
        surgery_start,
        anes_start,
        units = "mins"
      )
    ),
    
    anesthesia_duration_min = as.numeric(
      difftime(
        anes_end,
        anes_start,
        units = "mins"
      )
    ),
    
    cpb_duration_min = as.numeric(
      difftime(
        cpb_end,
        cpb_start,
        units = "mins"
      )
    )
  )

readr::write_rds(
  covars_with_times,
  covariates_with_times_file,
  compress = "xz"
)

# ------------------------------------------------------------
# 8. Locate individual hemodynamic files
# ------------------------------------------------------------

individual_files <- list.files(
  individual_directory,
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(individual_files) == 0L) {
  stop("No CSV files were found in the individual-data directory.")
}

# Extract the ID from the filename
ts_ids <- individual_files %>%
  basename() %>%
  tools::file_path_sans_ext() %>%
  as.character()

# Detect duplicate filenames/IDs
duplicate_file_ids <- tibble(id = ts_ids) %>%
  count(id, name = "n") %>%
  filter(n > 1L)

if (nrow(duplicate_file_ids) > 0L) {
  print(duplicate_file_ids)
  
  stop(
    paste0(
      "Multiple individual CSV files resolve to the same ID. ",
      "Check the filenames before combining them."
    )
  )
}

# ------------------------------------------------------------
# 9. Compare ID coverage
# ------------------------------------------------------------

covar_ids <- unique(covars$id)
or_ids <- unique(or_times_small$id)

id_coverage <- tibble(
  comparison = c(
    "Covariate IDs absent from OR data",
    "OR IDs absent from covariate data",
    "Covariate IDs absent from individual files",
    "Individual-file IDs absent from covariate data",
    "Individual-file IDs absent from OR data"
  ),
  number = c(
    sum(!covar_ids %in% or_ids),
    sum(!or_ids %in% covar_ids),
    sum(!covar_ids %in% ts_ids),
    sum(!ts_ids %in% covar_ids),
    sum(!ts_ids %in% or_ids)
  )
)

print(id_coverage)

# ------------------------------------------------------------
# 10. Read and combine individual hemodynamic files
# ------------------------------------------------------------

if (!file.exists(raw_timeseries_file) || force) {
  
  set_parallel_plan(ncores)
  
  all_raw <- furrr::future_map_dfr(
    individual_files,
    function(file_path) {
      
      subject_id <- file_path %>%
        basename() %>%
        tools::file_path_sans_ext() %>%
        as.character()
      
      readr::read_csv(
        file_path,
        show_col_types = FALSE
      ) %>%
        janitor::clean_names() %>%
        select(
          timestamp,
          contains("anesthesia"),
          ci,
          map,
          cvp
        ) %>%
        mutate(
          id = subject_id,
          timestamp = parse_datetime_utc(timestamp)
        ) %>%
        relocate(id, timestamp)
    },
    .options = furrr::furrr_options(
      packages = c(
        "readr",
        "dplyr",
        "lubridate",
        "janitor"
      )
    )
  )
  
  future::plan(future::sequential)
  
  readr::write_rds(
    all_raw,
    raw_timeseries_file,
    compress = "xz"
  )
}

hemo_ts_df <- readr::read_rds(raw_timeseries_file) %>%
  mutate(
    id = as.character(id),
    timestamp = parse_datetime_utc(timestamp)
  )

# ------------------------------------------------------------
# 11. Check duplicate hemodynamic timestamps
# ------------------------------------------------------------

duplicate_hemo_times <- hemo_ts_df %>%
  count(id, timestamp, name = "n") %>%
  filter(
    !is.na(timestamp),
    n > 1L
  )

if (nrow(duplicate_hemo_times) > 0L) {
  message(
    nrow(duplicate_hemo_times),
    " ID-timestamp combinations occur more than once."
  )
  
  # Do not automatically aggregate duplicate measurements.
  # Review duplicate_hemo_times before deciding whether to use
  # the first value, mean, median, or another rule.
}

# ------------------------------------------------------------
# 12. Function to create one patient's full time series
# ------------------------------------------------------------

create_full_ts <- function(subject) {
  
  times <- or_times_small %>%
    filter(id == subject)
  
  # ID not present in OR dataset
  if (nrow(times) == 0L) {
    return(NULL)
  }
  
  # Multiple operations cannot safely be resolved using ID alone
  if (nrow(times) > 1L) {
    warning(
      "ID ", subject,
      " has multiple OR records and was skipped."
    )
    return(NULL)
  }
  
  # Anesthesia boundaries are required
  if (
    is.na(times$anes_start[1]) ||
    is.na(times$anes_end[1])
  ) {
    warning(
      "ID ", subject,
      " has missing anesthesia start/end and was skipped."
    )
    return(NULL)
  }
  
  # Reject reversed anesthesia timestamps
  if (times$anes_end[1] < times$anes_start[1]) {
    warning(
      "ID ", subject,
      " has anesthesia end before anesthesia start and was skipped."
    )
    return(NULL)
  }
  
  subject_ts <- hemo_ts_df %>%
    filter(id == subject) %>%
    select(
      id,
      timestamp,
      ci,
      map,
      cvp
    )
  
  # Construct one timestamp per minute from anesthesia start to end
  timeseq <- seq(
    from = floor_date(times$anes_start[1], unit = "minute"),
    to = floor_date(times$anes_end[1], unit = "minute"),
    by = "1 min"
  )
  
  tibble(
    id = subject,
    time = timeseq,
    
    # Repeat event timestamps on every patient-minute row
    anes_ready = times$anes_ready[1],
    anes_start = times$anes_start[1],
    surgery_start = times$surgery_start[1],
    surgery_start_source = times$surgery_start_source[1],
    cpb_start = times$cpb_start[1],
    cpb_end = times$cpb_end[1],
    anes_end = times$anes_end[1]
  ) %>%
    mutate(
      minutes_from_anes_start = as.numeric(
        difftime(
          time,
          anes_start,
          units = "mins"
        )
      ),
      
      minutes_from_surgery_start = as.numeric(
        difftime(
          time,
          surgery_start,
          units = "mins"
        )
      ),
      
      anesthesia_phase = case_when(
        time >= anes_start & time <= anes_end ~
          "during-anesthesia",
        TRUE ~ NA_character_
      ),
      
      surgery_phase = case_when(
        is.na(surgery_start) ~
          "surgery-start-missing",
        time < surgery_start ~
          "pre-surgery",
        time >= surgery_start ~
          "surgery-started"
      ),
      
      cpb_phase = case_when(
        is.na(cpb_start) | is.na(cpb_end) ~
          "CPB-times-missing",
        time < cpb_start ~
          "pre-CPB",
        time >= cpb_start & time <= cpb_end ~
          "intra-CPB",
        time > cpb_end ~
          "post-CPB"
      )
    ) %>%
    left_join(
      subject_ts,
      by = c(
        "id",
        "time" = "timestamp"
      )
    )
}

# ------------------------------------------------------------
# 13. Create full minute-level dataset
# ------------------------------------------------------------

if (!file.exists(hemo_timeseries_file) || force) {
  
  # Only process IDs represented in the raw hemodynamic data
  ids <- unique(hemo_ts_df$id)
  
  set_parallel_plan(ncores)
  
  full_timeseries <- furrr::future_map(
    ids,
    create_full_ts,
    .options = furrr::furrr_options(
      packages = c(
        "dplyr",
        "tibble",
        "lubridate"
      )
    )
  )
  
  future::plan(future::sequential)
  
  full_timeseries_df <- full_timeseries %>%
    purrr::compact() %>%
    bind_rows()
  
  readr::write_rds(
    full_timeseries_df,
    hemo_timeseries_file,
    compress = "xz"
  )
}

# ------------------------------------------------------------
# 14. Load and inspect the completed dataset
# ------------------------------------------------------------

full_timeseries_df <- readr::read_rds(
  hemo_timeseries_file
)

glimpse(full_timeseries_df)

# Review the event timestamps and derived phases
full_timeseries_df %>%
  select(
    id,
    time,
    anes_start,
    surgery_start,
    surgery_start_source,
    cpb_start,
    cpb_end,
    anes_end,
    minutes_from_anes_start,
    minutes_from_surgery_start,
    surgery_phase,
    cpb_phase,
    ci,
    map,
    cvp
  ) %>%
  print(n = 20)

# ------------------------------------------------------------
# 15. Final quality-control summaries
# ------------------------------------------------------------

# Availability of surgery start
full_timeseries_df %>%
  distinct(
    id,
    surgery_start,
    surgery_start_source
  ) %>%
  count(
    surgery_start_source,
    name = "patients"
  ) %>%
  print()

# Anesthesia start-to-surgery start distribution
covars_with_times %>%
  summarise(
    patients = n(),
    surgery_start_available = sum(
      !is.na(surgery_start)
    ),
    surgery_start_missing = sum(
      is.na(surgery_start)
    ),
    median_minutes_anes_to_surgery = median(
      minutes_anes_start_to_surgery,
      na.rm = TRUE
    ),
    q1_minutes_anes_to_surgery = quantile(
      minutes_anes_start_to_surgery,
      0.25,
      na.rm = TRUE
    ),
    q3_minutes_anes_to_surgery = quantile(
      minutes_anes_start_to_surgery,
      0.75,
      na.rm = TRUE
    )
  ) %>%
  print()