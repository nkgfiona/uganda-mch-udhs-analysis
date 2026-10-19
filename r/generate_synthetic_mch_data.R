# generate_synthetic_mch_data.R
# ---------------------------------------------------------------------------
# Generates synthetic mother/child records for two survey years (2016, 2022),
# calibrated to REAL, cited Uganda maternal & child health indicators.
#
# WHAT HAPPENS IN THIS SCRIPT
# Using the UDHS regional tables (data/raw/udhs_broad_region_targets.csv), this
# script calibrates EACH REGION x YEAR cell to its real published rate, and
# uses the real births N of each cell as its sample size (so regional shares
# are the real ones too).
#
# WHAT IS REAL vs SYNTHETIC
#   REAL (from UDHS 2016 / 2022 published tables, hand-extracted by you):
#     - Skilled birth attendance (SBA) rate for each of 4 regions x 2 years
#     - Stunting rate for each of 4 regions x 2 years (Table 11.1, 0-59 months)
#     - Births N per region x year (used as the synthetic sample size)
#   REAL (cited national cross-checks, used for verification only):
#     - SBA 74.2% (2016) -> 88.4% (2022)        Source: UBOS SDG Data Fact Sheet
#     - Stunting 28.9% (2016) -> 24.4% (2022)    Source: UDHS Table 11.1 "Total" row
#     - ANC4+ visits 60.3% (2016)                Source: PLOS One, "Continuum of
#                                                care for maternal health in Uganda"
#     - ANC4+ visits 67.8% (2022)                Source: UDHS 2022 Table 9.2
#   LITERATURE-ADAPTED (NOT Uganda-specific estimates):
#     - Effects of wealth, mother's education, urban residence and ANC4+ on SBA,
#       and of wealth/education/urban residence on stunting. Adapted from a
#       comparable published DHS analysis (Nigeria NDHS 2018: aOR ~6.27 for
#       wealth, ~7.01 for education, richest/highest vs poorest/none). Uganda's
#       own coefficients need DHS microdata, which cannot legally be
#       redistributed (DHS Terms of Use: "micro-level data will not be
#       re-distributed, either directly or within any tool/dashboard").
#   PURELY SYNTHETIC:
#     - Every individual row, and how wealth/education/residence are distributed
#       inside a region (identical across regions by construction).
#
# Run from the r/ folder. Writes ../data/raw/mch_synthetic_data.csv
# ---------------------------------------------------------------------------

set.seed(42)

RAW_DIR      <- "../data/raw"
targets_path <- file.path(RAW_DIR, "udhs_broad_region_targets.csv")
if (!file.exists(targets_path)) {
  stop("Can't find ", targets_path, ". Run from the r/ folder, and put ",
       "udhs_broad_region_targets.csv in data/raw/.")
}

REGIONS        <- c("Central", "Eastern", "Northern", "Western")
region_targets <- read.csv(targets_path, stringsAsFactors = FALSE)
stopifnot(all(region_targets$region %in% REGIONS),
          all(c("survey_year", "region", "sba_pct", "stunted_pct", "births_n") %in% names(region_targets)))

# National cross-check values (verification only, not used to generate data)
NATIONAL_CHECK <- list("2016" = list(sba = 74.2, stunting = 28.9),
                       "2022" = list(sba = 88.4, stunting = 24.4))
# stunting figures are the Table 11.1 "Total" row (0-59 months), both years

# ANC4+ has no regional target, so it is calibrated nationally per year.
# 2016: PLOS One continuum-of-care study. 2022: UDHS 2022 Table 9.2.
TARGETS_ANC4 <- c("2016" = 0.599, "2022" = 0.678)

# ---- Helpers ---------------------------------------------------------------

# Finds the intercept that makes the mean predicted probability equal a target.
calibrate_intercept <- function(linear_predictor_no_intercept, target_rate) {
  f <- function(b0) mean(plogis(linear_predictor_no_intercept + b0)) - target_rate
  uniroot(f, interval = c(-15, 15))$root
}

# Calibrates a SEPARATE intercept for each region so that every region's mean
# predicted probability equals its real published rate. `targets` is a named
# vector of proportions (names = regions).
calibrate_by_region <- function(lp, region, targets) {
  p <- numeric(length(lp))
  for (r in names(targets)) {
    idx  <- region == r
    b0   <- calibrate_intercept(lp[idx], targets[[r]])
    p[idx] <- plogis(lp[idx] + b0)
  }
  p
}

# For each value of `groups`, draws a category using that group's probabilities.
draw_by_group <- function(groups, prob_list, categories) {
  out <- character(length(groups))
  for (g in names(prob_list)) {
    idx <- groups == g
    out[idx] <- sample(categories, sum(idx), replace = TRUE, prob = prob_list[[g]])
  }
  out
}

# ---- Effect sizes (log-odds), literature-adapted (see header) ----------------
wealth_levels <- c("Poorest", "Poorer", "Middle", "Richer", "Richest")
educ_levels   <- c("None", "Primary", "Secondary", "Higher")

wealth_effect <- c(Poorest = 0, Poorer = 0.4, Middle = 0.9, Richer = 1.4, Richest = 1.84)
educ_effect   <- c(None = 0, Primary = 0.5, Secondary = 1.1, Higher = 1.95)
urban_effect  <- c(Urban = 0.6, Rural = 0)
anc4_effect   <- c(Yes = 1.33, No = 0)

# Wealth quintile skewed by residence (urban skews richer); education skewed by
# wealth. These reflect the known structure of DHS wealth-index data, not
# restricted microdata.
wealth_by_residence <- list(Urban = c(0.05, 0.10, 0.20, 0.30, 0.35),
                            Rural = c(0.28, 0.25, 0.20, 0.15, 0.12))
educ_by_wealth <- list(
  Poorest = c(0.35, 0.45, 0.18, 0.02),
  Poorer  = c(0.25, 0.45, 0.25, 0.05),
  Middle  = c(0.15, 0.40, 0.35, 0.10),
  Richer  = c(0.08, 0.30, 0.42, 0.20),
  Richest = c(0.03, 0.15, 0.42, 0.40)
)

# ---- One survey year -------------------------------------------------------
generate_cohort <- function(year) {
  rt <- region_targets[region_targets$survey_year == year, ]
  stopifnot(nrow(rt) == length(REGIONS))

  # Real births N per region = number of synthetic records in that region.
  region <- sample(rep(rt$region, times = rt$births_n))
  n      <- length(region)

  residence        <- sample(c("Urban", "Rural"), n, replace = TRUE, prob = c(0.25, 0.75))
  wealth_quintile  <- draw_by_group(residence, wealth_by_residence, wealth_levels)
  mother_education <- draw_by_group(wealth_quintile, educ_by_wealth, educ_levels)
  mother_age       <- pmin(49, pmax(15, round(rnorm(n, mean = 27, sd = 6))))
  child_sex        <- sample(c("Male", "Female"), n, replace = TRUE)
  child_age_months <- sample(0:59, n, replace = TRUE)

  w <- unname(wealth_effect[wealth_quintile])
  e <- unname(educ_effect[mother_education])
  u <- unname(urban_effect[residence])

  # ---- ANC4+ (national calibration per year, real target both years) ----
  lp_anc4     <- w + e + u
  anc4_target <- TARGETS_ANC4[[as.character(year)]]
  p_anc4      <- plogis(lp_anc4 + calibrate_intercept(lp_anc4, anc4_target))
  anc4_visits <- ifelse(runif(n) < p_anc4, "Yes", "No")

  # ---- Skilled birth attendance: each region calibrated to its REAL rate ----
  lp_sba  <- w + e + u + unname(anc4_effect[anc4_visits])
  p_sba   <- calibrate_by_region(lp_sba, region, setNames(rt$sba_pct / 100, rt$region))
  skilled_birth_attendance <- ifelse(runif(n) < p_sba, "Yes", "No")

  # ---- Stunting: each region calibrated to its REAL rate. Deliberately NOT
  # linked to SBA: stunting reflects longer-term nutrition, not the birth event.
  lp_stunt <- -(w * 0.6) - (e * 0.5) - (u * 0.5)
  p_stunt  <- calibrate_by_region(lp_stunt, region, setNames(rt$stunted_pct / 100, rt$region))
  stunted  <- ifelse(runif(n) < p_stunt, "Yes", "No")

  data.frame(
    record_id = paste0("MCH", year, sprintf("%05d", seq_len(n))),
    survey_year = year,
    region, residence, wealth_quintile, mother_education, mother_age,
    child_sex, child_age_months,
    anc4_visits, skilled_birth_attendance, stunted,
    data_origin = "SYNTHETIC",
    stringsAsFactors = FALSE
  )
}

full_data <- rbind(generate_cohort(2016), generate_cohort(2022))

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(full_data, file.path(RAW_DIR, "mch_synthetic_data.csv"), row.names = FALSE)
cat(sprintf("\nWrote %d synthetic records to %s\n", nrow(full_data), file.path(RAW_DIR, "mch_synthetic_data.csv")))

# ---- Verification ----------------------------------------------------------
cat("\n=== Verification 1: region x year rates vs REAL UDHS targets ===\n")
cat("(Differences of about +/-1-2 points are sampling noise from drawing random records.)\n")
tmp <- full_data
tmp$sba   <- as.numeric(tmp$skilled_birth_attendance == "Yes")
tmp$stunt <- as.numeric(tmp$stunted == "Yes")
sim <- aggregate(cbind(sba, stunt) ~ survey_year + region, data = tmp, FUN = function(x) 100 * mean(x))
ver <- merge(sim, region_targets[, c("survey_year", "region", "sba_pct", "stunted_pct", "births_n")],
             by = c("survey_year", "region"))
ver <- data.frame(year = ver$survey_year, region = ver$region, n = ver$births_n,
                  sba_sim = round(ver$sba, 1), sba_real = ver$sba_pct, sba_diff = round(ver$sba - ver$sba_pct, 1),
                  stunt_sim = round(ver$stunt, 1), stunt_real = ver$stunted_pct,
                  stunt_diff = round(ver$stunt - ver$stunted_pct, 1))
print(ver, row.names = FALSE)

cat("\n=== Verification 2: national rates vs cited national figures ===\n")
for (yr in c(2016, 2022)) {
  sub <- tmp[tmp$survey_year == yr, ]
  nat <- NATIONAL_CHECK[[as.character(yr)]]
  cat(sprintf("  %d  SBA: %.1f%% (cited %.1f%%)   Stunting: %.1f%% (cited %.1f%%)\n",
              yr, 100 * mean(sub$sba), nat$sba, 100 * mean(sub$stunt), nat$stunting))
  anc_t <- TARGETS_ANC4[[as.character(yr)]]
  cat(sprintf("       ANC4+: %.1f%% (cited %.1f%%)\n", 100 * mean(sub$anc4_visits == "Yes"), 100 * anc_t))
}

cat("\n=== Verification 3: did the wealth gradient come through? (SBA rate by quintile) ===\n")
print(round(100 * tapply(tmp$sba, factor(tmp$wealth_quintile, levels = wealth_levels), mean), 1))
