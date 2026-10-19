# mch_aggregate_analysis.R
# ---------------------------------------------------------------------------
# What this script does, in order:
#   1. Loads both CSVs and VALIDATES your transcription (internal consistency)
#   2. Reshapes to long format: region / year / indicator / value / n
#   3. Builds a sub-region crosswalk (the 2016 and 2022 sub-region schemes differ)
#   4. Tests whether each 2016 -> 2022 change is statistically detectable
#   5. Draws the exploratory figures (titles are computed from the data)
#   6. Does SBA and stunting move together across sub-regions?
#   7. Saves tidy outputs for the R Markdown report and Power BI
#
# Run from the r/ folder (RStudio: Session > Set Working Directory > To Source
# File Location). Every number in this script is REAL published aggregate data;
# nothing here is synthetic.
# ---------------------------------------------------------------------------

# The four packages below are all part of tidyverse; library(tidyverse) works too.
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

RAW  <- "../data/raw"
PROC <- "../data/processed"
FIG  <- "../documentation/figures"
for (d in c(PROC, FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

region_path <- file.path(RAW, "udhs_broad_region_targets.csv")
sub_path    <- file.path(RAW, "udhs_subregion_indicators.csv")
if (!file.exists(region_path) || !file.exists(sub_path)) {
  stop("Can't find the UDHS CSVs in ", normalizePath(RAW, mustWork = FALSE),
       ". Are you running from the r/ folder, and are both files in data/raw/?")
}

region_raw <- read_csv(region_path, show_col_types = FALSE)
sub_raw    <- read_csv(sub_path,    show_col_types = FALSE)

# Published NATIONAL figures used only as a cross-check (same sources you cited
# in generate_synthetic_mch_data.R: UBOS SDG Data Fact Sheet; UDHS 2022 Key Findings)
PUBLISHED_NATIONAL <- tibble(
  survey_year = c(2016, 2022),
  sba_pub     = c(74.2, 88.4),
  stunt_pub   = c(28.9, 24.4)   # Table 11.1 "Total" row, both years (0-59 months)
)

# ---------------------------------------------------------------------------
# 1. VALIDATE THE TRANSCRIPTION
# ---------------------------------------------------------------------------
cat("\n=== 1. Validating hand-extracted data ===\n")
results <- logical(0)
check <- function(label, ok) {
  ok <- isTRUE(ok)
  results <<- c(results, ok)
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", label))
}

check("broad-region file has the expected columns",
      all(c("survey_year", "region", "sba_pct", "stunted_pct", "births_n", "region_share") %in% names(region_raw)))
check("sub-region file has the expected columns",
      all(c("survey_year", "subregion", "broad_region", "skilled_birth_attendance_pct",
            "births_n", "stunted_pct", "stunting_children_n") %in% names(sub_raw)))
check("all percentages are between 0 and 100",
      all(c(region_raw$sba_pct, region_raw$stunted_pct,
            sub_raw$skilled_birth_attendance_pct, sub_raw$stunted_pct) >= 0 &
          c(region_raw$sba_pct, region_raw$stunted_pct,
            sub_raw$skilled_birth_attendance_pct, sub_raw$stunted_pct) <= 100))
check("no duplicate year x region rows", !any(duplicated(region_raw[c("survey_year", "region")])))
check("no duplicate year x sub-region rows", !any(duplicated(sub_raw[c("survey_year", "subregion")])))
check("sub-region stunting children N: complete for both years (Table 11.1 basis)",
      !any(is.na(sub_raw$stunting_children_n)))

# Roll the sub-region rows up to broad regions and compare with the broad file.
sub_rollup <- sub_raw %>%
  mutate(n_stunt = coalesce(stunting_children_n, births_n)) %>%
  group_by(survey_year, region = broad_region) %>%
  summarise(births_n_sub = sum(births_n),
            sba_sub      = weighted.mean(skilled_birth_attendance_pct, births_n),
            stunt_sub    = weighted.mean(stunted_pct, n_stunt),
            children_n   = sum(stunting_children_n),   # NA in 2022 by design
            .groups = "drop")

recon <- inner_join(region_raw, sub_rollup, by = c("survey_year", "region"))
check("sub-region births add up EXACTLY to broad-region births", all(recon$births_n == recon$births_n_sub))
check("births-weighted sub-region SBA reproduces broad-region SBA (within 0.15 pp)",
      all(abs(recon$sba_sub - recon$sba_pct) < 0.15))
check("weighted sub-region stunting reproduces broad-region stunting (within 0.1 pp)",
      all(abs(recon$stunt_sub - recon$stunted_pct) < 0.1))
check("region_share equals births_n / total births", {
  s <- region_raw %>% group_by(survey_year) %>% mutate(chk = births_n / sum(births_n))
  all(abs(s$chk - s$region_share) < 0.0005)
})

national_check <- region_raw %>%
  group_by(survey_year) %>%
  summarise(sba_calc = sum(sba_pct * region_share), stunt_calc = sum(stunted_pct * region_share), .groups = "drop") %>%
  inner_join(PUBLISHED_NATIONAL, by = "survey_year")
check("regional SBA averages back to the published national SBA (within 0.2 pp)",
      all(abs(national_check$sba_calc - national_check$sba_pub) < 0.2))
check("regional stunting averages back to the published national figure (within 0.5 pp)",
      all(abs(national_check$stunt_calc - national_check$stunt_pub) < 0.5))
# Tolerance is looser for stunting than for SBA (0.2pp) because region_share is a
# births-based weight (from Table 9.6), but stunting's true weight is children N
# (Table 11.1). The two weighting bases are close but not identical.

if (!all(results)) stop("Validation failed: fix the CSV(s) in data/raw/ before continuing.")
cat(sprintf("  -> %d/%d checks passed. Your transcription is internally consistent.\n", sum(results), length(results)))

# ---------------------------------------------------------------------------
# 2. LONG FORMAT: region / year / indicator / value / n
# ---------------------------------------------------------------------------
# Sample size for the stunting comparison: use the reported children N where we
# have it (2016 sub-region table), otherwise births N as a PROXY (2022).
region_wide <- region_raw %>%
  left_join(select(sub_rollup, survey_year, region, children_n), by = c("survey_year", "region")) %>%
  mutate(n_stunting       = children_n,          # real Table 11.1 children N, both years
         n_stunting_basis = "children N (reported)")

region_long <- region_wide %>%
  select(region, survey_year, region_share, sba_pct, stunted_pct, births_n, n_stunting, n_stunting_basis) %>%
  pivot_longer(c(sba_pct, stunted_pct), names_to = "indicator_code", values_to = "value") %>%
  mutate(indicator = recode(indicator_code,
                            sba_pct     = "Skilled birth attendance",
                            stunted_pct = "Stunting (under 5)"),
         n       = if_else(indicator_code == "sba_pct", births_n, n_stunting),
         n_basis = if_else(indicator_code == "sba_pct", "births N", n_stunting_basis))

# National rows = regional values weighted by each region's share of births.
national_long <- region_long %>%
  group_by(survey_year, indicator) %>%
  summarise(value = sum(value * region_share), n = sum(n), n_basis = first(n_basis), .groups = "drop") %>%
  mutate(region = "National")

region_long_all <- bind_rows(
  select(region_long, region, survey_year, indicator, value, n, n_basis),
  select(national_long, region, survey_year, indicator, value, n, n_basis)
) %>% rename(unit = region)

# ---------------------------------------------------------------------------
# 3. SUB-REGION CROSSWALK
# ---------------------------------------------------------------------------
# UDHS 2016 reports 15 sub-regions; UDHS 2022 reports 14. To compare like with
# like we merge/rename the 2016 units to the 2022 scheme:
#   South Central + North Central (2016) -> Buganda (2022)
#   Bugisu (2016)                        -> Elgon (2022)
# These two are APPROXIMATE: confirm the district lists in each report's
# regional-definition table/annex and note any mismatch in the README.
sub_h <- sub_raw %>%
  mutate(n_stunt       = coalesce(stunting_children_n, births_n),
         n_stunt_basis = if_else(is.na(stunting_children_n), "births N (proxy)", "children N (reported)"),
         subregion_h   = case_when(subregion %in% c("South Central", "North Central") ~ "Buganda",
                                   subregion == "Bugisu" ~ "Elgon",
                                   TRUE ~ subregion)) %>%
  group_by(survey_year, subregion_h, broad_region) %>%
  summarise(sba_pct       = weighted.mean(skilled_birth_attendance_pct, births_n),
            stunted_pct   = weighted.mean(stunted_pct, n_stunt),
            births_n      = sum(births_n),
            n_stunt       = sum(n_stunt),
            n_stunt_basis = first(n_stunt_basis),
            .groups = "drop") %>%
  mutate(crosswalk = if_else(subregion_h %in% c("Buganda", "Elgon"),
                             "approximate (merged/renamed)", "same name both years"))

stopifnot(n_distinct(sub_h$subregion_h) == 14)   # 14 comparable units in both years

sub_long <- sub_h %>%
  pivot_longer(c(sba_pct, stunted_pct), names_to = "indicator_code", values_to = "value") %>%
  mutate(indicator = recode(indicator_code,
                            sba_pct     = "Skilled birth attendance",
                            stunted_pct = "Stunting (under 5)"),
         n       = if_else(indicator_code == "sba_pct", births_n, n_stunt),
         n_basis = if_else(indicator_code == "sba_pct", "births N", n_stunt_basis)) %>%
  select(unit = subregion_h, broad_region, crosswalk, survey_year, indicator, value, n, n_basis)

# ---------------------------------------------------------------------------
# 4. IS THE 2016 -> 2022 CHANGE STATISTICALLY DETECTABLE?
# ---------------------------------------------------------------------------
# For each unit and indicator: change in percentage points, and a two-proportion
# z-test. se = sqrt( p1(1-p1)/n1 + p2(1-p2)/n2 ).
# We report it twice:
#   naive  : treats the survey as a simple random sample (design effect = 1)
#   deff=2 : inflates the variance x2, a conservative stand-in for DHS clustering
V_ROBUST <- "Significant (robust to design effect = 2)"
V_NAIVE  <- "Significant only if design effect ignored"
V_NONE   <- "No statistically detectable change"

compare_two_years <- function(long_df) {
  id_cols <- intersect(c("unit", "broad_region", "crosswalk", "indicator"), names(long_df))
  long_df %>%
    pivot_wider(id_cols = all_of(id_cols), names_from = survey_year,
                values_from = c(value, n, n_basis), names_sep = "_") %>%
    mutate(
      change_pp      = value_2022 - value_2016,
      improvement_pp = if_else(indicator == "Stunting (under 5)", -change_pp, change_pp),
      p1 = value_2016 / 100,
      p2 = value_2022 / 100,
      se_pp          = 100 * sqrt(p1 * (1 - p1) / n_2016 + p2 * (1 - p2) / n_2022),
      ci_low         = change_pp - 1.96 * se_pp,
      ci_high        = change_pp + 1.96 * se_pp,
      p_naive        = 2 * pnorm(-abs(change_pp / se_pp)),
      se_pp_deff2    = sqrt(2) * se_pp,
      ci_low_deff2   = change_pp - 1.96 * se_pp_deff2,
      ci_high_deff2  = change_pp + 1.96 * se_pp_deff2,
      p_deff2        = 2 * pnorm(-abs(change_pp / se_pp_deff2)),
      verdict        = case_when(p_naive < 0.05 & p_deff2 < 0.05 ~ V_ROBUST,
                                 p_naive < 0.05                  ~ V_NAIVE,
                                 TRUE                            ~ V_NONE),
      n_caveat       = ""   
    ) %>%
    select(-p1, -p2)
}

region_change <- compare_two_years(region_long_all) %>%
  group_by(indicator) %>%
  mutate(improvement_rank = as.integer(rank(-if_else(unit == "National", NA_real_, improvement_pp),
                                            na.last = "keep"))) %>%
  ungroup() %>%
  arrange(indicator, improvement_rank)

sub_change <- compare_two_years(sub_long) %>%
  group_by(indicator) %>%
  mutate(p_holm_deff2    = p.adjust(p_deff2, method = "holm"),   # 14 tests per indicator
         improvement_rank = as.integer(rank(-improvement_pp))) %>%
  ungroup() %>%
  arrange(indicator, improvement_rank)

cat("\n=== 4a. Broad regions: change 2016 -> 2022 ===\n")
print(region_change %>%
        transmute(indicator, unit, from = round(value_2016, 1), to = round(value_2022, 1),
                  change_pp = round(change_pp, 1),
                  ci_naive = sprintf("[%.1f, %.1f]", ci_low, ci_high),
                  ci_deff2 = sprintf("[%.1f, %.1f]", ci_low_deff2, ci_high_deff2),
                  p_naive = round(p_naive, 4), p_deff2 = round(p_deff2, 4), verdict),
      n = Inf, width = Inf)

cat("\n=== 4b. Sub-regions: change 2016 -> 2022 (14 comparable units) ===\n")
print(sub_change %>%
        transmute(indicator, unit, from = round(value_2016, 1), to = round(value_2022, 1),
                  change_pp = round(change_pp, 1), p_naive = round(p_naive, 4),
                  p_deff2 = round(p_deff2, 4), p_holm_deff2 = round(p_holm_deff2, 4), verdict),
      n = Inf, width = Inf)
cat("p_holm_deff2 = design-effect p-value adjusted for running 14 tests per indicator (Holm).\n",
    "It is the strictest column: a change that survives it is very unlikely to be a fluke.\n", sep = "")

# ---------------------------------------------------------------------------
# 5. EXPLORATORY FIGURES  (titles are computed from the results, so they can't
#    drift out of sync with the numbers)
# ---------------------------------------------------------------------------
wrap <- function(x, width = 70) paste(strwrap(x, width), collapse = "\n")
SRC  <- "Source: UDHS 2016 and 2022 published tables (hand-extracted)."
theme_set(theme_minimal(base_size = 12) +
            theme(plot.title.position = "plot", legend.position = "bottom",
                  plot.caption = element_text(size = 8, colour = "grey40")))

nat <- region_change %>% filter(unit == "National")
nat_sba   <- nat %>% filter(indicator == "Skilled birth attendance")
nat_stunt <- nat %>% filter(indicator == "Stunting (under 5)")

# 5a. Trend lines per indicator ---------------------------------------------
p_trend <- region_long_all %>%
  mutate(is_national = unit == "National") %>%
  ggplot(aes(x = factor(survey_year), y = value, group = unit, colour = unit)) +
  geom_line(aes(linetype = is_national), linewidth = 1) +
  geom_point(size = 2.5) +
  scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"), guide = "none") +
  scale_colour_manual(values = c(Central = "#1b9e77", Eastern = "#d95f02", Northern = "#7570b3",
                                 Western = "#e7298a", National = "black")) +
  facet_wrap(~ indicator, scales = "free_y") +
  labs(title = wrap(sprintf(
         "Skilled birth attendance rose %.1f points nationally (%.1f%% to %.1f%%); stunting fell %.1f points (%.1f%% to %.1f%%)",
         nat_sba$change_pp, nat_sba$value_2016, nat_sba$value_2022,
         abs(nat_stunt$change_pp), nat_stunt$value_2016, nat_stunt$value_2022), 90),
       x = "UDHS survey year", y = "Percent", colour = NULL, caption = SRC)
ggsave(file.path(FIG, "01_trend_by_indicator.png"), p_trend, width = 10, height = 5.5, dpi = 150)

# 5b. Regional comparison bars (bars start at zero) --------------------------
r22 <- region_change %>% filter(unit != "National")
hi_stunt <- r22 %>% filter(indicator == "Stunting (under 5)") %>% slice_max(value_2022, n = 1)
hi_sba   <- r22 %>% filter(indicator == "Skilled birth attendance") %>% slice_max(value_2022, n = 1)
lo_sba   <- r22 %>% filter(indicator == "Skilled birth attendance") %>% slice_min(value_2022, n = 1)

p_bars <- region_long_all %>%
  filter(unit != "National") %>%
  ggplot(aes(x = unit, y = value, fill = factor(survey_year))) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(width = 0.85),
            vjust = -0.4, size = 3) +
  facet_wrap(~ indicator, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c(`2016` = "#9ecae1", `2022` = "#08519c")) +
  labs(title = wrap(sprintf(
         "In 2022, %s had the highest stunting (%.1f%%) and %s the lowest skilled birth attendance (%.1f%%)",
         hi_stunt$unit, hi_stunt$value_2022, lo_sba$unit, lo_sba$value_2022), 90),
       x = NULL, y = "Percent", fill = "Survey year", caption = SRC)
ggsave(file.path(FIG, "02_regional_comparison_bars.png"), p_bars, width = 10, height = 5.5, dpi = 150)

# 5c. Sub-region dumbbell plots ----------------------------------------------
verdict_cols <- setNames(c("#1b7837", "#f1a340", "grey70"), c(V_ROBUST, V_NAIVE, V_NONE))

dumbbell <- function(ind, title) {
  d <- sub_change %>% filter(indicator == ind) %>%
    mutate(unit_lab = if_else(crosswalk == "same name both years", unit, paste0(unit, "*")),
           unit_lab = reorder(unit_lab, value_2022 - value_2016))
  pts <- d %>% select(unit_lab, value_2016, value_2022) %>%
    pivot_longer(-unit_lab, names_to = "year", values_to = "value") %>%
    mutate(year = sub("value_", "", year))
  ggplot() +
    geom_segment(data = d, aes(x = value_2016, xend = value_2022, y = unit_lab, yend = unit_lab, colour = verdict),
                 linewidth = 1.3) +
    geom_point(data = pts, aes(x = value, y = unit_lab, fill = year), shape = 21, size = 3, stroke = 1) +
    scale_fill_manual(values = c(`2016` = "white", `2022` = "black")) +
    scale_colour_manual(values = verdict_cols, drop = FALSE) +
    labs(title = wrap(title, 72), x = "Percent", y = NULL, fill = "Survey year",
         colour = "2016 to 2022 change",
         caption = paste(SRC, "* Merged/renamed for comparability: South + North Central = Buganda; Bugisu = Elgon.",
                         sep = "\n")) +
    guides(colour = guide_legend(nrow = 3), fill = guide_legend(nrow = 2, override.aes = list(shape = 21)))
}

s_sba   <- sub_change %>% filter(indicator == "Skilled birth attendance")
s_stunt <- sub_change %>% filter(indicator == "Stunting (under 5)")

p_sub_sba <- dumbbell("Skilled birth attendance", sprintf(
  "Skilled birth attendance rose in %d of 14 sub-regions; %d of those changes hold up under a design-effect check",
  sum(s_sba$change_pp > 0), sum(s_sba$change_pp > 0 & s_sba$verdict == V_ROBUST)))
ggsave(file.path(FIG, "03_subregion_sba_change.png"), p_sub_sba, width = 9, height = 7, dpi = 150)

p_sub_stunt <- dumbbell("Stunting (under 5)", sprintf(
  "Stunting fell in %d of 14 sub-regions and rose in %d; only %d changes hold up under a design-effect check",
  sum(s_stunt$change_pp < 0), sum(s_stunt$change_pp > 0), sum(s_stunt$verdict == V_ROBUST)))
ggsave(file.path(FIG, "04_subregion_stunting_change.png"), p_sub_stunt, width = 9, height = 7, dpi = 150)

# ---------------------------------------------------------------------------
# 6. Do SBA and stunting move together across sub-regions?
# ---------------------------------------------------------------------------
# CAUTION: this is an ECOLOGICAL analysis (14 areas, not people). A correlation
# between areas does not show that individual mothers who deliver with a skilled
# attendant have less-stunted children. Spearman is used because n is tiny.
sub_wide <- sub_h %>%
  select(survey_year, subregion_h, sba_pct, stunted_pct) %>%
  pivot_wider(names_from = survey_year, values_from = c(sba_pct, stunted_pct))

run_cor <- function(x, y, label) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(comparison = label, n_subregions = length(x),
         spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}
correlations <- bind_rows(
  run_cor(sub_wide$sba_pct_2016, sub_wide$stunted_pct_2016, "2016 levels: SBA vs stunting"),
  run_cor(sub_wide$sba_pct_2022, sub_wide$stunted_pct_2022, "2022 levels: SBA vs stunting"),
  run_cor(sub_wide$sba_pct_2022 - sub_wide$sba_pct_2016,
          sub_wide$stunted_pct_2022 - sub_wide$stunted_pct_2016,
          "Change 2016 to 2022: SBA change vs stunting change")
)
cat("\n=== 6. Ecological correlations across the 14 sub-regions ===\n")
print(correlations %>% mutate(spearman_rho = round(spearman_rho, 2), p_value = round(p_value, 3)),
      n = Inf, width = Inf)

c22 <- correlations[2, ]
p_scatter <- ggplot(sub_wide, aes(x = sba_pct_2022, y = stunted_pct_2022)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey60", linewidth = 0.6, formula = y ~ x) +
  geom_point(size = 2.5, colour = "#08519c") +
  geom_text(aes(label = subregion_h), vjust = -0.8, size = 3, check_overlap = TRUE) +
  labs(title = wrap(sprintf(
         "Across 14 sub-regions in 2022, skilled birth attendance and stunting are %s correlated (Spearman rho = %.2f, p = %.2f)",
         if (c22$p_value < 0.05) "significantly" else "not significantly", c22$spearman_rho, c22$p_value), 90),
       x = "Skilled birth attendance, 2022 (%)", y = "Stunting under 5, 2022 (%)",
       caption = paste(SRC, "Area-level (ecological) association, not an individual-level effect.", sep = "\n"))
ggsave(file.path(FIG, "05_subregion_sba_vs_stunting_2022.png"), p_scatter, width = 8, height = 6, dpi = 150)

# ---------------------------------------------------------------------------
# 7. SAVE TIDY OUTPUTS (for the R Markdown report and Power BI)
# ---------------------------------------------------------------------------
write_csv(region_long_all, file.path(PROC, "region_long.csv"))
write_csv(sub_long,        file.path(PROC, "subregion_long_harmonized.csv"))
write_csv(region_change,   file.path(PROC, "region_change_summary.csv"))
write_csv(sub_change,      file.path(PROC, "subregion_change_summary.csv"))
write_csv(correlations,    file.path(PROC, "subregion_correlations.csv"))

cat("\n=== Done ===\n")
cat("Tables written to ", normalizePath(PROC), "\n", sep = "")
cat("Figures written to ", normalizePath(FIG), "\n", sep = "")

# ---------------------------------------------------------------------------
# 8. NATIONAL-ONLY INDICATORS: mortality, MMR, fertility rate
# ---------------------------------------------------------------------------
# These three do NOT belong in the regional analysis above:
#   - Neonatal/infant/under-5 mortality: DHS's regional breakdown (Table 8.2)
#     uses a 10-year reference period, while the national headline figure uses
#     a 5-year period -- they are not on the same time basis, so no attempt is
#     made to reconcile them here. National point estimates only.
#   - Maternal mortality ratio: no regional breakdown exists in any UDHS round
#     (the sisterhood estimation method doesn't support sub-national estimates).
#   - Total fertility rate: a summary rate across a woman's reproductive
#     lifespan, not a per-birth outcome -- national only, same as MMR.
#
# Significance testing: MMR has a REAL published 95% CI in the report, so its
# change is tested properly by converting each CI to an implied standard error
# (se = (upper - lower) / (2 * 1.96)) and running a two-sample z-test on the
# difference. Mortality and TFR have no published standard error in the
# tables extracted for this project 

natl_path <- file.path(RAW, "udhs_national_indicators.csv")
if (!file.exists(natl_path)) {
  stop("Can't find ", natl_path, ". Run from the r/ folder, with the file in data/raw/.")
}
natl_raw <- read_csv(natl_path, show_col_types = FALSE)

cat("\n=== 8a. Validating the national indicators file ===\n")
check("national file has the expected columns",
      all(c("survey_year", "indicator", "value", "ci_lower", "ci_upper", "se_basis") %in% names(natl_raw)))
check("exactly 5 indicators x 2 years = 10 rows", nrow(natl_raw) == 10)
check("MMR rows have a reported CI; mortality/TFR rows do not",
      all(!is.na(natl_raw$ci_lower[natl_raw$indicator == "Maternal mortality ratio"])) &&
      all(is.na(natl_raw$ci_lower[natl_raw$indicator != "Maternal mortality ratio"])))

natl_wide <- natl_raw %>%
  select(survey_year, indicator, value, ci_lower, ci_upper) %>%
  pivot_wider(names_from = survey_year, values_from = c(value, ci_lower, ci_upper), names_sep = "_")

# MMR: proper CI-derived z-test.
mmr <- natl_wide %>% filter(indicator == "Maternal mortality ratio") %>%
  mutate(
    change      = value_2022 - value_2016,
    se_2016     = (ci_upper_2016 - ci_lower_2016) / (2 * 1.96),
    se_2022     = (ci_upper_2022 - ci_lower_2022) / (2 * 1.96),
    se_diff     = sqrt(se_2016^2 + se_2022^2),
    z           = change / se_diff,
    p_value     = 2 * pnorm(-abs(z)),
    verdict     = if_else(p_value < 0.05, "Significant (based on published CIs)", "No statistically detectable change")
  )

# Mortality + TFR: point-estimate change only, no test (see header note).
other_natl <- natl_wide %>%
  filter(indicator != "Maternal mortality ratio") %>%
  mutate(change = value_2022 - value_2016,
         verdict = "No formal test: no published standard error available")

cat("\n=== 8b. Maternal mortality ratio: change 2016 -> 2022 (tested using published CIs) ===\n")
print(mmr %>% transmute(indicator, from = value_2016, to = value_2022, change = round(change, 0),
                        se_diff = round(se_diff, 1), z = round(z, 2), p_value = round(p_value, 4), verdict),
      width = Inf)

cat("\n=== 8c. Mortality and fertility: point-estimate change only (no formal test) ===\n")
print(other_natl %>% transmute(indicator, from = value_2016, to = value_2022, change = round(change, 1), verdict),
      width = Inf)

# ---- 8d. Trend figure ----
natl_change_all <- bind_rows(
  mmr %>% select(indicator, value_2016, value_2022, change),
  other_natl %>% select(indicator, value_2016, value_2022, change)
)
u5 <- natl_change_all %>% filter(indicator == "Under-5 mortality")
mmr_row <- natl_change_all %>% filter(indicator == "Maternal mortality ratio")

mmr_sig_phrase <- if (mmr$p_value < 0.05) "a statistically significant decline" else "a change not statistically distinguishable from no change"

p_national <- natl_raw %>%
  ggplot(aes(x = factor(survey_year), y = value, group = 1)) +
  geom_line(colour = "#08519c", linewidth = 1) +
  geom_point(size = 2.5, colour = "#08519c") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, na.rm = TRUE) +
  facet_wrap(~ indicator, scales = "free_y") +
  labs(title = wrap(sprintf(
         "Under-5 mortality fell from %.0f to %.0f per 1,000 live births; maternal mortality ratio fell from %.0f to %.0f per 100,000, %s (p = %.4f)",
         u5$value_2016, u5$value_2022, mmr_row$value_2016, mmr_row$value_2022,
         mmr_sig_phrase, mmr$p_value), 95),
       x = "UDHS survey year", y = NULL,
       caption = paste(SRC, "National point estimates. Error bars shown only where a published CI exists (MMR).", sep = "\n"))
ggsave(file.path(FIG, "06_national_mortality_mmr_tfr.png"), p_national, width = 10, height = 6, dpi = 150)

write_csv(mmr, file.path(PROC, "national_mmr_test.csv"))
write_csv(other_natl, file.path(PROC, "national_mortality_tfr_change.csv"))
cat("\nNational indicator tables and figure written alongside the regional outputs.\n")
