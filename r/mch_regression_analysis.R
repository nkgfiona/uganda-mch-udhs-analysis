# mch_regression_analysis.R
# ---------------------------------------------------------------------------
# Logistic regression of skilled birth attendance (SBA) and child stunting,
# using the synthetic dataset from generate_synthetic_mch_data.R.
#
# WHAT THE SCRIPT DOES
#   This is a METHOD DEMONSTRATION, not a set of new findings about Uganda.
#     * Wealth, education, urban and ANC4+ effects were BUILT IN from a
#       published Nigerian DHS analysis. The regression is expected to recover
#       them (see the "recovery check" below). That shows the code and method
#       work; it does not tell you what those effects are in Uganda.
#     * Region and survey-year differences are calibrated to REAL UDHS regional
#       rates, so they mirror real published gaps. They are inputs that have been
#       passed through a model, not independent discoveries.
#     * The synthetic data has no survey weights, clusters or strata. With real
#       DHS microdata you would use survey::svydesign() + survey::svyglm().
#   The REAL findings of this project come from mch_aggregate_analysis.R.
#
# Run from the r/ folder, after generate_synthetic_mch_data.R.
# Writes tidy tables to ../data/processed/
# ---------------------------------------------------------------------------

DATA_PATH <- "../data/raw/mch_synthetic_data.csv"
PROC      <- "../data/processed"
if (!file.exists(DATA_PATH)) stop("Run generate_synthetic_mch_data.R first (from the r/ folder).")
dir.create(PROC, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)

# ---- Reference categories ---------------------------------------------------
# Wealth, education and residence use the MOST DISADVANTAGED group as reference,
# so odds ratios read as "X times the odds versus the poorest / least educated /
# rural". Region has no single "most disadvantaged" group across both outcomes,
# so Central (the largest, best-served region) is the reference: region odds
# ratios read as "compared with Central".
df$wealth_quintile   <- factor(df$wealth_quintile, levels = c("Poorest", "Poorer", "Middle", "Richer", "Richest"))
df$mother_education  <- factor(df$mother_education, levels = c("None", "Primary", "Secondary", "Higher"))
df$region            <- factor(df$region, levels = c("Central", "Eastern", "Northern", "Western"))
df$residence         <- factor(df$residence, levels = c("Rural", "Urban"))
df$anc4_visits       <- factor(df$anc4_visits, levels = c("No", "Yes"))
df$survey_year_f     <- factor(df$survey_year, levels = c(2016, 2022))

df$sba_binary   <- ifelse(df$skilled_birth_attendance == "Yes", 1, 0)
df$stunt_binary <- ifelse(df$stunted == "Yes", 1, 0)

# ---- Helpers ----------------------------------------------------------------
or_table <- function(model) {
  co <- summary(model)$coefficients
  or <- exp(co[, "Estimate"])
  ci <- exp(confint.default(model))
  data.frame(
    term       = rownames(co),
    odds_ratio = round(or, 2),
    ci_lower   = round(ci[, 1], 2),
    ci_upper   = round(ci[, 2], 2),
    p_value    = round(co[, "Pr(>|z|)"], 4),
    row.names  = NULL
  )
}

# Area under the ROC curve via the rank (Mann-Whitney) formula: base R only.
auc <- function(y, p) {
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_stats <- function(model, y, label) {
  p <- fitted(model)
  data.frame(
    model             = label,
    n                 = length(y),
    mcfadden_r2       = round(1 - model$deviance / model$null.deviance, 3),
    auc               = round(auc(y, p), 3),
    accuracy          = round(mean((p > 0.5) == y), 3),
    baseline_accuracy = round(max(mean(y), 1 - mean(y)), 3)   # always guessing the commoner outcome
  )
}

# ---- MODEL 1: Skilled birth attendance ---------------------------------------
cat("=====================================================\n")
cat(" MODEL 1: Skilled Birth Attendance (odds ratios)\n")
cat("=====================================================\n")
model_sba <- glm(
  sba_binary ~ region + residence + wealth_quintile + mother_education + anc4_visits + survey_year_f,
  data = df, family = binomial()
)
or_sba <- or_table(model_sba)
print(or_sba)
cat("\nNote: anc4_visits sits partly ON the pathway between wealth/education and SBA.\n",
    "Its presence shrinks the wealth/education odds ratios (mediation): interpret as\n",
    "'holding ANC4+ constant', not as total effects.\n", sep = "")

# ---- MODEL 2: Stunting -------------------------------------------------------
cat("\n=====================================================\n")
cat(" MODEL 2: Child Stunting (odds ratios)\n")
cat("=====================================================\n")
model_stunt <- glm(
  stunt_binary ~ region + residence + wealth_quintile + mother_education + survey_year_f,
  data = df, family = binomial()
)
or_stunt <- or_table(model_stunt)
print(or_stunt)

# ---- Did regional gaps change between 2016 and 2022? -------------------------
# Likelihood-ratio test: does adding region x year interaction terms improve fit?
cat("\n=====================================================\n")
cat(" REGION x YEAR: did regional gaps change 2016 -> 2022?\n")
cat("=====================================================\n")
lrt_row <- function(m_main, label) {
  m_int <- update(m_main, . ~ . + region:survey_year_f)
  a <- anova(m_main, m_int, test = "Chisq")
  data.frame(outcome = label, df = a$Df[2], deviance_drop = round(a$Deviance[2], 1),
             p_value = signif(a$`Pr(>Chi)`[2], 3))
}
interaction_tests <- rbind(lrt_row(model_sba, "Skilled birth attendance"),
                           lrt_row(model_stunt, "Stunting"))
print(interaction_tests, row.names = FALSE)
cat("Reading: a small p-value means regions did NOT all change by the same amount (on the odds scale).\n",
    "Here that reflects the real regional rates the data was calibrated to, not a new discovery.\n", sep = "")

# ---- Model fit ---------------------------------------------------------------
cat("\n=====================================================\n")
cat(" MODEL FIT SUMMARY\n")
cat("=====================================================\n")
fit <- rbind(fit_stats(model_sba, df$sba_binary, "SBA"),
             fit_stats(model_stunt, df$stunt_binary, "Stunting"))
print(fit, row.names = FALSE)
cat("Accuracy only means something next to baseline_accuracy (always guessing the commoner outcome).\n",
    "AUC (0.5 = coin flip, 1 = perfect) is the fairer summary of how well each model separates outcomes.\n", sep = "")

# ---- Recovery check: did the regression find what was planted? --------------
# These are the log-odds effects planted in generate_synthetic_mch_data.R.
# For stunting the generator multiplies the SBA-style effects by -0.6 (wealth),
# -0.5 (education) and -0.5 (urban, whose base effect is 0.6, so -0.3 planted).
# KEEP IN SYNC with the generator if you change effect sizes there.
planted <- data.frame(
  term = c("wealth_quintilePoorer", "wealth_quintileMiddle", "wealth_quintileRicher", "wealth_quintileRichest",
           "mother_educationPrimary", "mother_educationSecondary", "mother_educationHigher",
           "residenceUrban", "anc4_visitsYes"),
  sba_logodds     = c(0.4, 0.9, 1.4, 1.84, 0.5, 1.1, 1.95, 0.6, 1.33),
  stunt_logodds   = c(-0.24, -0.54, -0.84, -1.104, -0.25, -0.55, -0.975, -0.3, NA)
)
recovery <- function(or_tab, planted_col) {
  m <- merge(planted[, c("term", planted_col)], or_tab, by = "term", sort = FALSE)
  m <- m[!is.na(m[[planted_col]]), ]
  data.frame(term = m$term,
             planted_OR = round(exp(m[[planted_col]]), 2),
             estimated_OR = m$odds_ratio, ci_lower = m$ci_lower, ci_upper = m$ci_upper,
             planted_inside_CI = exp(m[[planted_col]]) >= m$ci_lower & exp(m[[planted_col]]) <= m$ci_upper,
             row.names = NULL)
}
rec_sba   <- recovery(or_sba, "sba_logodds")
rec_stunt <- recovery(or_stunt, "stunt_logodds")
cat("\n=====================================================\n")
cat(" RECOVERY CHECK: planted vs estimated odds ratios\n")
cat("=====================================================\n")
cat("SBA model:\n");      print(rec_sba, row.names = FALSE)
cat("\nStunting model:\n"); print(rec_stunt, row.names = FALSE)
cat(sprintf("\nPlanted effect inside the 95%% CI: SBA %d/%d, stunting %d/%d (about 19 of 20 expected by chance).\n",
            sum(rec_sba$planted_inside_CI), nrow(rec_sba), sum(rec_stunt$planted_inside_CI), nrow(rec_stunt)))

# ---- Save tidy outputs -------------------------------------------------------
write.csv(or_sba,             file.path(PROC, "or_table_sba.csv"),             row.names = FALSE)
write.csv(or_stunt,           file.path(PROC, "or_table_stunting.csv"),        row.names = FALSE)
write.csv(fit,                file.path(PROC, "model_fit_summary.csv"),        row.names = FALSE)
write.csv(interaction_tests,  file.path(PROC, "region_year_interaction_tests.csv"), row.names = FALSE)
write.csv(rbind(cbind(model = "SBA", rec_sba), cbind(model = "Stunting", rec_stunt)),
          file.path(PROC, "recovery_check.csv"), row.names = FALSE)
cat("\nTables written to", normalizePath(PROC), "\n")
