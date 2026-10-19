# Uganda Maternal & Child Health: Progress 2016–2022
### Statistical Analysis of UDHS Indicators in R, with a Power BI Dashboard

**Fiona Nakabugo Kabuka** — Healthcare Data Analyst | M&E | Power BI | SQL | Data Management

A reproducible analysis of how maternal and child health changed in Uganda between the 2016 and 2022 Demographic and Health Surveys. It tests whether each reported change is statistically real, region by region, and it shows a regression workflow on a clearly labelled synthetic cohort. The real findings are kept strictly apart from the method demonstration.

---

## Overview

**Problem:** Survey reports give two numbers per indicator, such as 74% in 2016 and 88% in 2022. A programme manager needs to know whether a change like that is real, where it happened, and how far the comparison can be trusted.

**Process:** R. Published UDHS tables were extracted by hand and validated against the published national figures. Each 2016→2022 change was tested with a two-proportion z-test. The tests were adjusted for the survey's clustered design, and corrected for multiple comparisons across 14 sub-regions. A regression workflow was then built on a synthetic cohort calibrated to the real regional rates. Results are reported in an R Markdown document and a one-page Power BI dashboard.

**Outcome:** An HTML report in which every number is computed from the analysis files, and a dashboard that matches it. Both state what the data can and cannot support.

---

## Report and Dashboard Preview
https://nkgfiona.github.io/uganda-mch-udhs-analysis/r/mch_report.html
<img width="1726" height="935" alt="SBA" src="https://github.com/user-attachments/assets/5d6af2f1-02cc-4630-95a6-04b2db001cea" />
<img width="1720" height="924" alt="Stunting" src="https://github.com/user-attachments/assets/e78cd54f-36e3-4958-b175-96ca1bbc32d1" />

- **Report:** `r/mch_report.html` is the rendered report. Open it in any browser.
- **Dashboard:** `powerbi/mch.pbix` opens in [Power BI Desktop](https://www.microsoft.com/en-us/power-platform/products/power-bi/downloads) (free, Windows only). It contains only published aggregate figures. No individual-level data is embedded. See **The Dashboard** below.

---

## Tech Stack

| Layer | Tools |
|---|---|
| Data preparation and testing | R (dplyr, tidyr, readr) |
| Figures | ggplot2 |
| Regression | base R `glm()` (logistic) |
| Reporting | R Markdown (knitr, rmarkdown) |
| Dashboard | Power BI (DAX measures, indicator buttons, conditional colour) |

---

## Analysis Pipeline

```
UDHS 2016 and 2022 published tables (hand-extracted)
        │
        ▼
data/raw/udhs_*.csv   ── validated against published national figures
        │
        ├──► mch_aggregate_analysis.R ──► data/processed/*.csv + figures 01–05
        │        (tests, design effect, Holm correction, sub-region matching)
        │
        └──► generate_synthetic_mch_data.R ──► data/raw/mch_synthetic_data.csv
                                                   │
                                                   ▼
                                    mch_regression_analysis.R ──► data/processed/*.csv
                                                   │
                                                   ▼
        mch_report.Rmd (reads the saved files only) ──► mch_report.html
        powerbi/mch.pbix (reads the real aggregate tables only)
```

The report never re-runs the analysis. It reads the saved outputs and stops with a clear message if any are missing, so every figure in it traces back to a saved file.

---

## Data Sources & Disclosure

**Real and cited:** all findings in the National Trends and Regional Findings sections come from published UDHS tables, transcribed by hand and checked for internal consistency.

- Skilled birth attendance (SBA): Table 9.6 in both reports.
- Stunting: Table 11.1 in both reports (children aged 0–59 months, so both years use the same age base).
- Maternal mortality ratio and total fertility rate: Table 15.4 (2016) and Table 16.4 (2022).
- National neonatal, infant and under-five mortality: the published national figures.

**Synthetic and disclosed:** the regression section uses a cohort of 21,000 synthetic records.

- Each region and year was calibrated to its real published SBA and stunting rate, and the real births counts set the cohort size.
- The effects of wealth, maternal education, residence and antenatal care were **assumed**, adapted from a published analysis of Nigeria's 2018 DHS (odds ratios around 6.3 for wealth and 7.0 for education). They are not Uganda estimates.
- DHS's Terms of Use say its micro-level data "will not be re-distributed, either directly or within any tool/dashboard." This project used no DHS microdata. That is why the individual-level part is synthetic.

So the regression shows that the modelling code works. It does not estimate what drives these outcomes in Uganda. A recovery check supports this: 17 of 17 planted effects fell inside their fitted confidence intervals.

---

## Project Structure

```
02_m_and_e_dashboard/
├── data/
│   ├── raw/                    # UDHS reference tables, plus the generated synthetic cohort
│   └── processed/              # saved results read by the report and dashboard
├── r/
│   ├── mch_aggregate_analysis.R
│   ├── generate_synthetic_mch_data.R
│   ├── mch_regression_analysis.R
│   ├── mch_report.Rmd
│   └── mch_report.html
├── powerbi/
│   └── mch.pbix
├── documentation/
│   └── figures/                # figures 01–06
└── README.md
```

---

## How to Run It Yourself

In R, install the packages once:

```r
install.packages(c("dplyr", "tidyr", "readr", "ggplot2", "knitr", "rmarkdown"))
```

Then set the working directory to `r/` and run the scripts **in this order**:

```r
source("mch_aggregate_analysis.R")       # real-data tests, figures, saved results
source("generate_synthetic_mch_data.R")  # synthetic cohort (fixed seed, repeatable)
source("mch_regression_analysis.R")      # regression on the synthetic cohort
rmarkdown::render("mch_report.Rmd")      # builds mch_report.html
```

The generator uses a fixed seed, so the same data and results come out every time.

---

## The Dashboard

A single page, **UDHS 2016 vs 2022**, built only from the real published aggregates.

- **Indicator buttons** in the header switch the page between *Skilled birth attendance* and *Child stunting (under 5)*. The page opens on stunting, so it never shows the two indicators mixed together.
- **National cards** show each indicator's 2016 and 2022 values, the change, and a verdict. They do not respond to the buttons, so both national results stay visible.
- **Regional chart** compares 2016 and 2022 for the four regions. Its title and subtitle follow the selected indicator and state the direction of improvement ("higher is better" for skilled birth attendance, "lower is better" for stunting).
- **Sub-region chart** shows each sub-region's change in percentage points, coloured by verdict: green where the change is robust to the design effect, amber where it is significant only if the design effect is ignored, and grey where no change is detectable. Its subtitle counts how many of the 14 changes hold up under the design-effect check and how many survive Holm correction.
- **Sub-region table** shows both indicators side by side, one row per sub-region: 2016 and 2022 values, change, confidence interval, Holm-adjusted p-value and a colour-coded verdict. It does not depend on the buttons.
- **Notes on the page** give the legend, the source tables, the approximate Elgon and Buganda matches, and the different SBA reference periods (2016: births in the last 5 years; 2022: births in the last 2 years).

---

## Key Findings (real UDHS data)

- **Skilled birth attendance rose from 74.1% to 88.4%** (+14.3 points, p < 0.0001). It rose in all four regions and in all 14 sub-regions. The Western Region gained the most (+18.4 points). After Holm correction, 12 of the 14 sub-region gains remain statistically significant.
- **Child stunting fell from 28.8% to 24.0%** (−4.8 points, p = 0.0001). The national fall is robust, but the regional and sub-regional falls are weaker. Only two of the four regions hold up under the design-effect check, and no sub-region fall survives Holm correction. Stunting fell in 11 sub-regions and rose in 3.
- **Maternal mortality fell from 336 to 189 per 100,000 live births** (p = 0.0006, tested using the confidence intervals published in each report).
- **The regional gap in stunting widened** from 9.0 to 10.5 points, with the Western Region the highest, even as the national level fell.

---

## Limitations

- **The SBA comparison is not like-for-like in time.** The 2016 table covers births in the 5 years before the survey, and the 2022 table covers the 2 years before it. The two windows are roughly 7.5 years apart, not 6, so part of the measured gain reflects the different windows.
- **The design effect is assumed.** True design effects are not in the published summary tables, so tests are reported both naively and with a conservative design effect of 2.
- **Sub-region matching is approximate in two places.** The 2016 sub-regions South Central and North Central are combined into Buganda, and Bugisu is treated as Elgon.
- **Mortality and fertility are national only.** Regional mortality in the reports uses a 10-year window while the national figure uses 5 years, so the two cannot be reconciled. Maternal mortality has no regional estimate in any UDHS round. No standard errors are published for mortality and fertility, so only their point changes are reported.
- **Ecological patterns are not individual-level relationships.** Any association across sub-regions describes areas, not people.
- **The regression is a demonstration.** See Data Sources & Disclosure.

---

## Skills Demonstrated

- **M&E analysis:** turning published survey indicators into tested, interpretable change estimates, region by region.
- **Statistical judgement:** design-effect-adjusted tests, Holm correction across 14 simultaneous comparisons, a confidence-interval-based test where no sample size exists, and a clear refusal to run a test the data cannot support.
- **Data integrity:** hand-extracted tables validated against published totals, and a reference-window mismatch caught and disclosed.
- **Reproducible reporting:** an R Markdown report that computes every number from saved files.
- **Honest use of synthetic data:** real findings and a method demonstration kept visibly separate, with the legal reason stated.
- **Power BI:** DAX measures driven by indicator buttons, titles that change with the selection, verdict-based colour formatting, and a side-by-side table that stays independent of the buttons.

---

## About

**Fiona Nakabugo Kabuka**
Healthcare Data Analyst | M&E | Power BI | SQL | Data Management

[LinkedIn — www.linkedin.com/in/nakabugo-fiona-8334b638] · [Email — fnakabugo@gmail.com]

---

## Data Source Citations

- Uganda Bureau of Statistics (UBOS) and ICF. *Uganda Demographic and Health Survey 2016.* [FR333](https://dhsprogram.com/pubs/pdf/FR333/FR333.pdf)
- Uganda Bureau of Statistics (UBOS) and ICF. *Uganda Demographic and Health Survey 2022.* [Final report](https://www.ubos.org/wp-content/uploads/publications/UDHS-2022-Report.pdf)
- The DHS Program. [Datasets Terms of Use](https://dhsprogram.com/data/terms-of-use.cfm)
- Effect sizes for wealth and education were adapted from a published analysis of Nigeria's 2018 DHS: [Determinants of Skilled Birth Attendance in Nigeria](https://www.medrxiv.org/content/10.64898/2026.04.23.26350432.full.pdf) (preprint).
