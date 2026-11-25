rm(list=ls())
options(scipen=999)
select <- dplyr::select
vars <- dplyr::vars
coalesce <- dplyr::coalesce

getwd()

library(tidyverse)
library(haven)
library(survey)
library(brms)
library(cmdstanr)
library(bayesplot)
library(expss)
library(lavaan)
library(weights)
library(semTools)
library(modelsummary)
library(DescTools)
library(scales)
library(ggh4x)
library(openxlsx)
library(kableExtra)

data.0 <- foreign::read.spss(paste0("C:/Users/timgr/OneDrive - Wilfrid Laurier University/Academic/UK Nuclear Energy-NIA/",
  "NIA - Nuclear Opinions Survey - SPSS 2026 UK.sav"),
   to.data.frame = TRUE, use.value.labels = TRUE)

summary(data.0$weight)

data.1 <- data.0 %>%
  select(-contains("OT", ignore.case = FALSE)) %>%
  mutate(support_nuclear_energy = 6 - as.integer(q1b)) %>%
  mutate(support_new_nuclear = 6 - as.integer(q3a)) %>%
  mutate(Female = case_when(
    profile_gender == "Female" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(Age = fct_relevel(as_factor(age_cb), "65+", "16-24", "25-34", "35-44", "45-54", "55-64")) %>%
  mutate(Country = as_factor(case_when(
    profile_GOR %in% c("Scotland", "Wales", "Northern Ireland") ~ profile_GOR,
    TRUE ~ "England"
  ))) %>%
  mutate(Country = fct_relevel(Country, "England", "Wales", "Scotland", "Northern Ireland")) %>%
  mutate(CountryWales = case_when(
    Country == "Wales" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(CountryNorthernIreland = case_when(
    Country == "Northern Ireland" ~ 1,
    TRUE ~ 0
  )) %>%
  mutate(Vote_2024 = as_factor(case_when(
    Country == "Scotland" & pastvote_ge_yks_4 == "Yes" ~ "SNP",
    Country == "Wales" & pastvote_ge_yks_5 == "Yes" ~ "Plaid Cymru",
    pastvote_ge_yks_1 == "Yes" ~ "Conservative",
    pastvote_ge_yks_2 == "Yes" ~ "Labour",
    pastvote_ge_yks_3 == "Yes" ~ "Liberal Democrat",
    pastvote_ge_yks_6 == "Yes" ~ "Reform UK",
    pastvote_ge_yks_7 == "Yes" ~ "Green",
    pastvote_ge_yks_8 == "Yes" ~ "Other",
    pastvote_ge_yks_9 == "Yes" ~ "DNV/Ref",
    TRUE ~ "DNV/Ref"
  ))) %>%
  mutate(Vote_2024 = fct_relevel(Vote_2024, "Labour", "Conservative", "Liberal Democrat", "Reform UK",
    "Green", "SNP", "Plaid Cymru", "Other")) %>%
  mutate(Scotland_Vote_2024 = as_factor(case_when(
    Country == "Scotland" & Vote_2024 == "SNP" ~ "Scotland-SNP",
    Country == "Scotland" ~ "Scotland-Other",
    TRUE ~ "Outside Scotland"
  ))) %>%
  mutate(Scotland_Vote_2024 = fct_relevel(Scotland_Vote_2024, "Outside Scotland", "Scotland-SNP", "Scotland-Other"))

xtabs(~ q1b + support_nuclear_energy, addNA = TRUE, data = data.1)
xtabs(~ q3a + support_new_nuclear, addNA = TRUE, data = data.1)
xtabs(~ profile_GOR + Country, addNA = TRUE, data = data.1)
xtabs(~ profile_gender + Female, addNA = TRUE, data = data.1)
xtabs(~ Age, addNA = TRUE, data = data.1)
xtabs(~ Vote_2024, addNA = TRUE, data = data.1)
xtabs(~ Vote_2024 + Country, addNA = TRUE, data = data.1)
xtabs(~ Scotland_Vote_2024 + Country, addNA = TRUE, data = data.1)

wtd.cor(x = (select(data.1, support_nuclear_energy, support_new_nuclear)),
        weight = data.1$weight) %>%
  pluck(., "correlation")

svycralpha(~support_nuclear_energy + support_new_nuclear,
           (svydesign(id = ~1, weights = ~weight, data = data.1)), na.rm = TRUE)

range_std <- function(x) {
  x <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

wtd_summary <- function(x, w, na.rm = TRUE, digits = 5) {
  if (na.rm) {
    valid <- !is.na(x) & !is.na(w)
    x <- x[valid]
    w <- w[valid]
  }
  qt <- Hmisc::wtd.quantile(x, w, probs = c(0, .25, .5, .75, 1))
  m <- Hmisc::wtd.mean(x, w)
  out <- c(qt[1:3], Mean = m, qt[4:5])
  names(out) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  round(out, digits)
}

data.2 <- data.1 %>%
  mutate(support_nuclear = support_nuclear_energy + support_new_nuclear) %>%
  mutate(support_nuclear = range_std(support_nuclear) - 0.5)

wtd_summary(data.2$support_nuclear, data.2$weight)
hist(data.2$support_nuclear)

## Regression models
mean_y <- weighted.mean(x = data.2$support_nuclear, weights = data.2$weight, na.rm = TRUE)
sd_y <- sqrt(Hmisc::wtd.var(data.2$support_nuclear, weights = data.2$weight, na.rm = TRUE))

priors <- c(set_prior(sprintf("normal(%f, %f)", mean_y, sd_y), class = "Intercept"),
            set_prior(sprintf("normal(0, %f)", 0.5 * sd_y), class = "b"),
            set_prior(sprintf("exponential(%f)", 1 / sd_y), class = "sigma"))

m1 <- brm(formula = support_nuclear | weights(weight) ~ Country + Female + Age + Vote_2024,
  data = data.2,
  family = gaussian,
  prior = priors,
  backend = "cmdstanr",
  chains = 5,
  cores = 5,
  threads = 3,
  iter = 12000,
  warmup = 2000,
  refresh = 1000,
  seed = 20261125,
  save_pars = save_pars(all = TRUE))
print(m1, digits = 3)
prior_summary(m1)
m1 <- add_criterion(m1, "loo")
posterior_m1 <- as_draws_df(m1)

m2 <- brm(formula = support_nuclear | weights(weight) ~ CountryWales + CountryNorthernIreland
  + Scotland_Vote_2024 + Female + Age,
  data = data.2,
  family = gaussian,
  prior = priors,
  backend = "cmdstanr",
  chains = 5,
  cores = 5,
  threads = 3,
  iter = 12000,
  warmup = 2000,
  refresh = 1000,
  seed = 20261125,
  save_pars = save_pars(all = TRUE))
print(m2, digits = 3)
prior_summary(m2)
m2 <- add_criterion(m2, "loo")
posterior_m2 <- as_draws_df(m2)

mean(posterior_m2$b_Scotland_Vote_2024ScotlandMSNP < 0)