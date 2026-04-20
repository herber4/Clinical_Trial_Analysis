# Clinical trial statistics demo (simulated data)
# Run this script from the repo root with: Rscript R/analysis.R

set.seed(20260420)

# Recommended packages. survival and nlme are part of standard R distributions.
suppressPackageStartupMessages({
  library(survival)
  library(nlme)
})

out_dir <- file.path(getwd(), "results")
fig_dir <- file.path(getwd(), "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

n_per_arm <- 120
subject_id <- sprintf("S%03d", seq_len(2 * n_per_arm))
trt <- factor(rep(c("Placebo", "Active"), each = n_per_arm), levels = c("Placebo", "Active"))
strata <- factor(sample(c("Low", "High"), 2 * n_per_arm, replace = TRUE, prob = c(0.55, 0.45)))
age <- round(rnorm(2 * n_per_arm, mean = 58, sd = 11))
sex <- factor(sample(c("Female", "Male"), 2 * n_per_arm, replace = TRUE))
baseline_biomarker <- round(rlnorm(2 * n_per_arm, meanlog = log(100), sdlog = 0.25), 1)

# Primary continuous endpoint: change from baseline at week 12.
true_trt_effect <- ifelse(trt == "Active", -8, 0)
true_strata_effect <- ifelse(strata == "High", -3, 0)
change_wk12 <- round(true_trt_effect + true_strata_effect + rnorm(2 * n_per_arm, 0, 9), 1)
week12 <- round(baseline_biomarker + change_wk12, 1)

# Repeated measures data for an MMRM-style analysis.
visits <- c(0, 4, 8, 12)
long_dat <- do.call(rbind, lapply(seq_along(subject_id), function(i) {
  visit_effect <- c(0, -2.5, -5.5, -8.0)
  trt_boost <- if (trt[i] == "Active") c(0, -1.0, -2.5, -4.0) else c(0, 0, 0, 0)
  y <- baseline_biomarker[i] + visit_effect + trt_boost + rnorm(length(visits), 0, 4)
  data.frame(
    subject_id = subject_id[i],
    trt = trt[i],
    strata = strata[i],
    age = age[i],
    sex = sex[i],
    visit = visits,
    biomarker = round(y, 1),
    baseline_biomarker = baseline_biomarker[i]
  )
}))
long_dat$visit_f <- factor(long_dat$visit, levels = visits)
long_dat$subject_id <- factor(long_dat$subject_id)

# Binary responder endpoint.
response <- factor(ifelse(change_wk12 <= -5, "Responder", "Non-responder"), levels = c("Non-responder", "Responder"))

# Time-to-event endpoint.
# Active arm gets a longer event-free time on average.
rate <- ifelse(trt == "Active", 0.065, 0.10)
time_to_event <- rexp(2 * n_per_arm, rate = rate)
censor_time <- runif(2 * n_per_arm, min = 6, max = 18)
event_time <- pmin(time_to_event, censor_time)
event <- as.integer(time_to_event <= censor_time)

# Adverse events.
ae_types <- c("Headache", "Nausea", "Fatigue", "Diarrhea")
ae_counts <- do.call(rbind, lapply(seq_along(subject_id), function(i) {
  n_ae <- rpois(1, lambda = ifelse(trt[i] == "Active", 1.2, 0.8))
  if (n_ae == 0) return(NULL)
  data.frame(
    subject_id = subject_id[i],
    trt = trt[i],
    ae = sample(ae_types, n_ae, replace = TRUE),
    grade = sample(1:3, n_ae, replace = TRUE, prob = c(0.65, 0.25, 0.10)),
    serious = rbinom(n_ae, 1, prob = 0.08) == 1,
    stringsAsFactors = FALSE
  )
}))
if (is.null(ae_counts) || nrow(ae_counts) == 0) {
  ae_counts <- data.frame(subject_id = character(), trt = factor(), ae = character(), grade = integer(), serious = logical())
}

trial <- data.frame(
  subject_id = subject_id,
  trt = trt,
  strata = strata,
  age = age,
  sex = sex,
  baseline_biomarker = baseline_biomarker,
  change_wk12 = change_wk12,
  week12 = week12,
  response = response,
  event_time = round(event_time, 2),
  event = event
)

# 1) Baseline characteristics table.
baseline_summary <- data.frame(
  Variable = c("Age, mean (SD)", "Female, n (%)", "High strata, n (%)", "Baseline biomarker, mean (SD)"),
  Placebo = c(
    sprintf("%.1f (%.1f)", mean(age[trt == "Placebo"]), sd(age[trt == "Placebo"])),
    sprintf("%d (%.1f%%)", sum(sex[trt == "Placebo"] == "Female"), 100 * mean(sex[trt == "Placebo"] == "Female")),
    sprintf("%d (%.1f%%)", sum(strata[trt == "Placebo"] == "High"), 100 * mean(strata[trt == "Placebo"] == "High")),
    sprintf("%.1f (%.1f)", mean(baseline_biomarker[trt == "Placebo"]), sd(baseline_biomarker[trt == "Placebo"]))
  ),
  Active = c(
    sprintf("%.1f (%.1f)", mean(age[trt == "Active"]), sd(age[trt == "Active"])),
    sprintf("%d (%.1f%%)", sum(sex[trt == "Active"] == "Female"), 100 * mean(sex[trt == "Active"] == "Female")),
    sprintf("%d (%.1f%%)", sum(strata[trt == "Active"] == "High"), 100 * mean(strata[trt == "Active"] == "High")),
    sprintf("%.1f (%.1f)", mean(baseline_biomarker[trt == "Active"]), sd(baseline_biomarker[trt == "Active"]))
  ),
  stringsAsFactors = FALSE
)
write.csv(baseline_summary, file.path(out_dir, "baseline_summary.csv"), row.names = FALSE)

# 2) ANCOVA for continuous endpoint.
ancova <- lm(change_wk12 ~ trt + baseline_biomarker + strata + age + sex, data = trial)
ancova_sum <- summary(ancova)
ancova_est <- coef(ancova)["trtActive"]
ancova_se <- ancova_sum$coefficients["trtActive", "Std. Error"]
ancova_ci <- ancova_est + c(-1, 1) * qnorm(0.975) * ancova_se
ancova_p <- ancova_sum$coefficients["trtActive", "Pr(>|t|)"]
continuous_result <- data.frame(
  Endpoint = "Week 12 change from baseline",
  Estimator = "ANCOVA (Active vs Placebo)",
  Estimate = round(ancova_est, 2),
  CI_low = round(ancova_ci[1], 2),
  CI_high = round(ancova_ci[2], 2),
  P_value = signif(ancova_p, 3)
)
write.csv(continuous_result, file.path(out_dir, "continuous_endpoint.csv"), row.names = FALSE)

# 3) Logistic regression for response.
logit <- glm(response ~ trt + baseline_biomarker + strata + age + sex, data = trial, family = binomial())
logit_coef <- coef(logit)["trtActive"]
logit_se <- summary(logit)$coefficients["trtActive", "Std. Error"]
logit_or <- exp(logit_coef)
logit_ci <- exp(logit_coef + c(-1, 1) * qnorm(0.975) * logit_se)
logit_p <- summary(logit)$coefficients["trtActive", "Pr(>|z|)"]
response_result <- data.frame(
  Endpoint = "Responder at week 12",
  Estimator = "Logistic regression OR (Active vs Placebo)",
  OR = round(logit_or, 2),
  CI_low = round(logit_ci[1], 2),
  CI_high = round(logit_ci[2], 2),
  P_value = signif(logit_p, 3)
)
write.csv(response_result, file.path(out_dir, "binary_endpoint.csv"), row.names = FALSE)

# 4) Time-to-event: Kaplan-Meier and Cox model.
km <- survfit(Surv(event_time, event) ~ trt, data = trial)
cox <- coxph(Surv(event_time, event) ~ trt + strata + age + sex, data = trial)
cox_coef <- coef(cox)["trtActive"]
cox_se <- summary(cox)$coefficients["trtActive", "se(coef)"]
cox_hr <- exp(cox_coef)
cox_ci <- exp(cox_coef + c(-1, 1) * qnorm(0.975) * cox_se)
cox_p <- summary(cox)$coefficients["trtActive", "Pr(>|z|)"]
cox_result <- data.frame(
  Endpoint = "Time to first event",
  Estimator = "Cox PH HR (Active vs Placebo)",
  HR = round(cox_hr, 2),
  CI_low = round(cox_ci[1], 2),
  CI_high = round(cox_ci[2], 2),
  P_value = signif(cox_p, 3)
)
write.csv(cox_result, file.path(out_dir, "time_to_event.csv"), row.names = FALSE)

png(file.path(fig_dir, "km_curve.png"), width = 900, height = 650)
plot(km, col = c("black", "blue"), lwd = 2, xlab = "Months", ylab = "Event-free probability", mark.time = TRUE)
legend("topright", legend = levels(trt), col = c("black", "blue"), lwd = 2, bty = "n")
dev.off()

# 5) Repeated measures model (MMRM-like via nlme::lme).
mmrm <- lme(
  biomarker ~ trt * visit_f + baseline_biomarker + strata + age + sex,
  random = ~ 1 | subject_id,
  data = long_dat,
  na.action = na.omit,
  method = "REML"
)
mmrm_sum <- summary(mmrm)
mmrm_result <- capture.output(mmrm_sum)
writeLines(mmrm_result, file.path(out_dir, "mmrm_summary.txt"))

# Mean trajectory plot.
mean_traj <- aggregate(biomarker ~ trt + visit, data = long_dat, FUN = mean)
sd_traj <- aggregate(biomarker ~ trt + visit, data = long_dat, FUN = sd)
mean_traj$sd <- sd_traj$biomarker
png(file.path(fig_dir, "longitudinal_means.png"), width = 900, height = 650)
plot(NA, xlim = range(visits), ylim = range(mean_traj$biomarker - mean_traj$sd, mean_traj$biomarker + mean_traj$sd),
     xlab = "Visit (weeks)", ylab = "Biomarker", main = "Mean biomarker over time")
cols <- c("black", "blue")
for (i in seq_along(levels(trt))) {
  tr <- levels(trt)[i]
  sub <- mean_traj[mean_traj$trt == tr, ]
  lines(sub$visit, sub$biomarker, type = "b", pch = 16, lwd = 2, col = cols[i])
  arrows(sub$visit, sub$biomarker - sub$sd, sub$visit, sub$biomarker + sub$sd, angle = 90, code = 3, length = 0.05, col = cols[i])
}
legend("topright", legend = levels(trt), col = cols, lwd = 2, pch = 16, bty = "n")
dev.off()

# 6) Subgroup analysis for forest plot (treatment effect on continuous endpoint by strata and sex).
subgroups <- list(
  strata = levels(strata),
  sex = levels(sex)
)
forest_dat <- data.frame(Subgroup = character(), Level = character(), Estimate = numeric(), Low = numeric(), High = numeric(), stringsAsFactors = FALSE)
for (sg_name in names(subgroups)) {
  for (lvl in subgroups[[sg_name]]) {
    idx <- trial[[sg_name]] == lvl
    fit <- lm(change_wk12 ~ trt + baseline_biomarker + age + sex, data = trial[idx, ])
    ci <- confint(fit)["trtActive", ]
    forest_dat <- rbind(forest_dat, data.frame(
      Subgroup = sg_name,
      Level = lvl,
      Estimate = coef(fit)["trtActive"],
      Low = ci[1],
      High = ci[2]
    ))
  }
}
write.csv(forest_dat, file.path(out_dir, "subgroup_forest.csv"), row.names = FALSE)

png(file.path(fig_dir, "forest_plot.png"), width = 950, height = 650)
par(mar = c(5, 9, 4, 2))
plot(NA, xlim = range(c(forest_dat$Low, forest_dat$High)), ylim = c(0.5, nrow(forest_dat) + 0.5),
     yaxt = "n", ylab = "", xlab = "Treatment effect (Active - Placebo)")
axis(2, at = seq_len(nrow(forest_dat)), labels = paste(forest_dat$Subgroup, forest_dat$Level, sep = ": "), las = 1)
abline(v = 0, lty = 2)
for (i in seq_len(nrow(forest_dat))) {
  segments(forest_dat$Low[i], i, forest_dat$High[i], i, lwd = 2)
  points(forest_dat$Estimate[i], i, pch = 19)
}
dev.off()

# 7) Adverse events summary.
ae_summary <- do.call(rbind, lapply(levels(trt), function(arm) {
  arm_subj <- subject_id[trt == arm]
  arm_ae <- ae_counts[ae_counts$trt == arm, , drop = FALSE]
  data.frame(
    trt = arm,
    subject_count = length(arm_subj),
    subjects_with_any_ae = length(unique(arm_ae$subject_id)),
    events_total = nrow(arm_ae),
    subjects_with_serious_ae = if (nrow(arm_ae) == 0) 0 else length(unique(arm_ae$subject_id[arm_ae$serious])),
    stringsAsFactors = FALSE
  )
}))
ae_summary$percent_any_ae <- round(100 * ae_summary$subjects_with_any_ae / ae_summary$subject_count, 1)
ae_summary$percent_serious_ae <- round(100 * ae_summary$subjects_with_serious_ae / ae_summary$subject_count, 1)
write.csv(ae_summary, file.path(out_dir, "adverse_events_summary.csv"), row.names = FALSE)

cat("Done. Files written to:\n")
cat(" -", out_dir, "\n")
cat(" -", fig_dir, "\n")
