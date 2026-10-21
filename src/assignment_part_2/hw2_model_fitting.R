library(survival)

# Step 1: Load Simulated Data
# The dataset contains 4 simulation scenarios:
# - Two sample sizes: n = 200 and n = 500 patients
# - Two loss-to-follow-up rates: 5% and 10%
# - Each scenario consists of R = 2000 independent simulation runs
# - Each patient observation has 6 variables:
#   days, age, comorbidities, etiology, observed time, and event indicator
data_dir <- if (dir.exists("src/assignment_part_2/data")) {
  "src/assignment_part_2/data"
} else {
  "data"
}
sim_data_path <- file.path(data_dir, "simulated_data.rds")
simulated_data <- readRDS(sim_data_path)

# Observed censoring rate per scenario (share of patients with status == 0)
censoring_rate <- vapply(simulated_data, function(d) 1 - mean(d[["status"]]), numeric(1))

# Observed censoring ratio per run (mean and std across 2000 runs per scenario)
censoring_per_run <- lapply(simulated_data, function(d) {
  tapply(d[["status"]] == 0, d[["run_id"]], mean)
})
censoring_mean <- vapply(censoring_per_run, mean, numeric(1))
censoring_std  <- vapply(censoring_per_run, sd, numeric(1))

# Observed loss to follow-up per run (censored before cutoff tau = 90)
loss_per_run <- lapply(simulated_data, function(d) {
  tapply(d[["status"]] == 0 & d[["time"]] < 90, d[["run_id"]], mean)
})
loss_mean <- vapply(loss_per_run, mean, numeric(1))
loss_std  <- vapply(loss_per_run, sd, numeric(1))

censoring_summary <- data.frame(
  scenario = names(simulated_data),
  loss_mean = loss_mean,
  loss_std = loss_std,
  total_mean = censoring_mean,
  total_std = censoring_std
)
cat("Observed censoring breakdown across simulation runs (mean and std):\n")
print(censoring_summary, digits = 4, row.names = FALSE)

# Helper function to fit a parametric model across all 2000 runs in a scenario.
# Stores coefficients, scale and standard errors for every run. Runs where
# survreg fails or does not converge are skipped (excluded) and counted in
# the "n_failed" attribute of the result.
fit_one_run <- function(run_data, model_formula, dist_name) {
  fit <- tryCatch(
    survreg(model_formula, data = run_data, dist = dist_name),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  se <- sqrt(diag(vcov(fit)))[names(coef(fit))]
  c(coef(fit), scale = fit[["scale"]], setNames(se, paste0("se_", names(se))))
}

fit_model_across_runs <- function(scenario_df, model_formula, dist_name) {
  runs <- split(scenario_df, scenario_df[["run_id"]])
  estimates_list <- lapply(runs, fit_one_run, model_formula, dist_name)
  failed <- vapply(estimates_list, is.null, logical(1))
  estimates_df <- as.data.frame(do.call(rbind, estimates_list[!failed]))
  names(estimates_df) <- sub("(Intercept)", "intercept", names(estimates_df), fixed = TRUE)
  attr(estimates_df, "n_failed") <- sum(failed)
  estimates_df
}

formula_full <- Surv(time, status) ~ days + age + nr_comor + damage
formula_no_age <- Surv(time, status) ~ days + nr_comor + damage

# Step 2: Fit Models for Research Question 1 (True Log-Logistic Model)
# Objective: Fit the correctly specified model across all 4 scenarios.
# For each of the 4 scenarios:
#   For each simulation run (1 to 2000):
#     Fit the log-logistic survival model using all covariates:
#       response: observed time and event indicator
#       predictors: days to treatment, age, comorbidities, and etiology
#     Extract and store the estimated regression coefficients
#     and shape parameter for this run
rq1_estimates <- lapply(
  simulated_data,
  fit_model_across_runs,
  model_formula = formula_full,
  dist_name = "loglogistic"
)

# Step 3: Fit Models for Research Question 2 (Misspecified Models)
# Objective: Fit alternative models to evaluate misspecification.

# Part A: Misspecified Baseline Hazard (Weibull Model)
#   For each simulation run:
#     Fit a Weibull survival model using all covariates:
#       response: observed time and event indicator
#       predictors: days to treatment, age, comorbidities, and etiology
#     Extract and store the estimated regression coefficients
rq2a_estimates <- lapply(
  simulated_data,
  fit_model_across_runs,
  model_formula = formula_full,
  dist_name = "weibull"
)

# Part B: Omitted Confounder (Leave Out Age)
#   For each simulation run:
#     Fit a log-logistic survival model omitting age:
#       response: observed time and event indicator
#       predictors: days to treatment, comorbidities, and etiology
#     Extract and store the estimated regression coefficients
rq2b_estimates <- lapply(
  simulated_data,
  fit_model_across_runs,
  model_formula = formula_no_age,
  dist_name = "loglogistic"
)

# Step 4: Save Simulation Results
# Save all extracted estimates across runs and scenarios to disk
model_estimates <- list(
  rq1_loglogistic = rq1_estimates,
  rq2a_weibull = rq2a_estimates,
  rq2b_no_age = rq2b_estimates,
  censoring_rate = censoring_rate,
  censoring_std = censoring_std,
  loss_mean = loss_mean,
  loss_std = loss_std,
  censoring_summary = censoring_summary
)

saveRDS(model_estimates, file = file.path(data_dir, "model_estimates.rds"))
