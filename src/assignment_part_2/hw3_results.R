data_dir <- if (dir.exists("src/assignment_part_2/data")) "src/assignment_part_2/data" else "data"
model_estimates <- readRDS(file.path(data_dir, "model_estimates.rds"))
models <- c(rq1_loglogistic = "M1 (correct)", rq2a_weibull = "M2 Weibull", rq2b_no_age = "M3 no age")
truth <- c(days = -0.1, age = -0.02, nr_comor = -0.25, damage = -0.7)

# Performance measures
rows <- list()
for (m in names(models)) {
  for (s in names(model_estimates[[m]])) {
    est <- model_estimates[[m]][[s]]
    for (p in intersect(names(truth), names(est))) {
      b <- est[[p]]
      se <- est[[paste0("se_", p)]]
      tr <- truth[[p]]
      rows[[length(rows) + 1]] <- data.frame(
        model = models[[m]], scenario = s, coef = p, true = tr,
        bias = mean(b) - tr,
        rel_bias = 100 * (mean(b) - tr) / abs(tr),
        ese = sd(b),
        mse = mean((b - tr)^2),
        coverage = 100 * mean(abs(b - tr) <= 1.96 * se),
        failed = attr(est, "n_failed")
      )
    }
  }
}
perf <- do.call(rbind, rows)

print(round(model_estimates$censoring_rate, 3))
print(perf, digits = 3, row.names = FALSE)

# Table 1: beta_E and beta_K at 5% loss
tab <- perf[perf$coef %in% c("damage", "nr_comor") & grepl("loss05", perf$scenario),
            c("model", "scenario", "coef", "bias", "ese", "mse", "coverage")]
print(tab, digits = 3, row.names = FALSE)

# Figure 1
cols <- c("#2a78d6", "#eb6834", "#1baf7a")
fig_models <- c(rq1_loglogistic = "M1 log-logistic (correct)",
                rq2a_weibull = "M2 Weibull (wrong)",
                rq2b_no_age = "M3 no age (wrong)")
fig_scen <- c(n200_loss05 = "n=200\n5% lost", n200_loss10 = "n=200\n10% lost",
              n500_loss05 = "n=500\n5% lost", n500_loss10 = "n=500\n10% lost")
titles <- c(damage = "(a) Etiology effect, true value -0.70",
            nr_comor = "(b) Comorbidity effect, true value -0.25")

draw_panel <- function(p) {
  vals <- list(); at <- c(); cl <- c()
  for (j in seq_along(fig_scen)) {
    for (i in seq_along(fig_models)) {
      vals[[length(vals) + 1]] <- model_estimates[[names(fig_models)[i]]][[names(fig_scen)[j]]][[p]]
      at <- c(at, (j - 1) * 4 + i)
      cl <- c(cl, cols[i])
    }
  }
  boxplot(vals, at = at, col = adjustcolor(cl, 0.35), border = cl, outline = FALSE,
          xaxt = "n", las = 1, boxwex = 0.8, lwd = 1.2, main = titles[[p]],
          cex.main = 1, font.main = 1, ylab = "Estimate", frame = FALSE)
  abline(h = truth[[p]], lty = 2, col = "grey30", lwd = 1.5)
  axis(1, at = (seq_along(fig_scen) - 1) * 4 + 2, labels = fig_scen,
       tick = FALSE, padj = 0.6, cex.axis = 0.85)
}

png("results_figure.png", width = 3000, height = 1260, res = 300)
layout(matrix(c(1, 1, 2, 3), 2, byrow = TRUE), heights = c(0.12, 1))
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = c(fig_models, "True value"),
       fill = c(adjustcolor(cols, 0.35), NA), border = c(cols, NA),
       lty = c(NA, NA, NA, 2), col = c(NA, NA, NA, "grey30"),
       horiz = TRUE, bty = "n", cex = 0.9, merge = TRUE)
par(mar = c(4.2, 4.2, 2, 1))
draw_panel("damage")
draw_panel("nr_comor")
invisible(dev.off())
