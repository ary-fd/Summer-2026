## ---------------------------------------------------------------------------
## compute_benchmark_metrics.R
##
## Reads the 300 saved posterior samples + true parameter values from the
## ABC benchmarking run, and computes:
##   1. Prediction error metrics: bias, RMSE, coverage (per parameter)
##   2. SBC rank statistics + rank histograms (per parameter)
##
## Inputs (produced by run_abc_benchmark.R):
##   abc_results_10M/benchmark_posteriors.rds   - list of 300 posterior matrices
##   abc_results_10M/benchmark_params_true.rds  - data.frame, 300 rows x 4 params
##
## Outputs:
##   abc_results_10M/benchmark_metrics.csv      - summary table (bias/RMSE/coverage)
##   abc_results_10M/benchmark_plots.pdf        - scatter plots + SBC rank histograms
## ---------------------------------------------------------------------------

posteriors  <- readRDS("abc_results_10M/benchmark_posteriors.rds")
params_true <- readRDS("abc_results_10M/benchmark_params_true.rds")

param_names <- colnames(params_true)
n_test      <- nrow(params_true)
n_params    <- length(param_names)

cat("Loaded", n_test, "test cases,", n_params, "parameters:",
    paste(param_names, collapse = ", "), "\n")

## ---------------------------------------------------------------------------
## 1. Prediction error metrics: bias, RMSE, coverage
## ---------------------------------------------------------------------------

point_est <- matrix(NA, nrow = n_test, ncol = n_params,
                     dimnames = list(NULL, param_names))
ci_lower <- matrix(NA, nrow = n_test, ncol = n_params,
                    dimnames = list(NULL, param_names))
ci_upper <- matrix(NA, nrow = n_test, ncol = n_params,
                    dimnames = list(NULL, param_names))

for (i in seq_len(n_test)) {
    post_i <- posteriors[[i]]

    if (is.null(post_i) || nrow(post_i) == 0) {
        warning(sprintf("Test case %d has empty posterior — skipping", i))
        next
    }

    for (p in param_names) {
        draws <- post_i[, p]
        point_est[i, p] <- median(draws)
        ci_lower[i, p]  <- quantile(draws, 0.05)
        ci_upper[i, p]  <- quantile(draws, 0.95)
    }
}

metrics <- data.frame(
    parameter = param_names,
    bias      = NA_real_,
    rmse      = NA_real_,
    coverage_90 = NA_real_,
    correlation = NA_real_
)

for (j in seq_along(param_names)) {
    p <- param_names[j]
    truth <- params_true[[p]]
    est   <- point_est[, p]

    error <- est - truth
    metrics$bias[j]        <- mean(error, na.rm = TRUE)
    metrics$rmse[j]         <- sqrt(mean(error^2, na.rm = TRUE))
    metrics$coverage_90[j]  <- mean(truth >= ci_lower[, p] & truth <= ci_upper[, p],
                                     na.rm = TRUE)
    metrics$correlation[j]  <- cor(est, truth, use = "complete.obs")
}

print(metrics)
write.csv(metrics, "abc_results_10M/benchmark_metrics.csv", row.names = FALSE)
cat("Saved metrics table to abc_results_10M/benchmark_metrics.csv\n")

## ---------------------------------------------------------------------------
## 2. SBC: rank statistics
## ---------------------------------------------------------------------------

post_sizes <- sapply(posteriors, function(x) if (is.null(x)) NA else nrow(x))
L <- min(post_sizes, na.rm = TRUE)
cat("Posterior sizes range from", min(post_sizes, na.rm = TRUE),
    "to", max(post_sizes, na.rm = TRUE), "- subsampling all to L =", L, "\n")

ranks <- matrix(NA, nrow = n_test, ncol = n_params,
                 dimnames = list(NULL, param_names))

set.seed(42)
for (i in seq_len(n_test)) {
    post_i <- posteriors[[i]]
    if (is.null(post_i) || nrow(post_i) < L) next

    idx <- sample(seq_len(nrow(post_i)), L)
    post_i_sub <- post_i[idx, , drop = FALSE]

    for (p in param_names) {
        truth <- params_true[[p]][i]
        ranks[i, p] <- sum(post_i_sub[, p] < truth)
    }
}

ranks_df <- as.data.frame(ranks)
write.csv(ranks_df, "abc_results_10M/benchmark_sbc_ranks.csv", row.names = FALSE)
cat("Saved SBC ranks to abc_results_10M/benchmark_sbc_ranks.csv\n")

## ---------------------------------------------------------------------------
## 3. Plots: scatter (prediction error) + histogram (SBC) per parameter
## ---------------------------------------------------------------------------

pdf("abc_results_10M/benchmark_plots.pdf", width = 10, height = 5)

for (p in param_names) {
    par(mfrow = c(1, 2))

    truth <- params_true[[p]]
    est   <- point_est[, p]
    lims  <- range(c(truth, est), na.rm = TRUE)

    plot(truth, est,
         xlim = lims, ylim = lims,
         xlab = paste("True", p), ylab = paste("Estimated", p),
         main = paste("Prediction error:", p),
         pch = 19, col = rgb(0, 0, 1, 0.4))
    abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)

    r <- ranks_df[[p]]
    hist(r, breaks = seq(0, L, length.out = 21),
         main = paste("SBC ranks:", p),
         xlab = "Rank", col = "steelblue", border = "white")
    abline(h = length(na.omit(r)) / 20, col = "red", lwd = 2, lty = 2)
}

dev.off()
cat("Saved plots to abc_results_10M/benchmark_plots.pdf\n")

cat("\nDone. Summary:\n")
print(metrics)