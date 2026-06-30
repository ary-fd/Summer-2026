## ---------------------------------------------------------------------------
## compute_error_correlations.R
##
## Reads the saved benchmark posteriors and true parameters, computes
## residuals (theta_hat - theta_true) per parameter per test case, then
## looks at how errors co-vary across parameters.
##
## Outputs:
##   abc_results_10M/error_corr_matrix.csv   - correlation matrix of residuals
##   abc_results_10M/error_corr_plots.pdf    - pairs plot + correlation heatmap
## ---------------------------------------------------------------------------

posteriors  <- readRDS("abc_results_10M/benchmark_posteriors.rds")
params_true <- readRDS("abc_results_10M/benchmark_params_true.rds")

param_names <- colnames(params_true)
n_test      <- nrow(params_true)

## Build residual matrix: one row per test case, one column per parameter
residuals <- matrix(NA, nrow = n_test, ncol = length(param_names),
                    dimnames = list(NULL, param_names))

for (i in seq_len(n_test)) {
    post_i <- posteriors[[i]]
    if (is.null(post_i) || nrow(post_i) == 0) next
    for (p in param_names) {
        residuals[i, p] <- median(post_i[, p]) - params_true[[p]][i]
    }
}

residuals <- as.data.frame(residuals)

## Error correlation matrix
corr_mat <- cor(residuals, use = "complete.obs")
print(round(corr_mat, 3))
write.csv(corr_mat, "abc_results_10M/error_corr_matrix.csv")
cat("Saved error correlation matrix to abc_results_10M/error_corr_matrix.csv\n")

## ---------------------------------------------------------------------------
## Plots
## ---------------------------------------------------------------------------

pdf("abc_results_10M/error_corr_plots.pdf", width = 10, height = 10)

## 1. Pairs plot of residuals
## Off-diagonal: scatter of residual_i vs residual_j
## Diagonal: histogram of residuals for that parameter
## A strong off-diagonal pattern means those two parameters are being
## traded off against each other by ABC
pairs(
    residuals,
    main = "Residual pairs plot (theta_hat - theta_true)",
    pch  = 19,
    col  = rgb(0, 0, 1, 0.3),
    panel = function(x, y, ...) {
        points(x, y, pch = 19, col = rgb(0, 0, 1, 0.3))
        abline(h = 0, v = 0, col = "red", lty = 2)  # zero-error reference lines
        abline(lm(y ~ x), col = "darkred", lwd = 1.5)  # trend line
    }
)

## 2. Correlation heatmap
n <- length(param_names)
par(mar = c(6, 6, 4, 2))
image(
    1:n, 1:n, corr_mat[n:1, ],   # flip rows so diagonal reads top-left to bottom-right
    col  = colorRampPalette(c("#d73027", "white", "#4575b4"))(100),
    zlim = c(-1, 1),
    xaxt = "n", yaxt = "n",
    xlab = "", ylab = "",
    main = "Error correlation matrix"
)
axis(1, at = 1:n, labels = param_names, las = 2)
axis(2, at = 1:n, labels = rev(param_names), las = 1)

## Add correlation values as text
for (i in 1:n) {
    for (j in 1:n) {
        text(i, n + 1 - j, round(corr_mat[j, i], 2), cex = 1.2,
             col = ifelse(abs(corr_mat[j, i]) > 0.5, "white", "black"))
    }
}

dev.off()
cat("Saved plots to abc_results_10M/error_corr_plots.pdf\n")