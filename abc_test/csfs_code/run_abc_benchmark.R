library(abc)

# ---- Load fixed reference table ----
params   <- read.csv("abc_results_10M/params.csv")
sim_csfs <- read.csv("abc_results_10M/csfs_sim.csv")

# ---- Load the 300 test cases ----
params_test <- read.csv("abc_results_10M/params_test.csv")
csfs_test   <- read.csv("abc_results_10M/csfs_test.csv")

n_test <- nrow(params_test)
cat("Running ABC on", n_test, "test cases...\n")

posteriors <- vector("list", n_test)

for (i in 1:n_test) {
    target_i <- as.numeric(csfs_test[i, ])

    abc_out <- abc(
        target  = target_i,
        param   = params,
        sumstat = sim_csfs,
        tol     = 0.03,
        method  = "neuralnet"
    )

    posteriors[[i]] <- abc_out$adj.values

    if (i %% 20 == 0) {
        cat("Completed", i, "/", n_test, "\n")
        saveRDS(posteriors, "abc_results_10M/benchmark_posterior_0.03_checkpoint.rds")
    }
}

# Final save
saveRDS(posteriors,  "abc_results_10M/benchmark_posterior_0.03.rds")
saveRDS(params_test, "abc_results_10M/benchmark_params_true_0.03.rds")

cat("Done. Saved posteriors to benchmark_posterior_0.03.rds\n")