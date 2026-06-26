library(abc)

# Load data
params  <- read.csv("abc_results_10M/params.csv")
sim_csfs <- read.csv("abc_results_10M/csfs_sim.csv")
obs_csfs <- read.csv("abc_results_10M/csfs_observed.csv")

# Run ABC with neuralnet method (as in paper)
abc_out <- abc(
    target  = as.numeric(obs_csfs[1,]),
    param   = params,
    sumstat = sim_csfs,
    tol     = 0.005,       # keep closest 0.5% of simulations
    method  = "neuralnet"
)

# Summary of posterior
summary(abc_out)

# Plot posterior distributions
pdf("abc_results_10M/posterior_plots.pdf")
hist(abc_out, true = c(31250, 2150, 0.11, 25000))
dev.off()