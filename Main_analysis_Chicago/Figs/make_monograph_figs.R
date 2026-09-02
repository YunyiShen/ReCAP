library(ReCAP)
model_file <- paste0(
  "../../../Chicago-Deer-Culling/Population_simulation_scheme/",
  "6harv_4surv_3fec_equal.RData"
)
if (!file.exists(model_file)) {
  stop("App posterior result file not found: ", model_file)
}
load(model_file)

source("analysisLambda.R")
source("making_plot_harvest_rate.R")
source("plot_living_individuals.R")

source("population_change_vs_lambda.R")
source("plot_projected_schemes.R")
