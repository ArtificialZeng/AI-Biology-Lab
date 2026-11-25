options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("HonestDiD", quietly=TRUE)) {
  install.packages("HonestDiD", dependencies=c("Depends", "Imports", "LinkingTo"))
}
library(HonestDiD)
cat("INSTALLED_HONESTDID_VERSION:", as.character(packageVersion("HonestDiD")), "\n")
data("BCdata_EventStudy")
res_sd <- createSensitivityResults(
  betahat = BCdata_EventStudy$betahat,
  sigma = BCdata_EventStudy$sigma,
  numPrePeriods = 4,
  numPostPeriods = 4,
  Mvec = c(0, 0.01, 0.02, 0.05)
)
res_rm <- createSensitivityResults_relativeMagnitudes(
  betahat = BCdata_EventStudy$betahat,
  sigma = BCdata_EventStudy$sigma,
  numPrePeriods = 4,
  numPostPeriods = 4,
  Mbarvec = c(0.5, 1.0)
)
print(res_sd)
print(res_rm)
