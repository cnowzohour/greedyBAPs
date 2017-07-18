library(greedyBAPs)


mc.cores <- 20

# p <- 10
p <- 5
max.in.degree <- 2
n <- 1000
N <- 100
# R <- 100
R <- 10
equivalent.eps <- 1e-10
maxIter <- 10
maxSteps <- 100

res <- causalEffectsSimulation(
  N,
  p,
  n,
  mc.cores,
  max.in.degree,
  R,
  maxSteps,
  maxIter,
  equivalent.eps
)

roc <- plotROCCurve(res$CE.gt, res$CE.greedy, 1e-8, avg.auc = TRUE)
