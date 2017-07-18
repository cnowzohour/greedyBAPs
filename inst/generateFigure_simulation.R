library(greedyBAPs)


print(system.time(res <- causalEffectsSimulation(
  N = 100,
  p = 10,
  n = 1000,
  mc.cores = 20,
  max.in.degree = 2,
  n.restarts = 100,
  max.steps = 100,
  max.iter.ricf = 10,
  equivalent.eps = 1e-10
)))

postscript("simulation.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 5, width = 5)
roc <- plotROCCurve(res$CE.gt, res$CE.greedy, 1e-8, avg.auc = TRUE)
dev.off()
