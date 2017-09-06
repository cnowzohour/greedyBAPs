library(greedyBAPs)


set.seed(2017)

postscript("genomic.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 10)
par(mfrow=c(4,4))

for (i in 1:8) {
  data <- as.matrix(flowCytometry[flowCytometry$source == unique(flowCytometry$source)[i], 1:11])
  data <- scale(data, scale = FALSE)

  res.bap <- greedySearch(
    data,
    n.restarts = 100,
    max.iter.ricf = 10,
    max.in.degree = Inf,
    mc.cores = 20,
    dags.only = FALSE
  )

  res.dag <- greedySearch(
    data,
    n.restarts = 100,
    max.iter.ricf = 10,
    max.in.degree = Inf,
    mc.cores = 20,
    dags.only = TRUE
  )

  ylim <- range(c(unlist(res.bap$all.scores), unlist(res.dag$all.scores)))

  plotScoreCurves(
    res.bap$all.scores,
    res.bap$times,
    ylim = ylim,
    title.start = paste0("BAP; Dataset ", i, "; Score=")
  )

  plotScoreCurves(
    res.dag$all.scores,
    res.dag$times,
    ylim = ylim,
    title.start = paste0("DAG; Dataset ", i, "; Score=")
  )
}

dev.off()
