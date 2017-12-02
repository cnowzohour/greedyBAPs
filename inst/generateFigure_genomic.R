library(greedyBAPs)


############################################
### Code for pre-computation of datasets ###
############################################

set.seed(2017)

for (i in 1:14) {
  data <- as.matrix(flowCytometry[flowCytometry$source == unique(flowCytometry$source)[i], 1:11])

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

  save(res.bap, file = paste0("baps", i, ".RData"))
  save(res.dag, file = paste0("dags", i, ".RData"))
}


############################################
###     Code for generating figures      ###
############################################

plotTile <- function(i) {
  load(sprintf("baps%s.RData", i))
  load(sprintf("dags%s.RData", i))

  plotBAP(res.bap$final.bap, noframe = TRUE, vertex.label.color = "black")
  text(0, 0, sprintf("Dataset %i (BAP)", i), cex = 1.5)

  plotBAP(res.dag$final.bap, noframe = TRUE, vertex.label.color = "black")
  text(0, 0, sprintf("Dataset %i (DAG)", i), cex = 1.5)
}

pdf("genomic1.pdf", onefile = FALSE, paper = "special", width = 16, height = 16)
par(mfrow = c(4, 4), mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 1:8) plotTile(i)
dev.off()

pdf("genomic2.pdf", onefile = FALSE, paper = "special", width = 16, height = 12)
par(mfrow = c(3, 4), mar = c(0.1, 0.1, 0.1, 0.1))
for (i in 9:14) plotTile(i)
dev.off()
