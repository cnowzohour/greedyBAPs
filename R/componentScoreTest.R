# mg <- t(cbind(
#   c(0, 1, 100),
#   c(0, 0, 1),
#   c(100, 0, 0)
# ))
# rownames(mg) <- 1:3
# colnames(mg) <- 1:3

# set.seed(3)
#
# p <- 10
# mg <- GenerateMG4(p)[[1]]
#
# pl(mg)
#
# n <- 1000
# params <- GenerateParams(mg=mg)
# data <- GenerateData(n, params)
# covMat <- cov(data)
#
# comps <- connectedComponents(mg)
#
# print(system.time(for (i in 1:1e4) {
# for (comp in comps$comp) {
#   # computeComponentScore(comp, mg, covMat, n = n, method = "trace")
#   computeComponentScore(comp, mg, covMat, n = n, method = "eigenvalue")
# }
# }))
