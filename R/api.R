#' Greedy BAP search
#'
#' @param data n-by-p data matrix (n samples, p variables)
#' @param mg.start Starting graph(s). Must be either NULL, an adjacency matrix defining the mixed graph (if n.restarts is 0) or a list of adjacency matrices of length n.restarts+1. If NULL, a starting graph is uniformly sampled for each run.
#' @param n.restarts Number of restarts
#' @param mc.cores Parallelism (number of cores to be used). If set to 1, no parallelism is used. Should be at most n.restarts+1.
#' @param max.steps Max number of greedy steps per run
#' @param max.iter.ricf Max number of iterations in RICF
#' @param neighbourhood.size If set to a finite number, only this many (randomly sampled) neighbours are considered at each greedy step
#' @param eps.conv Minimal score improvement to continue search after each step
#' @param max.in.degree Max in-degree of considered graphs
#' @param dags.only Consider DAGs only?
#' @param direction Can be "forward" (only adding or changing edges), "backward" (only removing or changing edges), or "both" (adding, removing, and changing edges)
#' @param verbose Print extra log output?
#'
#' @importFrom foreach %dopar%
#' @export
greedySearch <- function(
  data,
  mg.start = NULL,
  n.restarts = 0,
  mc.cores = 1,
  max.steps = Inf,
  max.iter.ricf = 10,
  neighbourhood.size = Inf,  # max.pos
  eps.conv = 1e-12,
  max.in.degree = Inf,
  dags.only = FALSE,
  direction = "both",
  verbose = FALSE
) {

  cov.mat <- cov(data)
  p <- ncol(data)
  n <- nrow(data)

  bap.stash <- if (is.null(mg.start)) {
    p1 <- if (dags.only) 1 else 0.5
    GenerateMG4(p = p, N = n.restarts + 1, max.in.degree = max.in.degree, names = rownames(cov.mat), p1 = p1)
  } else {
    if (class(mg.start) == "list") {
      if (length(mg.start) == n.restarts + 1) mg.start else stop("Number of starting graphs must equal n.restarts+1")
    } else {
      if (class(mg.start) == "matrix" && all(dim(mg.start) == c(p, p))) list(mg.start) else stop("Invalid starting graphs")
    }
  }

  directionMap <- c(
    "forward" = 1,
    "backward" = 2,
    "both" = 3
  )

  oneRep <- function(i) fastGreedySearch(
    bap.stash[[i]],
    data = as.matrix(data),
    n = n,
    maxSteps = max.steps,
    direction = directionMap[direction],
    maxIter = max.iter.ricf,
    covMat = cov.mat,
    verbose = verbose,
    max.pos = neighbourhood.size,
    dags.only = dags.only,
    eps.conv = eps.conv
  )

  if (mc.cores > 1) {
    doParallel::registerDoParallel(cores = mc.cores)
    foreach::foreach(i = 1:(n.restarts + 1)) %dopar% { oneRep(i) }
  } else lapply(1:(n.restarts + 1), oneRep)
}
