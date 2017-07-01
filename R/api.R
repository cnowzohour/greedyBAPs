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
#' @return List containing:
#'   final.bap - Adjacency matrix of "winning" BAP
#'   final.score - Score of "winning" BAP
#'   all.baps - List of lists (one per greedy search run) of adjacency matrices of all visited BAPs
#'   all.scores - List of vector (one per greedy search run) of scores of all visited BAPs
#'   times - List of vectors (one per greedy search run) of cumulative running time for every step
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

  res <- if (mc.cores > 1) {
    doParallel::registerDoParallel(cores = mc.cores)
    foreach::foreach(i = 1:(n.restarts + 1)) %dopar% { oneRep(i) }
  } else lapply(1:(n.restarts + 1), oneRep)

  i.best <- which.max(sapply(res, function(obj) obj$score))
  list(
    final.bap = res[[i.best]]$mg,
    final.score = res[[i.best]]$score,
    all.baps = lapply(res, function(obj) lapply(obj$states, function(state) state$mg)),
    all.scores = lapply(res, function(obj) sapply(obj$states, function(state) state$score)),
    times = lapply(res, function(obj) sapply(obj$states, function(state) state$ct))
  )
}


#' Uniform BAP generation using MCMC
#'
#' @param p Number of vertices
#' @param N Number of graphs to generate
#' @param dag Whether to generate DAGs or more general BAPs
#' @param max.in.degree Maximal allowed in-degree
#' @param names Names to give to the rows / columns of the adjacency matrix
#' @param iter Number if MCMC iterations between samples
#' @return List of N p-by-p adjacency matrices
#'
#' @export
generateUniformBAP <- function(
  p,
  N = 1,
  dag = FALSE,
  max.in.degree=Inf,
  names=paste0("V", 1:p),
  iter=6*p^2
) {
  GenerateMG4(
    p,
    N = N,
    iter = iter,
    p1 = ifelse(dag, 1, 0.5),
    max.in.degree = max.in.degree,
    names = names
  )
}


#' Randomly generate parameters (edge weights) for a given BAP
#'
#' We consider the following linear recursive model:
#'   X = B^t X + eps,
#' where (B)_{i,j} is the edge weight of edge i->j and corresponds to the (i,j)-entry
#' in the adjacency matrix of the graph, and eps has a N(0, Omega) distribution.
#'
#' @param bap p-by-p adjacency matrix
#' @return List containing:
#'   B - p-by-p matrix of directed edge weights
#'   Omega - p-by-p matrix of bidirected edge weights
#'
#' @export
generateBAPPArameters <- function(
  bap
) {
  GenerateParams(
    mg = bap
  )
}


#' Randomly generate data for a given BAP parametrisation
#'
#' @param n Number of samples to generate
#' @param params Named list containing:
#'   B - p-by-p matrix of directed edge weights
#'   Omega - p-by-p matrix of bidirected edge weights
#' @return n-by-p matrix
#'
#' @export
generateData <- function(
  n,
  params
) {
  GenerateData(
    n,
    params
  )
}


#' Simulate BAPs and compute causal effects of greedy search results and ground truth
#'
#' @param N Number of simulation runs
#' @param p Number of vertices
#' @param n Number of samples to generate
#' @param mc.cores Parallelism (number of cores to be used). If set to 1, no parallelism is used. Should be at most N.
#' @param max.in.degree Maximal allowed in-degree
#' @param n.restarts Number of restarts
#' @param max.steps Max number of greedy steps per run
#' @param max.iter.ricf Max number of iterations in RICF
#' @param equivalent.eps Maximal allowed score difference for two graphs to be considered equivalent
#' @param verbose Print extra log output?
#' @return List containing:
#'   CE.gt - List of minimal absolute causal effect matrices E (one per simulation run)
#'     for the ground truths, where (E)_i,j is the causal effect on j on i
#'   CE.greedy - same for the greedy search results
#'
#' @importFrom foreach %dopar%
#' @export
causalEffectsSimulation <- function(
  N,
  p,
  n,
  mc.cores = 1,
  max.in.degree = Inf,
  n.restarts = 0,
  max.steps = Inf,
  max.iter.ricf = 10,
  equivalent.eps = 1e-10,
  verbose = FALSE
) {

  oneRep <- function(i) {
    print(i)
    seed <- sample(1e6, 1)
    set.seed(seed)
    cat("Seed: ", seed, "\n")
    causalEffects(
      p,
      max.in.degree,
      "snormal",
      1,
      n,
      FALSE,
      n.restarts,
      equivalent.eps,
      max.iter.ricf,
      max.steps,
      fast = TRUE,
      verbose = verbose
    )
  }

  res <- if (mc.cores > 1) {
    doParallel::registerDoParallel(cores = mc.cores)
    foreach::foreach(i = 1:N) %dopar% { oneRep(i) }
  } else lapply(1:N, oneRep)

  CE.minabs.gt <- lapply(res, function(obj) obj$CE.minabs.gt)
  CE.minabs.greedy <- lapply(res, function(obj) obj$CE.minabs.greedy)

  list(CE.gt = CE.minabs.gt, CE.greedy = CE.minabs.greedy)
}
