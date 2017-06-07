GenerateMG4 <- function(p, N=1, iter=6*p^2, p1=0.5, p2=0, p3=0, max.in.degree=Inf, names=paste("V", 1:p, sep="")) {
  # p1: P(directed | empty) = P(empty | directed)
  # 1-p1: P(bidirected | empty) = P(empty | bidirected)
  # p2: P(bidirected | directed) = P(directed | bidirected)
  # p3: P(direction switch | directed)

  # p2 must be < min(p1, 1-p1)
  # p3 must be < 1-p1-p2

  # DAG only version: p1=1

  # Initialize with empty matrix
  mg <- matrix(0, p, p)
  rownames(mg) <- names
  colnames(mg) <- names

  # MC-iteration
  res <- list()
  index <- 0
  for (i in 1:(N*iter)) {

    mg.old <- mg

    # Pick position uniformly at random
    pos <- sample(p^2, 1)
    pos.vec <- arrayInd(pos, .dim=c(p,p))
    pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

    # Compute in-degree at source node of position and its transpose
    low.indegree <- (length(which(mg[,pos.vec[2]] > 0)) < max.in.degree)
    low.indegree.t <- (length(which(mg[,pos.vec[1]] > 0)) < max.in.degree)

    # What is the current state at position?
    if ((mg[pos] == 0) && (mg[pos.trans] == 0) && (pos != pos.trans)) {

      # Empty - either add directed edge (and check if still acyclic) or add bidirected edge
      if (rbinom(1, 1, p1) == 1) {

        # Add directed edge (if acyclic and neighbourhood condition is fulfilled)
        if (low.indegree) {
          mg2 <- mg
          mg2[pos] <- 1
          tmp <- mg2
          tmp[tmp==100] <- 0
          if (ggm::isAcyclic(tmp)) mg <- mg2
        }
      } else {

        # Add bidirected edge (if neighbourhood condition is fulfilled)
        if (low.indegree && low.indegree.t) {
          mg[pos] <- 100
          mg[pos.trans] <- 100
        }
      }
    } else if (mg[pos] == 1) {

      ## Directed edge - remove or stay
      #if (rbinom(1, 1, p1) == 1) mg[pos] <- 0

      # Directed edge - remove, switch to bidirected, switch direction, or stay
      flag <- sample(4, 1, prob=c(p1, p2, p3, 1-p1-p2-p3))
      if (flag == 1) {
        mg[pos] <- 0
      } else if (flag == 2) {

        # Switch to bidirected (if neighbourhood condition is fulfilled)
        if (low.indegree.t) {
          mg[pos] <- 100
          mg[pos.trans] <- 100
        }
      } else if (flag == 3) {

        # Try to switch direction (if neighbourhood condition is fulfilled)
        if (low.indegree.t) {
          mg2 <- mg
          mg2[pos] <- 0
          mg2[pos.trans] <- 1
          tmp <- mg2
          tmp[tmp==100] <- 0
          if (ggm::isAcyclic(tmp)) mg <- mg2
        }
      }
    } else if (mg[pos] == 100) {

      # # Bidirected edge - remove or stay
      # if (rbinom(1, 1, 1-p1) == 1) {
      #   mg[pos] <- 0
      #   mg[pos.trans] <- 0
      # }

      # Bidirected edge - remove, switch to directed, or stay
      flag <- sample(3, 1, prob=c(1-p1, p2, p1-p2))
      if (flag == 1) {
        mg[pos] <- 0
        mg[pos.trans] <- 0
      } else if (flag == 2) {

        # Try to switch to directed
        mg2 <- mg
        mg2[pos] <- 1
        mg2[pos.trans] <- 0
        tmp <- mg2
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) mg <- mg2
      }
    }

    # Check if neighbourhood size condition is fulfilled
    #if (any(neighbourhoodSize(mg) > 1)) mg <- mg.old

    # "Harvest" if i is multiple of iter
    if (i %% iter == 0) {
      index <- index + 1
      res[[index]] <- mg
    }
  }

  return(res)
}


GenerateParams <- function(Bmask=NULL, Omegamask=NULL, mg=NULL, Bdist="snormal", Oscale=1) {
  # Randomly generate edge weights and error covariance matrix

  if (is.null(Bmask)) {
    Bmask <- mg
    Bmask[Bmask==100] <- 0
    Omegamask <- mg
    Omegamask[Omegamask!=100] <- 0
    Omegamask[Omegamask==100] <- 1
    diag(Omegamask) <- 1
  }

  # Fill edge weights with random numbers
  B <- Bmask
  indB <- which(B != 0)
  # Standard normal or uniform?
  if (Bdist=="snormal") {
    Bvals <- rnorm(length(indB))
  } else {
    Bvals <- runif(length(indB), 0.5, 0.9)
  }
  B[indB] <- Bvals

  # Repeat generating Omega until minimal eigenvalue is > 1e-6
  flag <- TRUE
  while (flag) {

    # Set covariances to standard normal random numbers
    Omega <- matrix(0, nrow(Omegamask), ncol(Omegamask))
    ind <- lower.tri(Omegamask) & (Omegamask==1)
    Omega[ind] <- rnorm(length(which(ind)))
    Omega <- Omega + t(Omega)
    colnames(Omega) <- colnames(Omegamask)
    rownames(Omega) <- rownames(Omegamask)

    # Set variances to rowsum of abs values plus chi^2(1)
    diag(Omega) <- rowSums(abs(Omega)) + rchisq(nrow(Omega), 1)

    # Check minimal eigenvalue
    if (eigen(Omega)$value[ncol(Omega)] > 1e-6) flag <- FALSE
  }

  # Rescale source nodes with Oscale
  indices <- which(sapply(1:ncol(Bmask), function(col) all(Bmask[,col]==0)))
  Omega[indices,] <- Omega[indices,] * sqrt(Oscale)
  Omega[,indices] <- Omega[,indices] * sqrt(Oscale)

  return(list(B=B, Omega=Omega))
}


isFaithful <- function(mg, Bhat, Ohat, Shat, faithful.eps) {
  B.ind <- which((abs(t(Bhat)) < faithful.eps) & (mg == 1))
  O.ind <- which((abs(Ohat) < faithful.eps) & (mg == 100))
  #S.ind <- which((abs(Shat) < faithful.eps) & (mg > 0))
  S.ind <- c()
  if ((length(B.ind > 0)) || (length(O.ind > 0)) || (length(S.ind > 0))) return(list(flag=FALSE, B.ind=B.ind, O.ind=O.ind, S.ind=S.ind))
  return(list(flag=TRUE))
}


GetSigma <- function(params) {
  # Returns true covariance matrix, given model parameters B and Omega

  p <- ncol(params$B)
  return(solve(diag(p)-t(params$B)) %*% params$Omega %*% t(solve(diag(p)-t(params$B))))
}


GenerateGroundTruth <- function(p, max.in.degree=Inf, Bdist="snormal", Oscale=1,
                                faithful.eps=0)
{
  res <- list()
  res$mg <- GenerateMG4(p, 1, max.in.degree=max.in.degree)[[1]]
  res$params <- GenerateParams(mg=res$mg, Bdist=Bdist, Oscale=Oscale)
  while (! isFaithful(res$mg, t(res$params$B), res$params$Omega, NULL, 10*faithful.eps)$flag) {
    print("Ground Truth not faithgful - regenerating...")
    res$params <- GenerateParams(mg=res$mg, Bdist=Bdist, Oscale=Oscale)
  }
  res$covMat <- GetSigma(res$params)
  return(res)
}


GenerateData <- function(n, params) {
  # Randomly sample from a path diagram with parameters given by params$B
  # and params$Omega

  eps <- mvtnorm::rmvnorm(n, sigma=params$Omega)
  data <- eps %*% t(solve(diag(nrow(params$B))-t(params$B)))
  colnames(data) <- colnames(params$B)

  return(data)
}
