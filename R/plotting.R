#' @export
plotROCCurve <- function(CE.true, CE.emp, tol.true, tol.emp, precision.recall=FALSE,
                         restrict=FALSE, avg.auc=FALSE)
{
  # Computes TPR and FPR from list of true effect matrices (CE.true) and list
  # of empirical effect matrices (CE.emp), and plots ROC curve. Only effects
  # bigger than the corresponding tol are counted.

  N <- length(CE.true)
  p <- nrow(CE.true[[1]])
  nmax <- p*(p-1)/2

  # Threshold effects
  CE.true.thresh <- lapply(CE.true, function(mat) abs(mat) > tol.true)
  CE.emp.thresh <- lapply(CE.emp, function(mat) abs(mat) > tol.emp)

  # Restrict to only cases with both effects and non-effects?
  if (restrict) {
    neffects <- sapply(CE.true.thresh, sum)
    indices <- which((neffects > 0) & (neffects < nmax))
  } else {
    indices <- 1:N
  }

  # True positives
  TP <- lapply(1:N, function(i) CE.true.thresh[[i]] & CE.emp.thresh[[i]])

  # False positives
  FP <- lapply(1:N, function(i) (!CE.true.thresh[[i]]) & CE.emp.thresh[[i]])

  # True positive rate
  tpr <- function(i, ord) {
    ordered <- sort(abs(CE.emp[[i]]), decreasing=TRUE)
    TP2 <- TP[[i]] & (abs(CE.emp[[i]]) > ordered[ord])
    if (! precision.recall) {
      val <- sum(TP2)/sum(CE.true.thresh[[i]])
      if (sum(CE.true.thresh[[i]]) == 0) val <- 0
    } else {
      FP2 <- FP[[i]] & (abs(CE.emp[[i]]) > ordered[ord])
      val <- sum(TP2) / (sum(TP2) + sum(FP2))
      if (sum(TP2) + sum(FP2) == 0) val <- 0
    }
    return(val)
  }
  TPR <- sapply(indices, function(i) c(sapply(1:nmax, function(ord) tpr(i, ord)), 1))
  if (precision.recall) TPR[nmax+1,] <- sapply(1:N, function(i) sum(TP[[i]])/nmax)
  TPR.mean <- rowMeans(TPR, na.rm=TRUE)

  # False positive rate
  fpr <- function(i, ord) {
    ordered <- sort(abs(CE.emp[[i]]), decreasing=TRUE)
    if (! precision.recall) {
      FP2 <- FP[[i]] & (abs(CE.emp[[i]]) > ordered[ord])
      # Since in this binary classification problem the samples are not completely
      # independent (there cannot be an (estimated) causal effect in both directions),
      # the denominator of the FPR is not the number of non-effects, but the maximum
      # number of estimatable effects, so that FPR is always in [0,1]
      val <- sum(FP2)/nmax
    } else {
      TP2 <- TP[[i]] & (abs(CE.emp[[i]]) > ordered[ord])
      FN <- sum(CE.true.thresh[[i]]) - sum(TP2)
      val <- sum(TP2) / (sum(TP2) + FN)
      if (sum(TP2) + FN == 0) val <- 0
    }
    return(val)
  }
  FPR <- sapply(indices, function(i) c(sapply(1:nmax, function(ord) fpr(i, ord)), 1))
  if (precision.recall) FPR[nmax+1,] <- sapply(1:N, function(i) 1)
  FPR.mean <- rowMeans(FPR, na.rm=TRUE)

  # Compute AUC
  if (avg.auc) {
    # Average of AUCs of individual curves
    auc <- mean(sapply(1:ncol(FPR), function(i) caTools::trapz(FPR[,i], TPR[,i])))
  } else {
    # AUC of averaged curves
    auc <- caTools::trapz(FPR.mean, TPR.mean)
  }
  print(paste("AUC:", auc))

  # Switch labels between FPR-TPR and precision-recall
  if (precision.recall) {
    main.part <- "Precision-Recall Curve"
    xlab <- "Mean Recall"
    ylab <- "Mean Precision"
  } else {
    main.part <- paste("ROC Curve (AUC=", round(auc, 4), ")", sep="")
    xlab <- "Mean FPR"
    ylab <- "Mean TPR"
  }
  # main.part <- "ROC Curves"

  # Plot ROC Curve
  if (length(indices) == 1) {
    title <- paste(main.part, " (", sum(CE.true.thresh[[indices[1]]]), " causal effects)", sep="")
  } else {
    title <- main.part
  }
  #dev.new()
  plot(0, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, type="n",
       main=title)
  abline(0,1)
  for (i in indices) lines(FPR[,i], TPR[,i], col="gray")
  lines(FPR.mean, TPR.mean, type="b")

  return(list(FPR=FPR, TPR=TPR, auc=auc))
}


#' @export
plotScoreCurves <- function(res, gt=NULL, data=NULL, n=NULL, maxIter=10, zoom=1,
                            equivalent.eps=0, factor=2, plot.steps=TRUE, plot.times=TRUE)
{
  # Plots score improvements for a run of greedy search with restarts (where
  # the last result is taken to be from forward search).
  #
  # Arguments:
  #   res - should be the return object of greedySearchRestart
  #   gt - ground truth object
  #   data - Data matrix
  #   pop.version - flag indicating whether population version was used
  #   maxIter - maxIter that was used

  R <- length(res) - 1

  if (! is.null(gt)) {
    # Compute neg. Entropy (score limit)
    score.limit <- -0.5 * log(det(2*pi*exp(1)*gt$covMat))

    # Compute score of ground truth
    score.gt <- computeGraphScore(gt$mg, gt$covMat, data, n, Inf)
  } else {
    score.limit <- 0
    score.gt <- 0
  }

  # Pull out scores
  scores <- lapply(res, function(obj) sapply(obj$states, function(state) state$score))
  score.greedy <- max(sapply(scores, function(scorevec) max(scorevec)))

  xmax <- max(sapply(1:(R+1), function(i) res[[i]]$it))
  ymin <- min(sapply(1:(R+1), function(i) min(c(score.limit, score.gt, scores[[i]][scores[[i]]>-Inf]))))
  ymax <- max(sapply(1:(R+1), function(i) max(c(score.limit, score.gt, scores[[i]]))))
  ymin <- ymax - zoom*(ymax-ymin)
  if (equivalent.eps > 0) {
    ymin <- score.gt - factor*equivalent.eps
    ymax <- score.gt + factor*equivalent.eps
  }

  # Plot scores vs steps
  if (plot.steps) {
    dev.new()
    plot(scores[[1]], xlim=c(0, xmax), ylim=c(ymin, ymax), type="b", xlab="Steps",
         ylab="Score", main=signif(score.gt-score.greedy, 4))
    if (R > 1) {
      for (i in 2:R) lines(scores[[i]], col=i)
      for (i in 2:R) points(scores[[i]], col=i)
    }
    lines(scores[[R+1]], col=R+1)
    points(scores[[R+1]], col=R+1, pch=2)
    abline(score.limit, 0, lty=2, col="red", lwd=2)
    abline(score.gt, 0, lwd=2)
    legend("bottomright", legend=c("Neg. Entropy", "Score of true BAP"), lty=c(2,1),
           col=c("red", "black"), lwd=2)
    if (equivalent.eps > 0) {
      # Dashed lines indicating epsilon-environment
      abline(score.gt - equivalent.eps, 0, lty=2)
      abline(score.gt + equivalent.eps, 0, lty=2)
    }
  }

  if (plot.times) {
    # Pull out times
    times <- lapply(res, function(obj) sapply(obj$states, function(state) state$ct))

    # Plot scores vs time
    dev.new()
    xmax <- max(sapply(1:(R+1), function(i) max(times[[i]])))
    plot(times[[1]], scores[[1]], xlim=c(0, xmax), ylim=c(ymin, ymax), type="b",
         xlab="Time", ylab="Score", main=signif(score.gt-score.greedy, 4))
    if (R > 1) {
      for (i in 2:R) points(times[[i]], scores[[i]], col=i)
      for (i in 2:R) lines(times[[i]], scores[[i]], col=i)
    }
    points(times[[R+1]], scores[[R+1]], col=R+1, pch=2)
    lines(times[[R+1]], scores[[R+1]], col=R+1)
    abline(score.limit, 0, lty=2, col="red", lwd=2)
    abline(score.gt, 0, lwd=2)
    legend("bottomright", legend=c("Neg. Entropy", "Score of true BAP"), lty=c(2,1),
           col=c("red", "black"), lwd=2)
    if (equivalent.eps > 0) {
      # Dashed lines indicating epsilon-environment
      abline(score.gt - equivalent.eps, 0, lty=2)
      abline(score.gt + equivalent.eps, 0, lty=2)
    }
  }

  return(list(score.limit=score.limit, score.gt=score.gt,
              score.greedy=score.greedy))
}


#' @export
pl <- function(mg, tcltk=FALSE, ...) ggm::plotGraph(mg, layout=igraph::layout.circle, tcltk=tcltk, ...)
