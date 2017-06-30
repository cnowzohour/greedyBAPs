#' Plots ROC curve from estimated and true causal effects
#'
#' @param CE.gt List of minimal absolute causal effect matrices E (one per simulation run)
#'   for the ground truths, where (E)_i,j is the causal effect on j on i (as returned by
#'   causalEffectsSimulation)
#' @param CE.greedy same for the greedy search results
#' @param tol Tolerance for causal effects (everythin smaller than this will be thresholded)
#' @param precision.recall Plot precision-recall instead of FPR / TPR?
#' @param avg.auc
#'   If TRUE: Average of AUCs of individual curves
#'   If FALSE: AUC of averaged curves
#' @return List of three elements:
#'   FPR - Matrix with false positive rates for each of the simulation runs (one column per run)
#'   TPR - Matrix with true positive rates for each of the simulation runs (one column per run)
#'   auc - Averaged AUC
#'
#' @export
plotROCCurve <- function(
  CE.gt,
  CE.greedy,
  tol = 1e-8,
  precision.recall = FALSE,
  avg.auc = TRUE
) {
  # Computes TPR and FPR from list of true effect matrices (CE.gt) and list
  # of empirical effect matrices (CE.greedy), and plots ROC curve. Only effects
  # bigger than the corresponding tol are counted.

  restrict <- FALSE

  N <- length(CE.gt)
  p <- nrow(CE.gt[[1]])
  nmax <- p*(p-1)/2

  # Threshold effects
  CE.gt.thresh <- lapply(CE.gt, function(mat) abs(mat) > tol)
  CE.greedy.thresh <- lapply(CE.greedy, function(mat) abs(mat) > tol)

  # Restrict to only cases with both effects and non-effects?
  if (restrict) {
    neffects <- sapply(CE.gt.thresh, sum)
    indices <- which((neffects > 0) & (neffects < nmax))
  } else {
    indices <- 1:N
  }

  # True positives
  TP <- lapply(1:N, function(i) CE.gt.thresh[[i]] & CE.greedy.thresh[[i]])

  # False positives
  FP <- lapply(1:N, function(i) (!CE.gt.thresh[[i]]) & CE.greedy.thresh[[i]])

  # True positive rate
  tpr <- function(i, ord) {
    ordered <- sort(abs(CE.greedy[[i]]), decreasing=TRUE)
    TP2 <- TP[[i]] & (abs(CE.greedy[[i]]) > ordered[ord])
    if (! precision.recall) {
      val <- sum(TP2)/sum(CE.gt.thresh[[i]])
      if (sum(CE.gt.thresh[[i]]) == 0) val <- 0
    } else {
      FP2 <- FP[[i]] & (abs(CE.greedy[[i]]) > ordered[ord])
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
    ordered <- sort(abs(CE.greedy[[i]]), decreasing=TRUE)
    if (! precision.recall) {
      FP2 <- FP[[i]] & (abs(CE.greedy[[i]]) > ordered[ord])
      # Since in this binary classification problem the samples are not completely
      # independent (there cannot be an (estimated) causal effect in both directions),
      # the denominator of the FPR is not the number of non-effects, but the maximum
      # number of estimatable effects, so that FPR is always in [0,1]
      val <- sum(FP2)/nmax
    } else {
      TP2 <- TP[[i]] & (abs(CE.greedy[[i]]) > ordered[ord])
      FN <- sum(CE.gt.thresh[[i]]) - sum(TP2)
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
    title <- paste(main.part, " (", sum(CE.gt.thresh[[indices[1]]]), " causal effects)", sep="")
  } else {
    title <- main.part
  }
  #dev.new()
  plot(0, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, type="n",
       main=title)
  abline(0,1)
  for (i in indices) lines(FPR[,i], TPR[,i], col="gray")
  lines(FPR.mean, TPR.mean, type="b")

  list(FPR = FPR, TPR = TPR, auc = auc)
}


#' Plots scores per greedy step against cumulative time
#'
#' @param scores List (one element per greedy search) with score vectors (as returned by greedySearch)
#' @param times List (one element per greedy search) with cumulative time vectors (as returned by greedySearch)
#' @param plot.steps Should a second plot of scores vs steps be generated?
#' @param ylim Vector (ymin, ymax) - if this is NULL, it is generated from scores
#' @param title.start String with title prefix (the max score is appended to this)
#'
#' @export
plotScoreCurves <- function(
  scores,
  times,
  plot.steps = FALSE,
  ylim = NULL,
  title.start = "Max Score: "
) {

  R <- length(scores)

  # Max score
  score.greedy <- max(sapply(scores, function(scorevec) max(scorevec)))

  # y-axis limit
  if (is.null(ylim))
    ylim <- c(
      min(sapply(1:R, function(i) min(scores[[i]][scores[[i]]>-Inf]))),
      max(sapply(1:R, function(i) max(scores[[i]])))
    )

  # Draw plot
  draw.plot <- function(x) {
    xmax <- max(sapply(1:R, function(i) max(x[[i]])))
    par(mar=c(2.5, 2.5, 2.5, 0.5))  # bottom, left, top, right
    plot(x[[1]], scores[[1]], xlim=c(0, xmax), ylim=ylim, type="b",
         xlab="", ylab="", main=paste(title.start, signif(score.greedy, 4), sep=""))
    if (R > 1) {
      for (i in 2:R) points(x[[i]], scores[[i]], col=i)
      for (i in 2:R) lines(x[[i]], scores[[i]], col=i)
    }
  }

  draw.plot(times)
  if (plot.steps) {
    dev.new()
    draw.plot(lapply(scores, function(scorevec) 1:length(scorevec)))
  }
}


#' Plots a BAP using circular layout
#'
#' @param bap p-by-p adjacency matrix
#' @param ... any extra arguments to pass on to ggm::plotGraph
#'
#' @export
plotBAP <- function(
  bap,
  tcltk=FALSE,
  ...
) ggm::plotGraph(bap, layout=igraph::layout.circle, tcltk=tcltk, ...)
