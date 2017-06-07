#' @import ggm
icfmagCustom <- function(mag, S, tol = 1e-06, maxIter=Inf) {
    if (!is.matrix(S)) {
        stop("Second argument is not a matrix!")
    }
    if (dim(S)[1] != dim(S)[2]) {
        stop("Second argument is not a square matrix!")
    }
    if (min(eigen(S)[[1]]) <= 0) {
        stop("Second argument is not a positive definite matrix!")
    }
    p <- nrow(S)
    temp <- unmakeMG(mag)
    mag.ug <- temp$ug
    mag.dag <- temp$dg
    mag.bg <- temp$bg
    if (p == 1) {
        return(list(Sigmahat = S, Omegahat = S, Bhat = NULL,
            Lambdahat = NULL, iterations = 1))
    }
    Omega <- diag(diag(S))
    dimnames(Omega) <- dimnames(S)
    B <- diag(p)
    dimnames(B) <- dimnames(S)
    UG.part <- (1:p)[0 == apply(mag.dag + mag.bg, 2, sum)]
    if (length(UG.part) > 0) {
        Lambda.inv <- fitConGraph(mag.ug[UG.part, UG.part, drop = FALSE],
            S[UG.part, UG.part, drop = FALSE], p + 1, tol = tol)$Shat
        Omega[UG.part, UG.part] <- Lambda.inv
    }
    pa.each.node <- function(amat) {
        p <- nrow(amat)
        b <- vector(p, mode = "list")
        ip <- 1:p
        for (i in 1:p) b[[i]] <- ip[amat[, i] == 1]
        b
    }
    spo <- pa.each.node(mag.bg)
    nsp <- pa.each.node(cmpGraph(mag.bg))
    pars <- pa.each.node(mag.dag)
    i <- 0
    repeat {
        i <- i + 1
        Omega.old <- Omega
        B.old <- B
        for (v in setdiff(1:p, UG.part)) {
            parv <- pars[[v]]
            spov <- spo[[v]]
            if (length(spov) == 0) {
                if (length(parv) != 0) {
                  if (i == 1) {
                    B[v, parv] <- -S[v, parv] %*% solve(S[parv,
                      parv])
                    Omega[v, v] <- S[v, v] + B[v, parv] %*% S[parv,
                      v]
                  }
                }
            }
            else {
                if (length(parv) != 0) {
                  O.inv <- matrix(0, p, p)
                  O.inv[-v, -v] <- solve(Omega[-v, -v])
                  Z <- O.inv[spov, -v] %*% B[-v, ]
                  lpa <- length(parv)
                  lspo <- length(spov)
                  XX <- matrix(0, lpa + lspo, lpa + lspo)
                  XX[1:lpa, 1:lpa] <- S[parv, parv]
                  XX[1:lpa, (lpa + 1):(lpa + lspo)] <- S[parv,
                    ] %*% t(Z)
                  XX[(lpa + 1):(lpa + lspo), 1:lpa] <- t(XX[1:lpa,
                    (lpa + 1):(lpa + lspo)])
                  XX[(lpa + 1):(lpa + lspo), (lpa + 1):(lpa +
                    lspo)] <- Z %*% S %*% t(Z)
                  YX <- c(S[v, parv], S[v, ] %*% t(Z))
                  temp <- YX %*% solve(XX)
                  B[v, parv] <- -temp[1:lpa]
                  Omega[v, spov] <- temp[(lpa + 1):(lpa + lspo)]
                  Omega[spov, v] <- Omega[v, spov]
                  temp.var <- S[v, v] - temp %*% YX
                  Omega[v, v] <- temp.var + Omega[v, spov] %*%
                    O.inv[spov, spov] %*% Omega[spov, v]
                }
                else {
                  O.inv <- matrix(0, p, p)
                  O.inv[-v, -v] <- solve(Omega[-v, -v])
                  Z <- O.inv[spov, -v, drop=FALSE] %*% B[-v, ]
                  XX <- Z %*% S %*% t(Z)
                  YX <- c(S[v, ] %*% t(Z))
                  Omega[v, spov] <- YX %*% solve(XX)
                  Omega[spov, v] <- Omega[v, spov]
                  temp.var <- S[v, v] - Omega[v, spov] %*% YX
                  Omega[v, v] <- temp.var + Omega[v, spov] %*%
                    O.inv[spov, spov] %*% Omega[spov, v]
                }
            }
        }
        if (sum(abs(Omega.old - Omega)) + sum(abs(B.old - B)) < tol) break
        if (i >= maxIter) {
          i <- -1
          break
        }
    }
    Sigma <- solve(B) %*% Omega %*% solve(t(B))
    Lambda <- matrix(0, p, p)
    if (length(UG.part) > 0) {
        Lambda[-UG.part, -UG.part] <- Omega[-UG.part, -UG.part]
    }
    Omega[UG.part, UG.part] <- 0
    return(list(Sigmahat = Sigma, Bhat = B, Omegahat = Omega,
        Lambdahat = Lambda, iterations = i))
}

#' @export
fitAncestralGraphCustom <- function(amat, S, n, tol = 1e-06, maxIter=Inf) {
    nam <- rownames(S)
    nod <- rownames(amat)
    if (is.null(nod)) {
        stop("The adjacency matrix has no labels!")
    }
    if (!all(is.element(nod, nam)))
        stop("The nodes of the graph do not match the names of the variables")
    else sek <- intersect(nam, nod)
    S <- S[sek, sek, drop = FALSE]
    amat <- amat[sek, sek, drop = FALSE]
    temp <- icfmagCustom(amat, S, tol, maxIter=maxIter)
    p <- ncol(S)
    df <- p * (p - 1)/2 - sum(In(amat + t(amat)))/2

    SK = S %*% solve(temp$Sigmahat)
    dev <- ((sum(diag(SK)) - log(det(SK))) - p) * n

    if (is.null(temp$Bhat)) {
        Beta <- NULL
    }
    else {
        Beta <- temp$Bhat
    }
    return(list(Shat = temp$Sigmahat, Lhat = temp$Lambdahat,
        Bhat = Beta, Ohat = temp$Omegahat, dev = dev, df = df,
        it = temp$iterations))
}
