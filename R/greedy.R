componentHash <- function(comp, mg) {
  # Perfect hash implementation for components
  #
  # Each entry in mg gets mapped to 0, 1, or 2 (=100) plus a binary p-vector representing comp$cnodes

  # Get subgraph restricted to nodeset of component (plus parent-edges)
  nodes <- sort(c(comp$cnodes, comp$parents))
  sub <- mg
  sub[, comp$parents] <- 0  # remove edges between parents and from component to parents
  sub <- sub[nodes, nodes, drop=FALSE]  # restrict to node set

  hash <- sub
  dim(hash) <- NULL
  hash[hash==100] <- 2
  hash2 <- rep(0, nrow(mg))
  hash2[comp$cnodes] <- 1
  hash3 <- rep(0, nrow(mg))
  hash3[comp$parents] <- 1
  hash <- paste(c(hash, hash2, hash3), collapse="")
  return(hash)
}


connectedComponents <- function(mg) {
  # Returns a list of the connected components of mg, only considering
  # bidirected edges

  if (is.null(colnames(mg))) {
    nnames <- 1:ncol(mg)
  } else {
    nnames <- colnames(mg)
  }

  # Loop over nodes
  visited <- array()
  components <- list()
  node.comp <- rep(NA, nrow(mg))  # Array to save component number for each node
  comp.index <- 0
  for (i in 1:nrow(mg)) {
    # Check if node was already visited
    if (! i %in% visited) {
      # Put node i in queue
      queue <- c(i)
      # Increase component counter
      comp.index <- comp.index + 1
      # Loop over queue while there are elements in it
      while (length(queue) > 0) {
        # Check if first node in queue was already visited
        if (! queue[1] %in% visited) {
          # Add neighbours to queue
          queue <- c(queue, which(mg[,queue[1]]==100))
          # Add node to visited
          visited[length(visited)+1] <- queue[1]
          # Add node to components list
          if (length(components) < comp.index) {
            components[[comp.index]] <- queue[1]
          } else {
            components[[comp.index]] <- c(components[[comp.index]], queue[1])
          }
          node.comp[queue[1]] <- comp.index
        }
        # Remove node from queue
        queue <- queue[-1]
      }
    }
  }

  for (i in 1:length(components)) {

    # Find parent set of component (outside of component)
    pa <- sort(setdiff(unique(unlist(sapply(components[[i]], function(node) which(mg[,node]==1)))), components[[i]]))
    components[[i]] <- list(cnodes=sort(components[[i]]))
    components[[i]]$parents <- pa

    # Component name (UID)
    names(components)[i] <- componentHash(components[[i]], mg)

    # Node names
    names(components[[i]]$cnodes) <- nnames[components[[i]]$cnodes]
    names(components[[i]]$parents) <- nnames[components[[i]]$parents]
  }

  return(list(comp=components, node.comp=node.comp))
}


computeComponentScore <- function(comp, mg, covMat, data=NULL, n=NULL, maxIter=10,
                                  edge.penalty=1, faithful.eps=0, maxItInf=FALSE,
                                  method="trace")
{

  nodes <- sort(c(comp$cnodes, comp$parents))
  pa <- comp$parents

  # Get subgraph restricted to nodeset of component (plus parent-edges)
  sub <- mg
  sub[,pa] <- 0  # remove edges between parents and from component to parents
  sub <- sub[nodes, nodes, drop=FALSE]  # restrict to node set
  indices <- match(pa, nodes)  # New indices of parents in sub-matrix

  # Fit subgraph
  err.msg <- try(res <- fitAncestralGraphCustom(sub, covMat[nodes, nodes, drop=FALSE], Inf, maxIter=maxIter))
  if (class(err.msg) == "try-error") return(-Inf)
  if ((maxItInf) && (res$it == -1)) return(-Inf)
  n.edges <- length(which(sub==1)) + length(which(sub==100))/2

  if (method == "eigenvalue") {
    # Likelihood score computation with parent-adjustment
    CM <- solve(res$Shat)  # Estimated concentration matrix
    if (length(pa) > 0) {
      # Correct score for parents: subtract the corresponding inverse variances
      # from the diagonal entries
      if (length(pa) == 1) {
        CM[indices,indices] <- CM[indices,indices] - 1/res$Shat[indices,indices]
      } else {
        diag(CM[indices,indices]) <- diag(CM[indices,indices]) - 1/diag(res$Shat[indices,indices])
      }
    }
    if (! is.null(data)) {
      n <- nrow(data)
      e <- eigen(CM)
      score <- -0.5 * ( length(comp$cnodes) * log(2*pi) + log(det(res$Shat)/prod(diag(res$Shat[indices,indices,drop=FALSE]))) +
                          as.numeric(sum((data[,nodes,drop=FALSE] %*% e$vectors %*% diag(sqrt(as.complex(e$values)), nrow=length(e$values)))^2))/n ) -
        edge.penalty * log(n)/n * (n.edges + ncol(sub))
    } else {
      # Population version
      if (n == Inf) {
        penalty <- 0
      } else {
        penalty <- edge.penalty * log(n)/n * (n.edges + ncol(sub))
      }
      # Check if distribution is faithful
      if (is.null(res$Bhat)) res$Bhat <- cbind(c(1))
      if ((faithful.eps > 0) && (! isFaithful(sub, res$Bhat, res$Ohat, res$Shat, faithful.eps)$flag)) {
        score <- -Inf
      } else {
        score <- -0.5 * ( length(comp$cnodes) * log(2*pi) +
                            log(det(res$Shat)/prod(diag(covMat[pa,pa,drop=FALSE]))) +
                            sum(diag(CM %*% covMat[nodes, nodes, drop=FALSE])) ) - penalty
      }
    }
  } else {

    ### Trace method

    if (is.null(n)) n <- nrow(data)

    score <- -0.5 * ( length(comp$cnodes) * log(2*pi) + log(det(res$Shat)/prod(diag(res$Shat[indices,indices,drop=FALSE]))) +
                        (n-1)/n*(sum(diag(solve(res$Shat) %*% covMat[nodes,nodes,drop=FALSE]))-length(pa)) ) - edge.penalty * log(n)/n * (n.edges + ncol(sub))
  }

  return(score)
}


getSkeleton <- function(mg) {
  # Transforms mixed graph adjacency matrix into symmetric skeleton matrix

  s <- matrix(0, nrow(mg), ncol(mg))
  s[mg + t(mg) != 0] <- 1

  return(s)
}


getVStructures <- function(mg) {
  # Returns list of v-structures (unshielded triples) of mixed graph mg as a
  # string, with nodes ordered.

  p <- ncol(mg)
  sk <- getSkeleton(mg)

  vstruct <- ""
  res <- list()

  # Loop over (sink) nodes
  counter <- 1
  for (sink in 1:p) {

    # Find source nodes
    snodes <- which(mg[,sink] != 0)

    # Loop over source nodes
    if (length(snodes) > 1) {
      for (i in 1:(length(snodes)-1)) {

        snode1 <- snodes[i]

        # Candidates for other source nodes
        candnodes <- snodes[(i+1):length(snodes)]

        # Find other source nodes, which are not connected to first source node
        snodes2 <- candnodes[which(sk[snode1,candnodes] == 0)]

        # Add found v-structures to string
        if (length(snodes2) > 0) {
          vstruct <- paste(vstruct, paste(sink, ":", snode1, ":", snodes2, ";", collapse="", sep=""), sep="")
          for (j in 1:length(snodes2)) {
            res[[counter]] <- c(sink, snode1, snodes2[j])
            counter <- counter + 1
          }
        }
      }
    }
  }

  res <- sapply(res, cbind)

  return(list(string=vstruct, mat=res))
}


n.edges <- function(mg) length(which(mg==1)) + length(which(mg==100)) / 2


getColliderInvariantModels <- function(mg) {
  # Finds all models that have the same skeleton and collider structure as mg
  # (but not including mg).

  # Stage 1: Fix arrow heads and tails that cannot be changed without adding /
  #          removing colliders
  #  - Fix all columns with more than one entry in {1,10}
  #  - For each fixed tail with free head, fix its head
  #  - Create list of columns L with exactly 1 fixed head. Loop over L:
  #    - fix all tails and their heads
  #    - check if any of the new heads fulfill condition for L, append to L if so.

  # Stage 2: Try any combination for the free edge ends
  #  - no new colliders: number of entries in {1,10} per column stays the same,
  #    if this number was bigger than 1
  #  - no "lost" edges: number of edges stays the same

  # Get number of edges
  e <- n.edges(mg)

  # Change graph representation:
  #   0: No edge
  #   1: Arrow head (free)
  #   10: Arrow head (fixed)
  #   2: Arrow tail (free)
  #   20: Arrow tail (fixed)
  mg.old <- mg
  mg[t(mg)==1] <- 2
  mg[mg==100] <- 1

  # Stage 1
  # Fix all columns with more than one entry in {1,10}
  p <- ncol(mg)
  indices <- matrix(colSums((mg==1) | (mg==10)) > 1, p, p, byrow=TRUE) & (mg < 10)
  mg[indices] <- 10 * mg[indices]

  # Fix free heads of fixed tails
  indices <- which((mg == 1) & (t(mg) == 20))
  mg[indices] <- 10 * mg[indices]

  # Loop over columns with precisely 1 fixed head
  L <- which(colSums(mg == 10) == 1)
  i <- 0
  while (i < length(L)) {
    i <- i + 1
    l <- L[i]
    # Fix free tails
    indices <- which(mg[, l] == 2)
    mg[indices, l] <- 20
    # Fix corresponding heads
    mg[l, indices] <- 10
    # Check if any of the new columns needs to be added to L
    L <- append(L, setdiff(which(colSums(mg == 10) == 1), L))
  }

  # Stage 2
  mg.list <- list()
  # Position indices of free heads / tails
  pos.array <- which((mg == 1) | (mg == 2))
  if (length(pos.array) == 0) return(list(mg.old))
  # Matrix, where each column is a possible combination of heads / tails
  combinations <- sfsmisc::digitsBase(0:(2^length(pos.array)-1), 2)
  combinations[combinations==0] <- 2
  # Get column counts of entries in {1, 10}
  counts <- colSums((mg==1) | (mg==10))
  # Loop over combinations and check corresponding graph
  for (j in 1:ncol(combinations)) {
    # Build new adjacency matrix
    mg.new <- mg
    mg.new[pos.array] <- combinations[,j]
    # Check no new colliders
    counts.new <- colSums((mg.new==1) | (mg.new==10))
    if (any(counts.new[counts.new > 1] != counts[counts.new > 1])) next
    # Convert adjacency matrix to usual style
    mg.new[mg.new==2] <- 0
    mg.new[mg.new==20] <- 0
    mg.new[mg.new==10] <- 1
    mg.new[(mg.new==1) & (t(mg.new)==1)] <- 100
    # Check number of edges stays the same
    if (n.edges(mg.new) != e) next
    # Add graph to list
    mg.list[[length(mg.list)+1]] <- mg.new
  }
  return(mg.list)
}


getCausalEffects <- function(B) {
  # Returns causal effect matrix E, where (E)_i,j is the causal effect on j on i
  diag(B) <- 0
  return(t(solve(diag(ncol(B))-t(B))))
}


fastGreedySearch <- function(mg.start, data=NULL, n=NULL, maxSteps=Inf, direction=3,
                             maxIter=10, edge.penalty=1, verbose=TRUE, covMat=NULL,
                             faithful.eps=0, max.pos=Inf, dags.only=FALSE, eps.conv=1e-12)
{
  # 1) Find all connected components of mg.start
  # 2) Compute score for each component
  # 3) Progressively score each new candidate component and store scores for later use
  # 4) At each step, choose edges whiche maximizes score improvement
  #    a) if directed edge - improvement is straightforward (only the component of the sink node changes)
  #    b) if undirected edge - join two components and compute improvement accordingly

  # Take time
  t <- proc.time()[3]

  p <- ncol(mg.start)
  if (is.null(covMat)) covMat <- cov(data)
  if (is.null(colnames(mg.start))) {
    nnames <- 1:ncol(mg.start)
  } else {
    nnames <- colnames(mg.start)
  }

  # Find all connected components of mg.start
  res <- connectedComponents(mg.start)

  # Compute score for each component
  scores <- lapply(1:length(res$comp), function(i) computeComponentScore(res$comp[[i]], mg.start, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps))
  names(scores) <- names(res$comp)

  ## Greedy search
  i <- 1

  # Define state variables
  state <- list()
  state$mg <- mg.start
  state$comp <- res$comp
  state$node.comp <- res$node.comp
  state$score <- sum(unlist(scores))
  states <- list()

  # Take time for first "step"
  state$t <- proc.time()[3] - t
  state$ct <- state$t  # Cumulative time
  t <- proc.time()[3]

  while (i <= maxSteps) {

    len.old <- length(scores)

    # Explore candidate space

    # Forward steps
    cand.add <- list()
    i.a <- 0
    if (direction != 2) {

      # Loop over positions which are zero in both the adj mat and its
      # transpose (i.e. there is no arrow in either direction); we WLOG only
      # consider positions on the lower triangle
      positions <- which((state$mg+t(state$mg) == 0) & lower.tri(state$mg))
      if (length(positions) > 1) positions <- sample(positions, min(length(positions), max.pos))
      for (pos in positions) {

        # Get position corrdinates and transpose index
        pos.vec <- arrayInd(pos, c(p,p))
        pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

        # Check which component(s) the edge lies in
        c1 <- state$node.comp[pos.vec[1]]
        c2 <- state$node.comp[pos.vec[2]]

        if (! dags.only) {
          # Add bidirected edge (this is always possible)
          newstate <- state
          newstate$mg[pos] <- 100
          newstate$mg[pos.trans] <- 100
          # New component index: if components are joined, indices could shift
          c1.new <- c1 - as.numeric(c2 < c1)
          if (c1 != c2) {  # Need to join components
            # Join components (add nodes and parents of 2nd component to 1st one)
            newstate$comp[[c1]]$cnodes <- sort(c(state$comp[[c1]]$cnodes, state$comp[[c2]]$cnodes))
            newstate$comp[[c1]]$parents <- sort(setdiff(c(state$comp[[c1]]$parents, state$comp[[c2]]$parents), newstate$comp[[c1]]$cnodes))
            # Delete 2nd component (WARNING: this changes the indices of comp!)
            newstate$comp[[c2]] <- NULL
            # Update node.comp: shift all indices after c2 one down (since c2 was deleted)
            newstate$node.comp[newstate$node.comp > c2] <- newstate$node.comp[newstate$node.comp > c2] - 1
            # Update node.comp: change all indices of 2nd component to new index of 1st component
            newstate$node.comp[state$node.comp==c2] <- c1.new
          }
          # Update component name
          cname1 <- componentHash(newstate$comp[[c1.new]], newstate$mg)
          names(newstate$comp)[c1.new] <- cname1
          # Compute score of new component if necessary
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1.new]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c1]]] - as.numeric(c1!=c2)*scores[[names(state$comp)[c2]]] + scores[[cname1]]
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }

        # Test if adding a directed edge in the lower triangle still results in an acyclic graph
        newstate <- state
        newstate$mg[pos] <- 1
        tmp <- newstate$mg
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) {
          # c2 is the component of the sink of the edge - change this component
          if (c1 != c2) newstate$comp[[c2]]$parents <- sort(unique(c(newstate$comp[[c2]]$parents, pos.vec[1])))
          cname1 <- componentHash(newstate$comp[[c2]], newstate$mg)
          names(newstate$comp)[c2] <- cname1
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c2]]] + scores[[cname1]]
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }

        # Do same for upper triangle
        newstate <- state
        newstate$mg[pos.trans] <- 1
        tmp <- newstate$mg
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) {
          # c1 is the component of the sink of the edge - change this component
          if (c1 != c2) newstate$comp[[c1]]$parents <- sort(unique(c(newstate$comp[[c1]]$parents, pos.vec[2])))
          cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
          names(newstate$comp)[c1] <- cname1
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }
      }
    }

    # Backward steps
    cand.del <- list()
    i.d <- 0
    if (direction != 1) {

      # Loop over positions which are non-zero in either the adj mat or its
      # transpose; we WLOG only consider positions on the lower triangle
      positions <- which((state$mg+t(state$mg) != 0) & lower.tri(state$mg))
      if (length(positions) > 1) positions <- sample(positions, min(length(positions), max.pos))
      for (pos in positions) {

        # Get position corrdinates and transpose index
        pos.vec <- arrayInd(pos, c(p,p))
        pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

        # Delete edge
        newstate <- state
        newstate$mg[pos] <- 0
        newstate$mg[pos.trans] <- 0

        c1 <- state$node.comp[pos.vec[1]]

        # Check edge type
        if (state$mg[pos] == 100) {  # Bidirected edge

          # We might need to split up the component
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new components if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }
        } else {  # Directed edge

          c2 <- state$node.comp[pos.vec[2]]

          if (c1 == c2) {
            # Both edge-endpoints belong to the same component - update that
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          } else {
            # Edge-endpoints lie in different components
            if (state$mg[pos] == 1) {
              # c2 is the sink node component
              # Update parents of c2 component
              if (all(newstate$mg[pos.vec[1], newstate$comp[[c2]]$cnodes, drop=FALSE]==0)) newstate$comp[[c2]]$parents <- setdiff(newstate$comp[[c2]]$parents, pos.vec[1])
              # Update component name
              cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
              names(newstate$comp)[c2] <- cname2
              # Compute score if necessary
              if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
              # Set score of candidate model
              newstate$score <- state$score - scores[[names(state$comp)[c2]]] + scores[[cname2]]
            } else {
              # c1 is the sink node component
              # Update parents of c1 component
              if (all(newstate$mg[pos.vec[2], newstate$comp[[c1]]$cnodes, drop=FALSE]==0)) newstate$comp[[c1]]$parents <- setdiff(newstate$comp[[c1]]$parents, pos.vec[2])
              # Update component name
              cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
              names(newstate$comp)[c1] <- cname1
              # Compute score if necessary
              if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
              # Set score of candidate model
              newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
            }
          }
        }

        # Add new state to candidates
        i.d <- i.d + 1
        cand.del[[i.d]] <- newstate
      }
    }

    # Edge type changes
    cand.cha <- list()
    i.c <- 0
    #if (direction == 3) {
    if (direction <= 3) {

      # Loop over all '100's in the lower triangle of the adj mat
      positions <- which(lower.tri(state$mg) & (state$mg == 100))
      if (length(positions) > 1) positions <- sample(positions, min(length(positions), max.pos))
      for (pos in positions) {

        # Get position corrdinates and transpose index
        pos.vec <- arrayInd(pos, c(p,p))
        pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

        # Check which component the edge lies in
        c1 <- state$node.comp[pos.vec[1]]

        # Check if 1 in lower triangle is admissible
        newstate <- state
        newstate$mg[pos] <- 1
        newstate$mg[pos.trans] <- 0
        tmp <- newstate$mg
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) {
          # We might need to split up the component
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }
          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }

        # Check if 1 in upper triangle is admissible
        newstate <- state
        newstate$mg[pos] <- 0
        newstate$mg[pos.trans] <- 1
        tmp <- newstate$mg
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) {
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }
          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }
      }

      # Loop over all '1's in adj mat
      positions <- which(state$mg == 1)
      if (length(positions) > 1) positions <- sample(positions, min(length(positions), max.pos))
      for (pos in positions) {

        # Get position coordinates and transpose index
        pos.vec <- arrayInd(pos, c(p,p))
        pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

        c1 <- state$node.comp[pos.vec[1]]
        c2 <- state$node.comp[pos.vec[2]]

        # Check if 1 in transposed position is possible
        newstate <- state
        newstate$mg[pos] <- 0
        newstate$mg[pos.trans] <- 1
        tmp <- newstate$mg
        tmp[tmp==100] <- 0
        if (ggm::isAcyclic(tmp)) {
          if (c1 == c2) {
            # Both edge-endpoints belong to the same component - update that
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          } else {
            # Edge-endpoints lie in different components
            # c1 is the component of the old source node (new sink node)
            # c2 is the component of the old sink node (new source node)
            # Update parents of c1 component: add node
            newstate$comp[[c1]]$parents <- sort(unique(c(newstate$comp[[c1]]$parents, pos.vec[2])))
            # Update parents of c2 component: remove node
            if (all(newstate$mg[pos.vec[1], newstate$comp[[c2]]$cnodes, drop=FALSE]==0)) newstate$comp[[c2]]$parents <- setdiff(newstate$comp[[c2]]$parents, pos.vec[1])
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute score if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] - scores[[names(state$comp)[c2]]] + scores[[cname1]] + scores[[cname2]]
          }

          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }

        if (! dags.only) {
          # Add model with bidirected edge
          newstate <- state
          newstate$mg[pos] <- 100
          newstate$mg[pos.trans] <- 100
          # New component index: if components are joined, indices could shift
          c1.new <- c1 - as.numeric(c2 < c1)
          if (c1 != c2) {  # Need to join components
            # Join components (add nodes and parents of 2nd component to 1st one)
            newstate$comp[[c1]]$cnodes <- sort(c(state$comp[[c1]]$cnodes, state$comp[[c2]]$cnodes))
            newstate$comp[[c1]]$parents <- sort(setdiff(c(state$comp[[c1]]$parents, state$comp[[c2]]$parents), newstate$comp[[c1]]$cnodes))
            # Delete 2nd component (WARNING: this changes the indices of comp!)
            newstate$comp[[c2]] <- NULL
            # Update node.comp: shift all indices after c2 one down (since c2 was deleted)
            newstate$node.comp[newstate$node.comp > c2] <- newstate$node.comp[newstate$node.comp > c2] - 1
            # Update node.comp: change all indices of 2nd component to new index of 1st component
            newstate$node.comp[state$node.comp==c2] <- c1.new
          }
          # Update component name
          cname1 <- componentHash(newstate$comp[[c1.new]], newstate$mg)
          names(newstate$comp)[c1.new] <- cname1
          # Compute score of new component if necessary
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1.new]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c1]]] - as.numeric(c1!=c2)*scores[[names(state$comp)[c2]]] + scores[[cname1]]
        }

        i.c <- i.c + 1
        cand.cha[[i.c]] <- newstate
      }
    }

    # Compare candidates and select best one
    if (i.a > 0) {
      best.a <- which.max(sapply(1:i.a, function(j) cand.add[[j]]$score))
      if (length(best.a) > 0) {
        best.a.score <- cand.add[[best.a]]$score
      } else {
        best.a.score <- -Inf
      }
    } else {
      best.a <- NULL
      best.a.score <- -Inf
    }
    if (i.d > 0) {
      best.d <- which.max(sapply(1:i.d, function(j) cand.del[[j]]$score))
      if (length(best.d) > 0) {
        best.d.score <- cand.del[[best.d]]$score
      } else {
        best.d.score <- -Inf
      }
    } else {
      best.d <- NULL
      best.d.score <- -Inf
    }
    if (i.c > 0) {
      best.c <- which.max(sapply(1:i.c, function(j) cand.cha[[j]]$score))
      if (length(best.c) > 0) {
        best.c.score <- cand.cha[[best.c]]$score
      } else {
        best.c.score <- -Inf
      }
    } else {
      best.c <- NULL
      best.c.score <- -Inf
    }

    states[[i]] <- state

    # Break out of loop, if no improvement of score is possible
    if (max(c(best.a.score, best.d.score, best.c.score)) <= state$score + eps.conv) break

    # Otherwise change state to best candidate
    if ((best.a.score > best.d.score) && (best.a.score > best.c.score)) {
      state <- cand.add[[best.a]]
    } else {
      if (best.d.score > best.c.score) {
        state <- cand.del[[best.d]]
      } else {
        state <- cand.cha[[best.c]]
      }
    }

    # Take time for step
    state$t <- proc.time()[3] - t
    state$ct <- state$ct + state$t
    t <- proc.time()[3]

    i <- i + 1

    # Info output
    action <- c("ADD", "DEL", "CHA")[which.max(c(best.a.score, best.d.score, best.c.score))]
    lengths <- paste("(", length(cand.add), ", ", length(cand.del), ", ", length(cand.cha), ")", sep="")
    newscores <- round((length(scores)-len.old)/(length(cand.add)+length(cand.del)+length(cand.cha)), 4)
    if (verbose) print(paste(i, action, "; (adds, dels, changes) =", lengths, "; %(new scores) =",
                             newscores, "; steptime =", round(state$t, 1), "; #(components) =",
                             length(state$comp), "; max(component.size) =", max(sapply(state$comp, function(comp) length(comp$cnodes))),
                             "; score =", round(state$score, 4)))
    if ((action == "CHA") && verbose) print(paste("Change in v-structures?", getVStructures(states[[i-1]]$mg)$string == getVStructures(state$mg)$string))
  }

  if (i >= maxSteps) print("MAX STEPS ACHIEVED")

  return(list(mg=state$mg, score=state$score, it=i, states=states, cand.add=cand.add,
              cand.del=cand.del, cand.cha=cand.cha, scores=scores))
}


findEquivalentModels <- function(state, mgs.seen, scores, score, epsilon,
                                 depth=1, depth.max=1, cumtime=0, time.max=10,
                                 data=NULL, n=NULL, direction=3, maxIter=10,
                                 edge.penalty=1, verbose=TRUE, covMat=NULL,
                                 faithful.eps=0)
{
  # Search neighbourhood of mg (recursively) for models within epsilon band
  # of the original score.
  #
  # state: State containing the starting model
  # mgs.seen: Hashes of all mgs so far encountered
  # scores: All component scores so far computed
  # score: Reference score
  # epsilon: Width of score band to explore
  # depth: Current depth (=recursion level)
  # depth.max: Limit for depth
  # time.max: Limit for running time

  # Check if max depth or max time has been reached
  if ((depth > depth.max) || (cumtime > time.max)) return(list(res=list(), scores=scores, mgs.seen=mgs.seen, cumtime=cumtime))

  # Take time
  t <- proc.time()[3]

  p <- ncol(state$mg)

  # Make sure that current mg is in mgs.seen, and state has all needed info
  if (depth==1) {
    hash <- componentHash(list(cnodes=1:p), state$mg)
    if (! hash %in% mgs.seen) mgs.seen[length(mgs.seen)+1] <- hash

    if (length(state)==1) {
      tmp <- connectedComponents(state$mg)
      scores <- lapply(1:length(tmp$comp), function(i)
        computeComponentScore(tmp$comp[[i]], state$mg, covMat,
                              data, n, maxIter))
      names(scores) <- names(tmp$comp)
      state <- list(mg=state$mg, comp=tmp$comp, node.comp=tmp$node.comp, scores=scores,
                    score=sum(unlist(scores)))
      score <- state$score
    }
  }

  # Forward steps
  cand.add <- list()
  i.a <- 0
  if ((direction == 1) || (direction == 3)) {

    # Loop over positions which are zero in both the adj mat and its
    # transpose (i.e. there is no arrow in either direction); we WLOG only
    # consider positions on the lower triangle
    for (pos in which((state$mg+t(state$mg) == 0) & lower.tri(state$mg))) {

      # Get position corrdinates and transpose index
      pos.vec <- arrayInd(pos, c(p,p))
      pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

      # Check which component(s) the edge lies in
      c1 <- state$node.comp[pos.vec[1]]
      c2 <- state$node.comp[pos.vec[2]]

      # Add bidirected edge (this is always possible)
      newstate <- state
      newstate$mg[pos] <- 100
      newstate$mg[pos.trans] <- 100
      hash <- componentHash(list(cnodes=1:p), newstate$mg)
      if (! (hash %in% mgs.seen)) {
        mgs.seen[length(mgs.seen)+1] <- hash
        # New component index: if components are joined, indices could shift
        c1.new <- c1 - as.numeric(c2 < c1)
        if (c1 != c2) {  # Need to join components
          # Join components (add nodes and parents of 2nd component to 1st one)
          newstate$comp[[c1]]$cnodes <- sort(c(state$comp[[c1]]$cnodes, state$comp[[c2]]$cnodes))
          newstate$comp[[c1]]$parents <- sort(setdiff(c(state$comp[[c1]]$parents, state$comp[[c2]]$parents), newstate$comp[[c1]]$cnodes))
          # Delete 2nd component (WARNING: this changes the indices of comp!)
          newstate$comp[[c2]] <- NULL
          # Update node.comp: shift all indices after c2 one down (since c2 was deleted)
          newstate$node.comp[newstate$node.comp > c2] <- newstate$node.comp[newstate$node.comp > c2] - 1
          # Update node.comp: change all indices of 2nd component to new index of 1st component
          newstate$node.comp[state$node.comp==c2] <- c1.new
        }
        # Update component name
        cname1 <- componentHash(newstate$comp[[c1.new]], newstate$mg)
        names(newstate$comp)[c1.new] <- cname1
        # Compute score of new component if necessary
        if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1.new]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
        newstate$score <- state$score - scores[[names(state$comp)[c1]]] - as.numeric(c1!=c2)*scores[[names(state$comp)[c2]]] + scores[[cname1]]
        i.a <- i.a + 1
        cand.add[[i.a]] <- newstate
      }


      # Test if adding a directed edge in the lower triangle still results in an acyclic graph
      newstate <- state
      newstate$mg[pos] <- 1
      tmp <- newstate$mg
      tmp[tmp==100] <- 0
      if (ggm::isAcyclic(tmp)) {
        hash <- componentHash(list(cnodes=1:p), newstate$mg)
        if (! (hash %in% mgs.seen)) {
          mgs.seen[length(mgs.seen)+1] <- hash
          # c2 is the component of the sink of the edge - change this component
          if (c1 != c2) newstate$comp[[c2]]$parents <- sort(unique(c(newstate$comp[[c2]]$parents, pos.vec[1])))
          cname1 <- componentHash(newstate$comp[[c2]], newstate$mg)
          names(newstate$comp)[c2] <- cname1
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c2]]] + scores[[cname1]]
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }
      }

      # Do same for upper triangle
      newstate <- state
      newstate$mg[pos.trans] <- 1
      tmp <- newstate$mg
      tmp[tmp==100] <- 0
      if (ggm::isAcyclic(tmp)) {
        hash <- componentHash(list(cnodes=1:p), newstate$mg)
        if (! (hash %in% mgs.seen)) {
          mgs.seen[length(mgs.seen)+1] <- hash
          # c1 is the component of the sink of the edge - change this component
          if (c1 != c2) newstate$comp[[c1]]$parents <- sort(unique(c(newstate$comp[[c1]]$parents, pos.vec[2])))
          cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
          names(newstate$comp)[c1] <- cname1
          if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
          newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          i.a <- i.a + 1
          cand.add[[i.a]] <- newstate
        }
      }
    }
  }

  # Backward steps
  cand.del <- list()
  i.d <- 0
  if ((direction == 2) || (direction == 3)) {

    # Loop over positions which are non-zero in either the adj mat or its
    # transpose; we WLOG only consider positions on the lower triangle
    for (pos in which((state$mg+t(state$mg) != 0) & lower.tri(state$mg))) {

      # Get position corrdinates and transpose index
      pos.vec <- arrayInd(pos, c(p,p))
      pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

      # Delete edge
      newstate <- state
      newstate$mg[pos] <- 0
      newstate$mg[pos.trans] <- 0
      hash <- componentHash(list(cnodes=1:p), newstate$mg)
      if (! (hash %in% mgs.seen)) {
        mgs.seen[length(mgs.seen)+1] <- hash

        c1 <- state$node.comp[pos.vec[1]]

        # Check edge type
        if (state$mg[pos] == 100) {  # Bidirected edge

          # We might need to split up the component
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new components if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }
        } else {  # Directed edge

          c2 <- state$node.comp[pos.vec[2]]

          if (c1 == c2) {
            # Both edge-endpoints belong to the same component - update that
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          } else {
            # Edge-endpoints lie in different components
            if (state$mg[pos] == 1) {
              # c2 is the sink node component
              # Update parents of c2 component
              if (all(newstate$mg[pos.vec[1], newstate$comp[[c2]]$cnodes, drop=FALSE]==0)) newstate$comp[[c2]]$parents <- setdiff(newstate$comp[[c2]]$parents, pos.vec[1])
              # Update component name
              cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
              names(newstate$comp)[c2] <- cname2
              # Compute score if necessary
              if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
              # Set score of candidate model
              newstate$score <- state$score - scores[[names(state$comp)[c2]]] + scores[[cname2]]
            } else {
              # c1 is the sink node component
              # Update parents of c1 component
              if (all(newstate$mg[pos.vec[2], newstate$comp[[c1]]$cnodes, drop=FALSE]==0)) newstate$comp[[c1]]$parents <- setdiff(newstate$comp[[c1]]$parents, pos.vec[2])
              # Update component name
              cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
              names(newstate$comp)[c1] <- cname1
              # Compute score if necessary
              if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
              # Set score of candidate model
              newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
            }
          }
        }

        # Add new state to candidates
        i.d <- i.d + 1
        cand.del[[i.d]] <- newstate
      }
    }
  }

  # Edge type changes
  cand.cha <- list()
  i.c <- 0
  if ((direction == 3) || (direction == 4)) {

    # Loop over all '100's in the lower triangle of the adj mat
    for (pos in which(lower.tri(state$mg) & (state$mg == 100))) {

      # Get position corrdinates and transpose index
      pos.vec <- arrayInd(pos, c(p,p))
      pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

      # Check which component the edge lies in
      c1 <- state$node.comp[pos.vec[1]]

      # Check if 1 in lower triangle is admissible
      newstate <- state
      newstate$mg[pos] <- 1
      newstate$mg[pos.trans] <- 0
      tmp <- newstate$mg
      tmp[tmp==100] <- 0
      if (ggm::isAcyclic(tmp)) {
        hash <- componentHash(list(cnodes=1:p), newstate$mg)
        if (! (hash %in% mgs.seen)) {
          mgs.seen[length(mgs.seen)+1] <- hash
          # We might need to split up the component
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }

          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }
      }

      # Check if 1 in upper triangle is admissible
      newstate <- state
      newstate$mg[pos] <- 0
      newstate$mg[pos.trans] <- 1
      tmp <- newstate$mg
      tmp[tmp==100] <- 0
      if (ggm::isAcyclic(tmp)) {
        hash <- componentHash(list(cnodes=1:p), newstate$mg)
        if (! (hash %in% mgs.seen)) {
          mgs.seen[length(mgs.seen)+1] <- hash
          # Compute sub-components of changed components
          res <- connectedComponents(newstate$mg[state$comp[[c1]]$cnodes, state$comp[[c1]]$cnodes, drop=FALSE])
          if (length(res$comp) > 1) {
            # There is more than one -> split up
            # Express node indices in computed components in terms of original indices
            res$comp[[1]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[1]]$cnodes, drop=FALSE]
            res$comp[[2]]$cnodes <- state$comp[[c1]]$cnodes[res$comp[[2]]$cnodes, drop=FALSE]
            # Re-compute parents of new components
            res$comp[[1]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[1]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[1]]$cnodes))
            res$comp[[2]]$parents <- sort(setdiff(unique(unlist(sapply(res$comp[[2]]$cnodes, function(node) which(newstate$mg[,node]==1)))), res$comp[[2]]$cnodes))
            # Update component list
            newstate$comp[[c1]] <- res$comp[[1]]
            c2 <- length(newstate$comp) + 1
            newstate$comp[[c2]] <- res$comp[[2]]
            # Update node.comp: Assign nodes in new 2nd component correct index
            newstate$node.comp[res$comp[[2]]$cnodes] <- c2
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute scores of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]] + scores[[cname2]]
          } else {
            # There is still only one component
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          }

          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }
      }
    }

    # Loop over all '1's in adj mat
    for (pos in which(state$mg == 1)) {

      # Get position coordinates and transpose index
      pos.vec <- arrayInd(pos, c(p,p))
      pos.trans <- (pos.vec[1]-1)*p + pos.vec[2]

      c1 <- state$node.comp[pos.vec[1]]
      c2 <- state$node.comp[pos.vec[2]]

      # Check if 1 in transposed position is possible
      newstate <- state
      newstate$mg[pos] <- 0
      newstate$mg[pos.trans] <- 1
      tmp <- newstate$mg
      tmp[tmp==100] <- 0
      if (ggm::isAcyclic(tmp)) {
        hash <- componentHash(list(cnodes=1:p), newstate$mg)
        if (! (hash %in% mgs.seen)) {
          mgs.seen[length(mgs.seen)+1] <- hash
          if (c1 == c2) {
            # Both edge-endpoints belong to the same component - update that
            # Update component name
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            # Compute score of new component if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] + scores[[cname1]]
          } else {
            # Edge-endpoints lie in different components
            # c1 is the component of the old source node (new sink node)
            # c2 is the component of the old sink node (new source node)
            # Update parents of c1 component: add node
            newstate$comp[[c1]]$parents <- sort(unique(c(newstate$comp[[c1]]$parents, pos.vec[2])))
            # Update parents of c2 component: remove node
            if (all(newstate$mg[pos.vec[1], newstate$comp[[c2]]$cnodes, drop=FALSE]==0)) newstate$comp[[c2]]$parents <- setdiff(newstate$comp[[c2]]$parents, pos.vec[1])
            # Update component names
            cname1 <- componentHash(newstate$comp[[c1]], newstate$mg)
            cname2 <- componentHash(newstate$comp[[c2]], newstate$mg)
            names(newstate$comp)[c1] <- cname1
            names(newstate$comp)[c2] <- cname2
            # Compute score if necessary
            if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            if (! cname2 %in% names(scores)) scores[[cname2]] <- computeComponentScore(newstate$comp[[c2]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
            # Set score of candidate model
            newstate$score <- state$score - scores[[names(state$comp)[c1]]] - scores[[names(state$comp)[c2]]] + scores[[cname1]] + scores[[cname2]]
          }

          i.c <- i.c + 1
          cand.cha[[i.c]] <- newstate
        }
      }

      # Add model with bidirected edge
      newstate <- state
      newstate$mg[pos] <- 100
      newstate$mg[pos.trans] <- 100
      hash <- componentHash(list(cnodes=1:p), newstate$mg)
      if (! (hash %in% mgs.seen)) {
        mgs.seen[length(mgs.seen)+1] <- hash
        # New component index: if components are joined, indices could shift
        c1.new <- c1 - as.numeric(c2 < c1)
        if (c1 != c2) {  # Need to join components
          # Join components (add nodes and parents of 2nd component to 1st one)
          newstate$comp[[c1]]$cnodes <- sort(c(state$comp[[c1]]$cnodes, state$comp[[c2]]$cnodes))
          newstate$comp[[c1]]$parents <- sort(setdiff(c(state$comp[[c1]]$parents, state$comp[[c2]]$parents), newstate$comp[[c1]]$cnodes))
          # Delete 2nd component (WARNING: this changes the indices of comp!)
          newstate$comp[[c2]] <- NULL
          # Update node.comp: shift all indices after c2 one down (since c2 was deleted)
          newstate$node.comp[newstate$node.comp > c2] <- newstate$node.comp[newstate$node.comp > c2] - 1
          # Update node.comp: change all indices of 2nd component to new index of 1st component
          newstate$node.comp[state$node.comp==c2] <- c1.new
        }
        # Update component name
        cname1 <- componentHash(newstate$comp[[c1.new]], newstate$mg)
        names(newstate$comp)[c1.new] <- cname1
        # Compute score of new component if necessary
        if (! cname1 %in% names(scores)) scores[[cname1]] <- computeComponentScore(newstate$comp[[c1.new]], newstate$mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
        newstate$score <- state$score - scores[[names(state$comp)[c1]]] - as.numeric(c1!=c2)*scores[[names(state$comp)[c2]]] + scores[[cname1]]

        i.c <- i.c + 1
        cand.cha[[i.c]] <- newstate
      }
    }
  }

  # Find neighbours, that lie in epsilon band around score
  neighbours <- list()
  if (i.a > 0) neighbours <- cand.add[which(sapply(1:i.a, function(j) abs(cand.add[[j]]$score-score) < epsilon))]
  if (i.d > 0) neighbours <- c(neighbours, cand.del[which(sapply(1:i.d, function(j) abs(cand.del[[j]]$score-score) < epsilon))])
  if (i.c > 0) neighbours <- c(neighbours, cand.cha[which(sapply(1:i.c, function(j) abs(cand.cha[[j]]$score-score) < epsilon))])

  cumtime <- cumtime + proc.time()[3] - t

  # Recursively search neighbourhood of each of the neighbours
  res <- list()
  index <- 0
  for (neighbour in neighbours) {
    index <- index + 1
    tmp <- findEquivalentModels(neighbour, mgs.seen, scores, score, epsilon,
                                depth=depth+1, depth.max=depth.max,
                                cumtime=cumtime, time.max=time.max, data=data,
                                n=n, direction=direction, maxIter=maxIter,
                                edge.penalty=edge.penalty, verbose=verbose,
                                covMat=covMat, faithful.eps=faithful.eps)
    cumtime <- tmp$cumtime
    res[[length(res)+1]] <- tmp$res
    mgs.seen <- tmp$mgs.seen
    scores <- tmp$scores
  }

  # # Include starting state
  # if (depth==1) {
  #   neighbours <- c(list(state), neighbours)
  # }

  # Combine results of sub searches with neighbours found here and return
  return(list(res=c(neighbours, unlist(res, recursive=FALSE)), scores=scores,
              mgs.seen=mgs.seen, cumtime=cumtime))
}


fastFindEquivalentModels <- function(state, mgs.seen, scores, score, epsilon, depth.max,
                                     data=NULL, n=NULL, maxIter=10, edge.penalty=1,
                                     verbose=TRUE, covMat=NULL, faithful.eps=0)
{

  t <- proc.time()[3]

  # Populate EC with Collider invariant models
  models <- getColliderInvariantModels(state$mg)

  # Put them into mgs.seen
  ec <- list()
  for (mg in models) {
    hash <- componentHash(list(cnodes=1:nrow(mg)), mg)
    if (! (hash %in% mgs.seen)) mgs.seen[length(mgs.seen)+1] <- hash

    # Create state
    comp <- connectedComponents(mg)
    newscore <- 0
    for (i in 1:length(comp$comp)) {
      comphash <- componentHash(comp$comp[[i]], mg)
      compscore <- computeComponentScore(comp$comp[[i]], mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps)
      if (! (comphash %in% names(scores))) scores[[comphash]] <- compscore
      newscore <- newscore + compscore
    }
    newscores <- sum(unlist(lapply(1:length(comp$comp), function(i) computeComponentScore(comp$comp[[i]], mg, covMat, data, n, maxIter, edge.penalty, faithful.eps=faithful.eps))))
    newstate <- list()
    newstate$mg <- mg
    newstate$comp <- comp$comp
    newstate$node.comp <- comp$node.comp
    newstate$score <- newscore
    ec[[length(ec)+1]] <- newstate
  }

  cumtime <- proc.time()[3] - t

  # Start greedy equivalence search from each of these models
  for (newstate in ec) {

    res <- findEquivalentModels(newstate, mgs.seen, scores, newstate$score, epsilon,
                                depth.max=depth.max, time.max=Inf, data=data,
                                n=n, direction=4, maxIter=maxIter,
                                edge.penalty=edge.penalty, verbose=verbose,
                                covMat=covMat, faithful.eps=faithful.eps)
    mgs.seen <- res$mgs.seen
    scores <- res$scores
    ec <- append(ec, res$res)
    cumtime <- cumtime + res$cumtime
  }


  return(list(res=ec, scores=scores, mgs.seen=mgs.seen, cumtime=cumtime))
}


causalEffects <- function(p, max.in.degree, Bdist, Oscale, n, pop.version, R,
                          equivalent.eps, maxIter, maxSteps, depth.max=p*(p-1)/2,
                          time.max=Inf, faithful.eps=0, verbose=TRUE, max.pos=Inf,
                          mc.cores=1, forward=TRUE, fast=FALSE)
{

  # Generate ground truth
  gt <- GenerateGroundTruth(p, max.in.degree, Bdist, Oscale, faithful.eps=faithful.eps)

  # Create data and covariance matrices
  if (pop.version) {
    data <- NULL
    covMat <- gt$covMat
  } else {
    data <- GenerateData(n, gt$params)
    covMat <- cov(data)
  }

  # Run R greedy searches
  res.greedy <- greedySearch(
    data,
    n.restarts = R,
    max.iter.ricf = maxIter,
    max.steps = maxSteps,
    max.in.degree = max.in.degree,
    verbose = verbose,
    mc.cores = mc.cores
  )

  # Find highest-scoring model
  i <- which.max(sapply(res.greedy$scores, function(scores) scores[length(scores)]))

  # Find equivalent models of greedy result
  if (fast){
    tmp <- fastFindEquivalentModels(list(mg=res.greedy$final.bap), c(), list(),
                                    res.greedy$final.score, equivalent.eps, depth.max=depth.max,
                                    data=data, n=n, maxIter=maxIter, covMat=covMat,
                                    faithful.eps=faithful.eps)
  } else {
    tmp <- findEquivalentModels(list(mg=res.greedy$final.bap), c(), list(),
                                res.greedy$final.score, equivalent.eps, depth.max=depth.max,
                                time.max=time.max, data=data, n=n, maxIter=maxIter,
                                covMat=covMat, faithful.eps=faithful.eps)
  }
  scores <- tmp$scores
  equiv.models.greedy <- tmp$res
  if (tmp$cumtime >= time.max) print("Max time during EC search exceeded!")
  print(paste("Number of models in emp. EC:", length(tmp$res)))

  # Find equivalent models of ground truth
  if (fast) {
    tmp <- fastFindEquivalentModels(list(mg=gt$mg), c(), scores, NA, equivalent.eps,
                                    depth.max=depth.max, data=data, n=n, maxIter=maxIter,
                                    covMat=covMat, faithful.eps=faithful.eps)
  } else {
    tmp <- findEquivalentModels(list(mg=gt$mg), c(), scores, NA, equivalent.eps,
                                depth.max=depth.max, time.max=time.max, data=data,
                                n=n, maxIter=maxIter, covMat=covMat, faithful.eps=faithful.eps)
  }
  equiv.models.gt <- tmp$res
  if (tmp$cumtime >= time.max) print("Max time during EC search exceeded!")
  print(paste("Number of models in true EC:", length(tmp$res)))

  # Compute causal effects of greedy result
  ### TODO: Wrap these up in try(...)
  CE.greedy <- lapply(equiv.models.greedy,
                      function(model) getCausalEffects(
                        fitAncestralGraphCustom(model$mg, covMat, n, maxIter=maxIter)$Bhat))
  CE.minabs.greedy <- do.call("pmin", lapply(CE.greedy, function(mat) abs(mat))) - diag(p)

  # Compute causal effects of ground truth
  CE.gt <- lapply(equiv.models.gt,
                  function(model) getCausalEffects(
                    fitAncestralGraphCustom(model$mg, covMat, n, maxIter=maxIter)$Bhat))
  CE.minabs.gt <- do.call("pmin", lapply(CE.gt, function(mat) abs(mat))) - diag(p)

  return(list(gt=gt, covMat=covMat, res.greedy=res.greedy, i=i, equiv.models.gt=equiv.models.gt,
              equiv.models.greedy=equiv.models.greedy, CE.gt=CE.gt, CE.greedy=CE.greedy,
              CE.minabs.gt=CE.minabs.gt, CE.minabs.greedy=CE.minabs.greedy, data=data,
              n=n))
}
