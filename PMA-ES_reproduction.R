###############################################################################
#  PMA-ES: Probabilistic Morphological Analysis with Entropy Screening
#  R reproduction of the analysis pipeline of the paper (Faye & Cardebat).
#  Single source of all figures, tables and figure series reported in the paper.
#
#  Computations in base R. Optional final Excel export (writexl or openxlsx).
#  The script reproduces:
#    A. Real case "global competitiveness of wine industries" (Faye, Cardebat &
#       Hocquard, 2023, Food Systems no. 8, horizon 2035): meso marginals,
#       entropy, disagreement between experts, correlation across experts (sign),
#       MDI strength, hybrid copula, and comparison of independence, all-positive
#       and hybrid dependence (concentration, leading scenario, Spearman), plus
#       Janus factors.
#    A bis. Scenario families under dependence and their robustness across
#       classification procedures (latent classes, adaptive tree, copula
#       ranking): Table 6 of the paper.
#    A ter. Sensitivity of the concentration to rho_max (Section 5.2).
#    A quater. Robustness to the choice of copula: Student copula (Section 5.2).
#    B. Synthetic case: simulation study (validation, Section 5.1), grid of
#       panel size E x coupling strength (pivots, sign, RMSE, false couplings)
#       and robustness to phi, following Morris, White & Crowther (2019).
#    B bis. Benchmark: distance of the reconstructed joint distribution to the
#       true one, PMA-ES versus independence and structure-only dependence
#       (Table 3 of the paper).
#  A final Excel export gathers all these results, including the series
#  plotted in Figures 1 to 3 (sheets Fig1_series, Fig2_series, Fig3_series).
#
#  Run:  source("PMA-ES_reproduction.R")   (about 30 to 60 minutes, mostly the
#  simulations in parts B and B bis. Lower REPS_SIM and REPS_BENCH for a quick
#  test. Part B bis uses two cores on Linux and macOS, one on Windows. Each
#  replication has its own seed, so results do not depend on the number of cores.)
#  Monte Carlo: N_MC = 1e6 copula draws, seed 1.
###############################################################################

## ===========================================================================
## 0. CORE FUNCTIONS
## ===========================================================================

# Shannon entropy normalized in base k: 0 = certainty, 1 = equiprobability.
norm_entropy <- function(p) {
  p <- p[p > 0]
  k <- length(p)
  if (k <= 1) return(0)
  -sum(p * log(p)) / log(k)
}

# Projection onto the nearest positive semidefinite correlation matrix
# (simplified Higham algorithm).
nearest_corr <- function(A, iter = 300) {
  Y <- A ; dS <- matrix(0, nrow(A), ncol(A))
  for (i in seq_len(iter)) {
    R <- Y - dS
    e <- eigen(R, symmetric = TRUE)
    w <- pmax(e$values, 0)
    X <- e$vectors %*% diag(w) %*% t(e$vectors)
    dS <- X - R ; Y <- X ; diag(Y) <- 1
  }
  (Y + t(Y)) / 2
}

# Returns a positive definite correlation matrix (strict floor on the
# eigenvalues, so that the Cholesky decomposition succeeds).
make_pd <- function(C, floor_ = 1e-6) {
  if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) >= floor_) return(C)
  C <- nearest_corr(C)
  e <- eigen(C, symmetric = TRUE)
  w <- pmax(e$values, floor_)
  C <- e$vectors %*% diag(w) %*% t(e$vectors)
  d <- sqrt(diag(C)); C <- C / outer(d, d)
  (C + t(C)) / 2
}

# Ordinal Gaussian copula: sample of scenarios (levels 0..k-1) that exactly
# preserves the marginals P (list of vectors) under the latent correlation C.
sample_copula <- function(P, C, N = 1e6, seed = 0) {
  set.seed(seed)
  n <- length(P)
  U <- chol(C)                                   # U'U = C (U upper triangular)
  Z <- matrix(rnorm(N * n), N, n) %*% U
  levels_ <- matrix(0L, N, n)
  for (i in seq_len(n)) {
    k <- length(P[[i]])
    thresholds <- qnorm(cumsum(P[[i]])[1:(k - 1)])   # k-1 latent thresholds
    levels_[, i] <- findInterval(Z[, i], thresholds)  # -> 0..k-1
  }
  levels_
}

# Exact joint mass under independence (product of marginals) over the whole space.
joint_indep <- function(P) {
  grid <- expand.grid(lapply(P, function(p) seq_along(p) - 1L))
  apply(grid, 1L, function(cell) prod(mapply(function(p, m) p[m + 1L], P, cell)))
}


# ---------------------------------------------------------------------------
# Families under dependence and their robustness across classification
# procedures: latent classes (LCA), adaptive partitioning tree
# (gain = H(v) + sum of MI(v, .)), and direct ranking of the joint distribution.
# ---------------------------------------------------------------------------

# Latent class analysis by EM, weighted by the mass of each cell.
lca_em <- function(X, w, K, Km, n_start = 5, iter = 500, tol = 1e-9, seed = 1) {
  set.seed(seed); n <- nrow(X); J <- ncol(X); best <- NULL
  for (s in seq_len(n_start)) {
    pi_ <- rep(1 / K, K)
    th  <- lapply(seq_len(J), function(j) { m <- matrix(runif(K * Km[j]), K); m / rowSums(m) })
    ll_old <- -Inf
    for (it in seq_len(iter)) {
      lp <- matrix(log(pi_), n, K, byrow = TRUE)
      for (j in seq_len(J)) lp <- lp + t(log(th[[j]][, X[, j] + 1L, drop = FALSE]))
      mx <- apply(lp, 1, max); pr <- exp(lp - mx); s_ <- rowSums(pr); post <- pr / s_
      ll <- sum(w * (log(s_) + mx))
      wp <- post * w; pi_ <- colSums(wp) / sum(w)
      for (j in seq_len(J)) {
        m <- matrix(vapply(0:(Km[j] - 1), function(v) colSums(wp[X[, j] == v, , drop = FALSE]),
                           numeric(K)), K)
        th[[j]] <- (m + 1e-10) / rowSums(m + 1e-10)
      }
      if (abs(ll - ll_old) < tol * abs(ll)) break
      ll_old <- ll
    }
    if (is.null(best) || ll > best$ll) best <- list(ll = ll, pi = pi_, theta = th, post = post)
  }
  best
}

# Families under dependence and robustness.
#   P       : list of marginals. Km : number of levels per factor
#   sample_ : sample of scenarios under the copula (N x n, levels 0..k-1)
#   p_ind   : exact masses under independence, in expand.grid order
#   labels  : display labels of the factors. symbols : list, one vector of
#             level symbols per factor
triangulation <- function(P, Km, sample_, p_ind, labels,
                          symbols = lapply(Km, function(k) if (k == 2) c("~", "+") else c("-", "~", "+")),
                          lca_threshold = 0.05, k_max = 8, struct_threshold = 0.60,
                          min_leaf_weight = 0.10, tree_depth = 3, copula_top = 50,
                          n_top_indep = 3, verbose = TRUE) {
  n <- length(P)
  # 1. Joint distribution under the copula, on the observed cells
  mult <- cumprod(c(1, head(Km, -1)))
  tab  <- table(as.integer(sample_ %*% mult))
  keys <- as.integer(names(tab)); w <- as.numeric(tab) / sum(tab)
  decode <- function(k) { v <- integer(n); for (i in seq_len(n)) { v[i] <- k %% Km[i]; k <- k %/% Km[i] }; v }
  X <- t(vapply(keys, decode, integer(n)))
  p_dep <- numeric(prod(Km)); p_dep[keys + 1L] <- w
  rank_dep <- rank(-p_dep, ties.method = "min")
  rank_ind <- rank(-p_ind, ties.method = "min")

  # 2. LCA: add a class as long as it captures >= lca_threshold of total dependence
  fits <- list(lca_em(X, w, 1, Km, n_start = 1))
  TC <- sum(w * log(w)) - fits[[1]]$ll          # total correlation (nats)
  K <- 1; trace <- data.frame(K = 1, share_captured = 0, gain = NA)
  if (TC > 1e-6) for (k in 2:k_max) {
    fits[[k]] <- lca_em(X, w, k, Km)
    gain <- (fits[[k]]$ll - fits[[k - 1]]$ll) / TC
    trace <- rbind(trace, data.frame(K = k, share_captured = (fits[[k]]$ll - fits[[1]]$ll) / TC, gain = gain))
    if (gain < lca_threshold) break
    K <- k
  }
  lca <- fits[[K]]
  cls <- max.col(lca$post, ties.method = "first")

  # 3. Adaptive tree. A family is recovered if it dominates (> 50%) a leaf
  H_w <- function(x, ww, k) { p <- vapply(0:(k - 1), function(v) sum(ww[x == v]), 0)
                              p <- p / sum(p); p <- p[p > 0]; -sum(p * log(p)) }
  split_gain <- function(idx, v) {
    ww <- w[idx]; g <- 0
    for (j in seq_len(n)) {
      Hc <- 0
      for (m in 0:(Km[v] - 1)) { s <- X[idx, v] == m
        if (any(s)) Hc <- Hc + sum(ww[s]) / sum(ww) * H_w(X[idx[s], j], ww[s], Km[j]) }
      g <- g + H_w(X[idx, j], ww, Km[j]) - Hc
    }
    g
  }
  leaves <- list()
  grow <- function(idx, depth, used) {
    cand <- setdiff(seq_len(n), used)
    if (depth >= tree_depth || !length(cand)) { leaves[[length(leaves) + 1]] <<- idx; return(invisible()) }
    v <- cand[which.max(vapply(cand, function(v) split_gain(idx, v), 0))]
    children <- lapply(0:(Km[v] - 1), function(m) idx[X[idx, v] == m])
    if (any(vapply(children, function(e) sum(w[e]), 0) < min_leaf_weight)) {
      leaves[[length(leaves) + 1]] <<- idx; return(invisible()) }
    for (e in children) grow(e, depth + 1, c(used, v))
  }
  grow(seq_along(w), 0, integer(0))
  dom <- vapply(leaves, function(f) {
    pc <- tapply(w[f], factor(cls[f], levels = seq_len(K)), sum); pc[is.na(pc)] <- 0
    if (max(pc) / sum(pc) > 0.5) which.max(pc) else NA_integer_ }, integer(1))

  # 4. Table of families (representative = most probable scenario of the class)
  marg_mode <- vapply(P, which.max, integer(1))
  ord <- order(-lca$pi)
  compo <- function(cell, cl = NULL) paste(vapply(seq_len(n), function(j) {
    s <- symbols[[j]][cell[j] + 1L]
    if (!is.null(cl) && lca$theta[[j]][cl, cell[j] + 1L] >= struct_threshold) s <- paste0("[", s, "]")
    if (cell[j] + 1L != marg_mode[j]) s <- paste0(s, "*")
    paste0(labels[j], s) }, ""), collapse = " ")
  fam <- do.call(rbind, lapply(seq_along(ord), function(r) {
    cl <- ord[r]; s <- which(cls == cl); i <- s[which.max(w[s])]; key <- keys[i] + 1L
    m2 <- any(dom == cl, na.rm = TRUE); m3 <- rank_dep[key] <= copula_top
    data.frame(Family = paste0("F", r), Weight_pct = round(100 * lca$pi[cl], 1),
      Composition = compo(X[i, ], cl),
      Mass_dep_1e3 = round(1000 * p_dep[key], 2), Rank_dep = rank_dep[key],
      Mass_ind_1e3 = round(1000 * p_ind[key], 2), Rank_ind = rank_ind[key],
      Robustness = 1L + m2 + m3,
      Procedures = paste(c("LCA", if (m2) "Tree", if (m3) paste0("Copula top ", copula_top)), collapse = " + "),
      stringsAsFactors = FALSE) }))
  top <- order(-p_ind)[seq_len(n_top_indep)]
  contrast <- do.call(rbind, lapply(top, function(i) data.frame(
    Family = paste0("Top_indep_", rank_ind[i]), Weight_pct = NA, Composition = compo(decode(i - 1L)),
    Mass_dep_1e3 = round(1000 * p_dep[i], 2), Rank_dep = rank_dep[i],
    Mass_ind_1e3 = round(1000 * p_ind[i], 2), Rank_ind = rank_ind[i],
    Robustness = NA_integer_, Procedures = "", stringsAsFactors = FALSE)))
  table_ <- rbind(fam, contrast)

  if (verbose) {
    cat(sprintf("  Total dependence = %.3f nats. %d families capture %.1f%% of the dependence. Tree: %d leaves\n",
                TC, K, 100 * trace$share_captured[trace$K == K], length(leaves)))
    cat("  [x] = structuring level (P >= 0.60 within the family). * = differs from the marginal mode\n")
    print(table_, row.names = FALSE, right = FALSE)
  }
  list(table = table_, K = K, lca_trace = trace, total_dependence = TC,
       n_leaves = length(leaves), p_dep = p_dep)
}


## ===========================================================================
## A. REAL CASE: GLOBAL COMPETITIVENESS OF WINE INDUSTRIES
##    (Faye, Cardebat & Hocquard, 2023, Food Systems no. 8, horizon 2035,
##     10 meso factors retained after influence-dependence selection)
## ===========================================================================
cat("\n=====================================================================\n")
cat("A. WINE COMPETITIVENESS: entropy, signs, MDI, hybrid copula\n")
cat("=====================================================================\n")

# Probabilities by expert: list of matrices (levels x 13 experts).
MESO <- list(
 Identity = rbind(c(0.3,0.2,0.4,0.4,0.6,0.6,0.4,0.3,0.4,0.4,0.5,0.3,0.6),
                  c(0.7,0.8,0.6,0.6,0.4,0.4,0.6,0.7,0.6,0.6,0.5,0.7,0.4)),
 CSR = rbind(c(0.1,0.6,0.2,0.2,0.1,0.1,0.1,0.1,0.2,0.3,0.1,0.2,0.2),
             c(0.3,0.4,0.5,0.2,0.3,0.3,0.2,0.3,0.5,0.4,0.3,0.2,0.3),
             c(0.6,0,0.3,0.6,0.6,0.6,0.7,0.6,0.3,0.3,0.6,0.6,0.5)),
 Financing = rbind(c(0.3,0.4,0.5,0.4,0.3,0.3,0.4,0.3,0.5,0.5,0.4,0.5,0.3),
                   c(0.7,0.6,0.5,0.6,0.7,0.7,0.6,0.7,0.5,0.5,0.6,0.5,0.7)),
 Support = rbind(c(0.7,0.5,0.4,0.4,0.5,0.5,0.3,0.3,0.4,0.6,0.6,0.6,0.6),
                 c(0.3,0.5,0.6,0.6,0.5,0.5,0.7,0.7,0.6,0.4,0.4,0.4,0.4)),
 Capacity = rbind(c(0.3,0.1,0.3,0.3,0.3,0.3,0.2,0.2,0.3,0.3,0.5,0.3,0.6),
                  c(0.6,0.6,0.5,0.5,0.5,0.5,0.7,0.2,0.5,0.5,0.4,0.4,0.3),
                  c(0.1,0.3,0.2,0.2,0.2,0.2,0.3,0.6,0.2,0.2,0.1,0.3,0.1)),
 Brand = rbind(c(0.1,0.6,0.2,0.3,0.3,0.3,0.1,0.2,0.2,0.4,0.5,0.3,0.1),
               c(0.7,0.2,0.4,0.5,0.6,0.6,0.5,0.6,0.4,0.3,0.2,0.4,0.4),
               c(0.2,0.2,0.4,0.2,0.1,0.1,0.4,0.2,0.4,0.3,0.3,0.3,0.5)),
 Agroeco = rbind(c(0.1,0.2,0.2,0.2,0.2,0.2,0.1,0.3,0.2,0.2,0.2,0.2,0.2),
                 c(0.7,0.6,0.6,0.3,0.4,0.4,0.3,0.4,0.6,0.5,0.5,0.3,0.4),
                 c(0.2,0.2,0.2,0.5,0.4,0.4,0.6,0.3,0.2,0.3,0.3,0.5,0.4)),
 Digital = rbind(c(0.1,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0),
                 c(0.6,0.3,0.4,0.2,0.3,0.3,0.2,0.1,0.4,0.4,0.3,0.3,0.2),
                 c(0.3,0.5,0.4,0.7,0.7,0.7,0.7,0.8,0.4,0.4,0.5,0.5,0.8)),
 R_D = rbind(c(0.1,0.5,0.1,0.2,0.4,0.4,0.1,0.1,0.1,0.2,0.2,0.3,0.1),
             c(0.6,0.4,0.5,0.5,0.4,0.4,0.6,0.4,0.5,0.4,0.5,0.3,0.3),
             c(0.3,0.1,0.4,0.3,0.2,0.2,0.4,0.5,0.4,0.4,0.3,0.4,0.6)),
 Tourism = rbind(c(0.1,0.1,0.3,0.3,0.5,0.5,0.2,0.1,0.3,0.2,0.5,0.2,0.3),
                 c(0.3,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.2,0.2),
                 c(0.6,0.6,0.5,0.5,0.3,0.3,0.6,0.7,0.5,0.5,0.2,0.6,0.5)))
facM <- names(MESO); nM <- length(facM)
labM <- c("Identity","CSR","Financing","Support","Capacity","Brand",
          "Agroeco","Digital","R&D","Tourism")   # display labels

# Normalization by expert -> marginals, entropy, disagreement (sample SD), leaning
NORM <- lapply(MESO, function(A) sweep(A, 2, colSums(A), "/"))
Pm   <- lapply(NORM, rowMeans)
Km   <- sapply(Pm, length)
Hm   <- sapply(Pm, norm_entropy)
Dm   <- sapply(NORM, function(A) mean(apply(A, 1, sd)) * 100)   # sd() = sample SD
lean <- t(sapply(NORM, function(A) {
  k <- nrow(A); ((seq_len(k) - 1)/(k - 1)) %*% A }))            # nM x 13
Rie  <- cor(t(lean))                                            # sign of dependence

cat(sprintf("\n%-12s %3s %8s %9s   distribution (%%)\n","Factor","k","entropy","disagr."))
for (i in seq_len(nM)) {
  dist <- paste(sprintf("%.0f", 100 * Pm[[i]]), collapse = " / ")
  cat(sprintf("%-12s %3d %8.2f %9.1f   %s\n", labM[i], Km[i], Hm[i], Dm[i], dist))
}

# Meso direct influence matrix (mean of experts, 0-3), influence of the ROW on
# the COLUMN, extracted from the full MDI (22 factors) of the source study.
MDI_meso <- rbind(
c(0,0.5,0.5,0.75,0.75,1.5,0.75,1.25,0.5,1.75),
c(2,0,1,1,1.25,1.75,2.25,1.5,1.5,1),
c(1.5,1,0,1.75,1.75,1.75,1,1,1,0.5),
c(2,1.25,0.75,0,1.25,1.5,1.5,0.75,1.25,1.25),
c(0.75,1.5,2.5,2.25,0,1.5,2.25,0,1.75,1.25),
c(2.75,2,1.5,1.25,1.75,0,2,1.75,1.5,1.75),
c(1.5,2.5,1.5,1.25,2,1.75,0,1,2.5,1.5),
c(1.5,1,2.25,1.25,1.25,2,0.75,0,2,1.75),
c(1.25,2,2.5,1.75,1.75,2.25,1.75,1.25,0,0.25),
c(2,1.75,1,1.25,1,1.5,1.25,1.75,0.75,0))
dimnames(MDI_meso) <- list(facM, facM)
# Coupling strength calibrated on the symmetrized meso MDI (rho_max = 0.6)
Bsym <- (MDI_meso + t(MDI_meso)) / 2 ; diag(Bsym) <- 0
RHO_MAX <- 0.6
MAG <- RHO_MAX * Bsym / max(Bsym)                 # positive strength

# --- HYBRID copula: MDI strength x sign across experts, 0 if |r| < 0.15 ---
ETA <- 0.15
Chyb <- matrix(0, nM, nM)
for (i in seq_len(nM)) for (j in seq_len(nM))
  if (i != j && abs(Rie[i, j]) >= ETA) Chyb[i, j] <- MAG[i, j] * sign(Rie[i, j])
diag(Chyb) <- 1
n_neg <- sum(Chyb[upper.tri(Chyb)] < 0)          # number of negative couplings (before projection)
Chyb <- make_pd(Chyb)

# All-positive reference
Cpos <- MAG ; diag(Cpos) <- 1 ; Cpos <- make_pd(Cpos)
cat(sprintf("\nHybrid copula: %d negative couplings out of 45 | min eigenvalue = %.3f\n",
            n_neg, min(eigen(Chyb, symmetric = TRUE, only.values = TRUE)$values)))

# --- Concentration: independence (exact) vs all-positive and hybrid (Monte Carlo) ---
N_MC <- 1e6
p_ind <- joint_indep(Pm)
smp_pos <- sample_copula(Pm, Cpos, N = N_MC, seed = 1)
smp_hyb <- sample_copula(Pm, Chyb, N = N_MC, seed = 1)

# joint mass (copula) by cell, via a table on a mixed-radix key
mixed_key <- function(sample_, K) {
  mult <- cumprod(c(1, head(K, -1)))
  as.integer(sample_ %*% mult)
}
p_from_sample <- function(sample_, K) {
  tab <- table(mixed_key(sample_, K)) / nrow(sample_)
  tab
}
top50 <- function(pv) sum(sort(pv, decreasing = TRUE)[1:50])
tp <- p_from_sample(smp_pos, Km); th <- p_from_sample(smp_hyb, Km)

cat(sprintf("\nCumulative mass of the top 50:\n"))
cat(sprintf("   independence         : %5.1f %%\n", 100 * top50(p_ind)))
cat(sprintf("   all-positive copula  : %5.1f %%\n", 100 * top50(as.numeric(tp))))
cat(sprintf("   hybrid copula        : %5.1f %%\n", 100 * top50(as.numeric(th))))
# Masses over the whole space (17,496 cells, expand.grid order) for Spearman
space_mass <- function(sample_, K) {
  v <- tabulate(mixed_key(sample_, K) + 1L, prod(K)); v / nrow(sample_) }
p_pos <- space_mass(smp_pos, Km); p_hyb <- space_mass(smp_hyb, Km)
spear <- function(a, b) cor(rank(a), rank(b))
cat(sprintf("\nLeading scenario (x1e-3): indep. %.1f | all-positive %.1f | hybrid %.1f\n",
            1000 * max(p_ind), 1000 * max(p_pos), 1000 * max(p_hyb)))
cat(sprintf("Spearman vs independence: all-positive %.2f | hybrid %.2f\n",
            spear(p_ind, p_pos), spear(p_ind, p_hyb)))

# --- Janus factors: strong MDI coupling but weak |r| across experts ---
cat("\nJanus factors (MDI strength vs correlation across experts):\n")
pairs_ <- list(c("Agroeco","R_D"), c("Brand","R_D"), c("Agroeco","CSR"), c("Brand","Agroeco"))
reading <- function(ar) {
  if (ar < ETA) "Janus (faces cancel out)"
  else if (ar < 0.35) "weak (not antagonistic)"
  else if (ar < 0.60) "comoving" else "strongly comoving"
}
for (pr in pairs_) {
  i <- match(pr[1], facM); j <- match(pr[2], facM)
  cat(sprintf("  %-18s MDI strength=%.2f  r=%+.2f  -> %s\n",
              paste(pr, collapse = "-"), MAG[i, j], Rie[i, j], reading(abs(Rie[i, j]))))
}

# Systematic rule (same as in the analyst tool): MDI strength >= rho_max/2 and
# |r| across experts < eta. The pairs above are a commented selection.
jj <- which(upper.tri(MAG) & MAG >= 0.5 * RHO_MAX & abs(Rie) < ETA, arr.ind = TRUE)
wine_Janus_rule <- data.frame(Factor_1 = labM[jj[, 1]], Factor_2 = labM[jj[, 2]],
  MDI_strength = round(MAG[jj], 2), r_across_experts = round(Rie[jj], 2), stringsAsFactors = FALSE)
wine_Janus_rule <- wine_Janus_rule[order(-wine_Janus_rule$MDI_strength, abs(wine_Janus_rule$r_across_experts)), ]
cat(sprintf("\nSystematic rule (MDI >= %.2f and |r| < %.2f): %d Janus pairs out of %d\n",
            0.5 * RHO_MAX, ETA, nrow(wine_Janus_rule), nM * (nM - 1) / 2))
print(wine_Janus_rule, row.names = FALSE)


## ===========================================================================
## A bis. SCENARIO FAMILIES UNDER DEPENDENCE (Table 6)
##   Families = latent classes of the joint distribution under the hybrid copula
##   (one more class as long as it captures >= 5% of total dependence).
##   Robustness = number of procedures that recover the family: LCA (by
##   construction), adaptive tree (the family dominates a leaf), copula top 50
##   (representative scenario).
## ===========================================================================
cat("\n=====================================================================\n")
cat("A bis. FAMILIES UNDER DEPENDENCE: LCA, adaptive tree, copula ranking\n")
cat("=====================================================================\n")
op <- options(width = 250)
fam_dep <- triangulation(Pm, Km, smp_hyb, p_ind, labM)
options(op)
wine_families_dep <- fam_dep$table

## ===========================================================================
## A ter. SENSITIVITY TO RHO_MAX (Section 5.2)
## ===========================================================================
cat("\nA ter. Sensitivity of the hybrid copula to rho_max:\n")
cell_txt <- function(i) {                          # composition of a cell (expand.grid order)
  k <- i - 1L; v <- integer(nM); for (j in seq_len(nM)) { v[j] <- k %% Km[j]; k <- k %/% Km[j] }
  paste(sapply(seq_len(nM), function(j) paste0(labM[j], (if (Km[j] == 2) c("~","+") else c("-","~","+"))[v[j] + 1L])), collapse = " ") }
wine_rho_sensitivity <- do.call(rbind, lapply(c(0.3, 0.45, 0.6, 0.75), function(rm) {
  MAGr <- rm * Bsym / max(Bsym); C <- matrix(0, nM, nM)
  for (i in seq_len(nM)) for (j in seq_len(nM))
    if (i != j && abs(Rie[i, j]) >= ETA) C[i, j] <- MAGr[i, j] * sign(Rie[i, j])
  diag(C) <- 1; C <- make_pd(C)
  p <- space_mass(sample_copula(Pm, C, N = N_MC, seed = 1), Km)
  data.frame(rho_max = rm, Top50_mass_pct = round(100 * top50(p), 1),
             Spearman_vs_indep = round(spear(p_ind, p), 2),
             Leading_scenario = cell_txt(which.max(p)), stringsAsFactors = FALSE) }))
print(wine_rho_sensitivity, row.names = FALSE, right = FALSE)


## ===========================================================================
## A quater. ROBUSTNESS TO THE CHOICE OF COPULA (Section 5.2)
##   Student copula on the same correlation matrix (4 and 8 degrees of
##   freedom). Latent thresholds are taken on the Student scale, so that the
##   marginals are preserved exactly.
## ===========================================================================
cat("\nA quater. Student copula versus Gaussian copula:\n")
sample_student <- function(P, C, df, N = 1e6, seed = 1) {
  set.seed(seed); n <- length(P); U <- chol(C)
  Z <- matrix(rnorm(N * n), N, n) %*% U
  Tm <- Z / sqrt(rchisq(N, df) / df)
  lv <- matrix(0L, N, n)
  for (i in seq_len(n)) { k <- length(P[[i]])
    lv[, i] <- findInterval(Tm[, i], qt(cumsum(P[[i]])[1:(k - 1)], df)) }
  lv }
wine_student <- list(); wine_student_families <- list()
for (df_ in c(Inf, 8, 4)) {
  smp <- if (is.infinite(df_)) smp_hyb else sample_student(Pm, Chyb, df_, N = N_MC, seed = 1)
  p <- space_mass(smp, Km)
  wine_student[[length(wine_student) + 1]] <- data.frame(
    Copula = if (is.infinite(df_)) "Gaussian" else paste0("Student df=", df_),
    Top50_mass_pct = round(100 * top50(p), 1), Leading_scenario_1e3 = round(1000 * max(p), 1),
    Spearman_vs_indep = round(spear(p_ind, p), 2), Spearman_vs_Gaussian = round(spear(p_hyb, p), 2),
    Same_leading_scenario = which.max(p) == which.max(p_hyb), stringsAsFactors = FALSE)
  if (!is.infinite(df_)) {
    cat(sprintf("\n  Families under the Student copula, df = %d:\n", df_))
    op <- options(width = 250)
    ft <- triangulation(Pm, Km, smp, p_ind, labM)$table
    options(op)
    wine_student_families[[length(wine_student_families) + 1]] <- cbind(Copula = paste0("Student df=", df_), ft)
  }
}
wine_student <- do.call(rbind, wine_student)
wine_student_families <- do.call(rbind, wine_student_families)
print(wine_student, row.names = FALSE)

cat("\n--- end of real case ---\n")


## ===========================================================================
## B. SIMULATION STUDY: SYNTHETIC CASE (validation, Section 5.1)
##    Following Morris, White & Crowther (2019). Recovery of pivots and of the
##    sign and strength of dependence by panel size E and coupling strength.
##    Robustness to violation of the covariation assumption (phi). Base R only.
##    REPS_SIM = 500 (as in the paper). Lower it to speed things up.
## ===========================================================================
cat("\n=====================================================================\n")
cat("B. SIMULATION STUDY: grid (E x strength) and robustness\n")
cat("=====================================================================\n")
set.seed(2026)
n_sim <- 10L; K_sim <- 3L; cc_sim <- c(-1, 0, 1); TEMP <- 0.55
REPS_SIM <- 500L

softmax_dist <- function(s) { e <- exp(-((s - cc_sim)^2) / TEMP); e / sum(e) }

make_truth <- function(strength) {                 # ground truth (marginals, C*, pivots)
  m <- numeric(n_sim); piv <- logical(n_sim)
  for (i in seq_len(n_sim)) {
    if (i <= n_sim %/% 2) m[i] <- sample(c(-1.2, 1.2), 1)   # concentrated core
    else { m[i] <- runif(1, -0.15, 0.15); piv[i] <- TRUE }  # flat pivot
  }
  C <- diag(n_sim); idx <- sample(n_sim); k <- 1L            # signed couplings in blocks
  while (k < n_sim) {
    i <- idx[k]; j <- idx[k + 1L]
    C[i, j] <- C[j, i] <- sample(c(-1, 1), 1) * strength; k <- k + 2L
  }
  e <- eigen(C, symmetric = TRUE); w <- pmax(e$values, 1e-6) # PSD projection
  C <- e$vectors %*% diag(w) %*% t(e$vectors)
  d <- sqrt(diag(C)); C <- C / outer(d, d)
  list(m = m, C = C, piv = piv)
}

one_rep <- function(E, strength, phi = 1, sigma = 0.9) {
  tr <- make_truth(strength); m <- tr$m; C <- tr$C; piv_true <- tr$piv
  U <- chol(C)                                     # U'U = C
  S <- matrix(m, E, n_sim, byrow = TRUE) +
       sigma * (phi * (matrix(rnorm(E * n_sim), E, n_sim) %*% U) +
                (1 - phi) * matrix(rnorm(E * n_sim), E, n_sim))
  ex <- array(0, c(E, n_sim, K_sim))               # vectorized softmax
  for (kk in seq_len(K_sim)) ex[, , kk] <- exp(-((S - cc_sim[kk])^2) / TEMP)
  den <- ex[, , 1] + ex[, , 2] + ex[, , 3]
  for (kk in seq_len(K_sim)) ex[, , kk] <- ex[, , kk] / den
  p_hat <- apply(ex, c(2, 3), mean)                # n x K
  w_lean <- (seq_len(K_sim) - 1) / (K_sim - 1)
  lean <- matrix(0, n_sim, E)
  for (kk in seq_len(K_sim)) lean <- lean + w_lean[kk] * t(ex[, , kk])
  Mt <- abs(C) * 3; diag(Mt) <- 0; Msum <- matrix(0, n_sim, n_sim)  # noisy true MDI
  for (e_ in seq_len(E))
    Msum <- Msum + pmax(pmin(round(Mt + matrix(rnorm(n_sim^2, 0, 0.6), n_sim, n_sim)), 3), 0)
  M <- Msum / E
  piv_hat <- apply(p_hat, 1, max) < 0.55
  Rie <- cor(t(lean)); B <- (M + t(M)) / 2; diag(B) <- 0; MAG <- 0.6 * B / max(B)
  eta_ad <- { t95 <- qt(0.975, E - 2); t95 / sqrt(E - 2 + t95^2) }   # 5% significance threshold
  Ch <- Cha <- matrix(0, n_sim, n_sim)
  for (i in seq_len(n_sim)) for (j in seq_len(n_sim)) if (i != j) {
    if (abs(Rie[i, j]) >= 0.15)   Ch[i, j]  <- MAG[i, j] * sign(Rie[i, j])   # fixed threshold
    if (abs(Rie[i, j]) >= eta_ad) Cha[i, j] <- MAG[i, j] * sign(Rie[i, j])   # adaptive threshold
  }
  up <- upper.tri(C); iu <- which(up, arr.ind = TRUE); v <- abs(C[up])
  coupled <- iu[v >= 0.15, , drop = FALSE]; nullp <- iu[v < 0.15, , drop = FALSE]
  c(piv      = mean(piv_hat == piv_true),
    sign     = if (nrow(coupled)) mean(sign(Ch[coupled])  == sign(C[coupled])) else NA,
    rmse     = if (nrow(coupled)) sqrt(mean((Ch[coupled] - C[coupled])^2))     else NA,
    false    = if (nrow(nullp))  mean(abs(Ch[nullp])  > 0) else 0,
    sign_ad  = if (nrow(coupled)) mean(sign(Cha[coupled]) == sign(C[coupled])) else NA,
    false_ad = if (nrow(nullp))  mean(abs(Cha[nullp]) > 0) else 0)
}

gl <- list()                                       # grid strength x E
for (strength in c(0.3, 0.5, 0.7)) for (E in c(14L, 30L, 50L)) {
  A <- replicate(REPS_SIM, one_rep(E, strength))
  mn <- rowMeans(A, na.rm = TRUE)
  se <- apply(A, 1L, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  gl[[length(gl) + 1L]] <- data.frame(
    strength = strength, E = E, pivots = round(mn["piv"], 3),
    sign = round(mn["sign"], 3), sign_SE = round(se["sign"], 3),
    RMSE = round(mn["rmse"], 3), false_coupl_fixed = round(mn["false"], 3),
    sign_adaptive = round(mn["sign_ad"], 3),
    false_coupl_adaptive = round(mn["false_ad"], 3), row.names = NULL)
}
synth_grid <- do.call(rbind, gl)
cat(sprintf("\nGrid (%d replications per condition):\n", REPS_SIM))
print(synth_grid, row.names = FALSE)

rl <- list()                                       # sign robustness vs phi (E=30, strength=0.5)
for (phi in c(1, 0.75, 0.5, 0.25, 0)) {
  A <- replicate(REPS_SIM, one_rep(30L, 0.5, phi = phi))
  rl[[length(rl) + 1L]] <- data.frame(phi = phi,
    sign = round(mean(A["sign", ], na.rm = TRUE), 3),
    sign_SE = round(sd(A["sign", ], na.rm = TRUE) / sqrt(REPS_SIM), 3), row.names = NULL)
}
synth_robustness <- do.call(rbind, rl)
cat("\nRobustness (sign vs phi, E=30, strength=0.5):\n"); print(synth_robustness, row.names = FALSE)


## ===========================================================================
## B bis. BENCHMARK: RECONSTRUCTION OF THE TRUE JOINT DISTRIBUTION (Table 3)
##   Same estimated marginals for every method, only the dependence differs.
##   Reference = copula with the true C*. Methods: independence, structure-only
##   copula (MDI strength, all signs positive), PMA-ES with the fixed threshold
##   (0.15) and with the 5% significance threshold. Criteria over the 3^10
##   scenarios: total variation distance and coverage efficiency of the fifty
##   leading scenarios. Common random numbers (1e5 draws per replication).
## ===========================================================================
cat("\n=====================================================================\n")
cat("B bis. BENCHMARK: reconstruction of the true joint distribution\n")
cat("=====================================================================\n")
REPS_BENCH <- 500L
mult_sim <- K_sim^(0:(n_sim - 1))
joint_from_Z <- function(Z, C, P) {
  X <- Z %*% chol(C); key <- integer(nrow(X))
  for (i in seq_len(n_sim)) key <- key + mult_sim[i] * findInterval(X[, i], qnorm(cumsum(P[i, ])[1:(K_sim - 1)]))
  tabulate(key + 1L, K_sim^n_sim) / nrow(X) }
bench_rep <- function(E, strength, phi = 1, sigma = 0.9, N = 1e5, top = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  tr <- make_truth(strength); m <- tr$m; C <- tr$C
  U <- chol(C)
  S <- matrix(m, E, n_sim, byrow = TRUE) + sigma * (phi * (matrix(rnorm(E * n_sim), E, n_sim) %*% U) +
       (1 - phi) * matrix(rnorm(E * n_sim), E, n_sim))
  ex <- array(0, c(E, n_sim, K_sim))
  for (kk in 1:K_sim) ex[, , kk] <- exp(-((S - cc_sim[kk])^2) / TEMP)
  den <- ex[, , 1] + ex[, , 2] + ex[, , 3]; for (kk in 1:K_sim) ex[, , kk] <- ex[, , kk] / den
  P <- apply(ex, c(2, 3), mean)
  lean <- matrix(0, n_sim, E); for (kk in 1:K_sim) lean <- lean + (kk - 1) / (K_sim - 1) * t(ex[, , kk])
  Mt <- abs(C) * 3; diag(Mt) <- 0; Msum <- matrix(0, n_sim, n_sim)
  for (e_ in 1:E) Msum <- Msum + pmax(pmin(round(Mt + matrix(rnorm(n_sim^2, 0, 0.6), n_sim, n_sim)), 3), 0)
  M <- Msum / E; R <- cor(t(lean)); B <- (M + t(M)) / 2; diag(B) <- 0; MAG <- 0.6 * B / max(B)
  eta_ad <- { t95 <- qt(0.975, E - 2); t95 / sqrt(E - 2 + t95^2) }
  Cpos <- MAG; diag(Cpos) <- 1
  Ch <- ifelse(abs(R) >= 0.15, MAG * sign(R), 0); diag(Ch) <- 1
  Ca <- ifelse(abs(R) >= eta_ad, MAG * sign(R), 0); diag(Ca) <- 1
  Z <- matrix(rnorm(N * n_sim), N, n_sim)
  pt <- joint_from_Z(Z, C, P)
  meth <- list(indep = diag(n_sim), allpos = make_pd(Cpos), hyb = make_pd(Ch), hyb_sig = make_pd(Ca))
  ttop <- order(-pt)[1:top]; best <- sum(pt[ttop])
  unlist(lapply(meth, function(Cm) {
    pm <- joint_from_Z(Z, Cm, P)
    mtop <- order(-pm)[1:top]
    c(tv = 0.5 * sum(abs(pm - pt)), overlap = mean(mtop %in% ttop), cover = sum(pt[mtop]) / best)
  })) }
run_reps <- function(f, n) {
  if (.Platform$OS.type == "unix") parallel::mclapply(seq_len(n), f, mc.cores = 2) else lapply(seq_len(n), f) }
bench_conds <- rbind(expand.grid(strength = c(0.3, 0.5, 0.7), E = c(14L, 30L, 50L), phi = 1),
                     data.frame(strength = 0.5, E = 30L, phi = c(0.5, 0)))
synth_benchmark <- do.call(rbind, lapply(seq_len(nrow(bench_conds)), function(r) {
  cd <- bench_conds[r, ]
  A <- do.call(cbind, run_reps(function(i) bench_rep(cd$E, cd$strength, cd$phi, seed = 2027L + 1000L * r + i), REPS_BENCH))
  mn <- rowMeans(A); se <- apply(A, 1, sd) / sqrt(REPS_BENCH)
  data.frame(cd, t(round(mn, 3)), max_MC_SE = round(max(se), 3), row.names = NULL) }))
names(synth_benchmark) <- sub("^allpos", "structure_only", sub("^hyb_sig", "PMAES_signif", sub("^hyb\\.", "PMAES_fixed.", names(synth_benchmark))))
cat(sprintf("\nBenchmark (%d replications per condition):\n", REPS_BENCH)); print(synth_benchmark, row.names = FALSE)


## ===========================================================================
## C. EXPORT OF RESULTS TO AN EXCEL WORKBOOK
##    writexl preferred (light), otherwise openxlsx, otherwise CSV (base R).
##    To install:  install.packages("writexl")   (or "openxlsx")
## ===========================================================================
output_xlsx <- "results_PMA-ES.xlsx"

# --- Real case (wine): marginals/entropy/disagreement, matrices, concentration, Janus ---
p1 <- sapply(seq_len(nM), function(i) if (Km[i] == 3) Pm[[i]][1] else NA_real_)
p2 <- sapply(seq_len(nM), function(i) Pm[[i]][Km[i] - 1])
p3 <- sapply(seq_len(nM), function(i) Pm[[i]][Km[i]])
wine_marginals <- data.frame(
  Factor = labM, k = Km,
  p_minus = round(p1, 3), p_tilde = round(p2, 3), p_plus = round(p3, 3),
  Entropy = round(Hm, 3), Disagreement_pts = round(Dm, 1),
  row.names = NULL, stringsAsFactors = FALSE)

mat_df <- function(M, labels) {           # matrix -> data.frame with a label column
  d <- as.data.frame(round(M, 3)); colnames(d) <- labels
  cbind(Factor = labels, d, stringsAsFactors = FALSE)
}
wine_MDI_meso      <- mat_df(MDI_meso, labM)
wine_MDI_sym       <- mat_df(Bsym,  labM)
wine_strength_MAG  <- mat_df(MAG,   labM)
wine_corr_experts  <- mat_df(Rie,   labM)
wine_hybrid_copula <- mat_df(Chyb,  labM)

wine_concentration <- data.frame(
  Configuration   = c("Independence", "All-positive copula", "Hybrid copula"),
  Top50_mass_pct  = round(100 * c(top50(p_ind), top50(as.numeric(tp)), top50(as.numeric(th))), 1),
  Leading_scenario_1e3 = round(1000 * c(max(p_ind), max(p_pos), max(p_hyb)), 1),
  Spearman_vs_indep = round(c(1, spear(p_ind, p_pos), spear(p_ind, p_hyb)), 2),
  stringsAsFactors = FALSE)

wine_Janus <- do.call(rbind, lapply(pairs_, function(pr) {
  i <- match(pr[1], facM); j <- match(pr[2], facM)
  data.frame(Pair = paste(pr, collapse = "-"),
             MDI_strength = round(MAG[i, j], 2), r_across_experts = round(Rie[i, j], 2),
             Reading = reading(abs(Rie[i, j])), stringsAsFactors = FALSE)
}))

# --- Series behind Figures 1 to 3 of the paper (synthetic case) ---
g <- synth_grid
fig1_series <- data.frame(E = c(14, 30, 50),
  Sign_rho_0.3 = g$sign[g$strength == 0.3], Sign_rho_0.5 = g$sign[g$strength == 0.5],
  Sign_rho_0.7 = g$sign[g$strength == 0.7], Pivots = g$pivots[g$strength == 0.5])
g5 <- g[g$strength == 0.5, ]
fig2_series <- data.frame(E = g5$E,
  Sign_fixed_0.15 = g5$sign, Sign_adaptive = g5$sign_adaptive,
  False_coupl_fixed = g5$false_coupl_fixed, False_coupl_adaptive = g5$false_coupl_adaptive)
fig3_series <- synth_robustness[, c("phi", "sign")]

sheets <- list(
  Wine_marginals      = wine_marginals,
  Wine_MDI_meso       = wine_MDI_meso,
  Wine_MDI_sym        = wine_MDI_sym,
  Wine_strength_MAG   = wine_strength_MAG,
  Wine_corr_experts   = wine_corr_experts,
  Wine_hybrid_copula  = wine_hybrid_copula,
  Wine_concentration  = wine_concentration,
  Wine_Janus          = wine_Janus,
  Wine_Janus_rule     = wine_Janus_rule,
  Wine_families_dep   = wine_families_dep,
  Wine_rho_sensitivity = wine_rho_sensitivity,
  Wine_student        = wine_student,
  Wine_student_families = wine_student_families,
  Synth_grid          = synth_grid,
  Synth_robustness    = synth_robustness,
  Synth_benchmark     = synth_benchmark,
  Fig1_series         = fig1_series,
  Fig2_series         = fig2_series,
  Fig3_series         = fig3_series)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(sheets, path = output_xlsx)
  cat(sprintf("-> Results written to %s (%d sheets)\n", output_xlsx, length(sheets)))
} else if (requireNamespace("openxlsx", quietly = TRUE)) {
  openxlsx::write.xlsx(sheets, file = output_xlsx, overwrite = TRUE)
  cat(sprintf("-> Results written to %s (%d sheets, via openxlsx)\n", output_xlsx, length(sheets)))
} else {
  folder <- "results_PMA-ES_csv"; dir.create(folder, showWarnings = FALSE)
  for (nm in names(sheets))
    write.csv(sheets[[nm]], file.path(folder, paste0(nm, ".csv")), row.names = FALSE)
  cat("-> Neither writexl nor openxlsx: results written as CSV files in",
      folder, "/\n   install writexl or openxlsx to get a real .xlsx file.\n")
}
