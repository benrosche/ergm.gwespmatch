library(ergm)
library(network)

# ergm's shared-partner cache miscomputes edgewise RTP, so reference
# comparisons disable it.  See test "RTP matches ergm's documented definition".
TO <- list(cache.sp = FALSE)

rand_undirected <- function(n, p = .3, k = 3, seed) {
  set.seed(seed)
  g <- network(matrix(rbinom(n * n, 1, p), n, n), directed = FALSE)
  g %v% "grp" <- sample(seq_len(k), n, replace = TRUE)
  g
}
rand_directed <- function(n, p = .4, k = 2, seed) {
  set.seed(seed)
  m <- matrix(rbinom(n * n, 1, p), n, n); diag(m) <- 0
  m <- pmax(m, t(m) * rbinom(n * n, 1, .6)); diag(m) <- 0   # ensure mutual dyads
  g <- network(m, directed = TRUE)
  g %v% "grp" <- sample(seq_len(k), n, replace = TRUE)
  g
}

# ---------------------------------------------------------------------------
# Grund & Densley (2015), Table 1, Graph B.
# K4-minus-an-edge among four same-group nodes, plus an out-group node tied to
# two of them.  Published: GWESP = 7 / 7.27 / 7.55 ; GWESP_Match = 5 / 5.22 / 5.39
# ---------------------------------------------------------------------------
test_that("triad reproduces Grund & Densley (2015) Table 1, Graph B", {
  nw <- network.initialize(5, directed = FALSE)
  for (e in list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(1,5), c(2,5)))
    nw[e[1], e[2]] <- 1
  nw %v% "grp" <- c(1, 1, 1, 1, 2)

  # sanity: the graph really is theirs
  expect_equal(as.numeric(summary(nw ~ triangle)), 3)
  expect_equal(as.numeric(summary(nw ~ nodematch("grp"))), 5)
  expect_equal(as.numeric(summary(nw ~ esp(1:3))), c(6, 0, 1))

  expect_equal(round(as.numeric(summary(nw ~ gwespmatch(0.00, match = "grp"))), 2), 5.00)
  expect_equal(round(as.numeric(summary(nw ~ gwespmatch(0.25, match = "grp"))), 2), 5.22)
  expect_equal(round(as.numeric(summary(nw ~ gwespmatch(0.50, match = "grp"))), 2), 5.39)
})

# ---------------------------------------------------------------------------
# triad == gwesp computed on the within-group subgraph, which is exactly how
# Hong et al. (2024) implement categorical transitivity: F(~., ~nodematch()).
# ---------------------------------------------------------------------------
test_that("triad equals ergm F(~gwesp, ~nodematch) on undirected networks", {
  for (s in 1:5) {
    g <- rand_undirected(12, seed = s)
    for (d in c(0, 0.4, 1.0))
      expect_equal(
        as.numeric(summary(g ~ gwespmatch(d, match = "grp"), term.options = TO)),
        as.numeric(summary(g ~ F(~gwesp(d, fixed = TRUE), ~nodematch("grp")),
                           term.options = TO)),
        tolerance = 1e-10)
  }
})

test_that("triad equals ergm F(~dgwesp, ~nodematch) for all five directed types", {
  for (s in 1:3) {
    g <- rand_directed(10, seed = 100 + s)
    for (ty in c("OTP", "ITP", "RTP", "OSP", "ISP")) for (d in c(0, 0.5))
      expect_equal(
        as.numeric(summary(g ~ gwespmatch(d, type = ty, match = "grp"),
                           term.options = TO)),
        as.numeric(summary(g ~ F(~dgwesp(d, fixed = TRUE, type = ty), ~nodematch("grp")),
                           term.options = TO)),
        info = paste(ty, d), tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# REGRESSION: decay = 0 used to return the SUM of ESP counts instead of the
# number of matched edges with at least one shared partner, making the
# statistic discontinuous at 0.  Grund & Densley fix tau = 0 in their models,
# so this mattered.
# ---------------------------------------------------------------------------
test_that("decay = 0 is continuous (regression)", {
  nw <- network.initialize(5, directed = FALSE)
  for (e in list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(1,5), c(2,5)))
    nw[e[1], e[2]] <- 1
  nw %v% "grp" <- c(1, 1, 1, 1, 2)

  for (hom in c("triad", "dyad")) {
    at_zero <- as.numeric(summary(nw ~ gwespmatch(0,    match = "grp", homophily = hom)))
    near_zero <- as.numeric(summary(nw ~ gwespmatch(1e-7, match = "grp", homophily = hom)))
    expect_equal(at_zero, near_zero, tolerance = 1e-5, info = hom)
    expect_equal(at_zero, 5)   # 5 matched edges, each with >= 1 shared partner
  }
})

# ---------------------------------------------------------------------------
# REGRESSION: the change statistic must agree with the summary statistic for
# every toggle.  Two bugs used to break this:
#   (a) c_gwespmatch returned early when the FOCAL edge was unmatched.  Under
#       "dyad" that is wrong: a heterophilous edge can still close a matched
#       edge, because the closing node need not match.
#   (b) the RTP change statistic enumerated only the edges that gain `tail` as
#       a shared partner, never those that gain `head`.
# ---------------------------------------------------------------------------
test_that("change statistic agrees with summary statistic for every toggle", {
  for (dir in c(FALSE, TRUE)) {
    for (hom in c("triad", "dyad")) {
      types <- if (dir) c("OTP", "ITP", "RTP", "OSP", "ISP") else "OTP"
      for (ty in types) {
        g <- if (dir) rand_directed(11, seed = 7) else rand_undirected(11, k = 2, seed = 7)
        f <- if (dir) ~gwespmatch(0.4, type = ty, match = "grp", homophily = hom)
             else     ~gwespmatch(0.4, match = "grp", homophily = hom)
        set.seed(42)
        for (it in 1:15) {
          d <- sample(1:11, 2)
          newval <- 1L - as.integer(g[d[1], d[2]])
          gf <- ergm.godfather(update(f, g ~ .),
                               changes = cbind(d[1], d[2], newval),
                               stats.start = TRUE, changes.only = FALSE)
          g[d[1], d[2]] <- newval
          expect_equal(as.numeric(gf)[2],
                       as.numeric(summary(update(f, g ~ .))),
                       tolerance = 1e-10,
                       info = paste(if (dir) "directed" else "undirected", ty, hom, it))
        }
      }
    }
  }
})

test_that("curved (fixed = FALSE) change statistics agree with summary", {
  for (hom in c("triad", "dyad")) {
    g <- rand_undirected(11, k = 2, seed = 3)
    f <- ~gwespmatch(0.5, fixed = FALSE, cutoff = 5, match = "grp", homophily = hom)
    set.seed(9)
    for (it in 1:15) {
      d <- sample(1:11, 2)
      newval <- 1L - as.integer(g[d[1], d[2]])
      gf <- ergm.godfather(update(f, g ~ .),
                           changes = cbind(d[1], d[2], newval),
                           stats.start = TRUE, changes.only = FALSE)
      g[d[1], d[2]] <- newval
      expect_equal(as.numeric(gf[2, ]),
                   as.numeric(summary(update(f, g ~ .))),
                   tolerance = 1e-10, info = paste(hom, it))
    }
  }
})

# ---------------------------------------------------------------------------
# gwespmatch implements ergm's DOCUMENTED RTP definition ("k is an RTP shared
# partner of (i,j) iff i<->k<->j").  ergm's own shared-partner CACHE gets this
# wrong; with cache.sp = FALSE ergm agrees with us.  This test pins our answer
# to the definition, not to the buggy cache.
# ---------------------------------------------------------------------------
test_that("RTP matches ergm's documented definition (ergm's SP cache does not)", {
  # 1->2, 2->1, 2->3, 3->1, 3->2, 4->1, 4->3
  # mutual dyads: 1<->2 and 2<->3.  Edge 3->1 has RTP partner k=2 (3<->2<->1).
  A <- matrix(c(0,1,0,0,
                1,0,1,0,
                1,1,0,0,
                1,0,1,0), 4, 4, byrow = TRUE)
  g <- network(A, directed = TRUE)
  g %v% "grp" <- rep(1, 4)

  expect_equal(as.numeric(summary(g ~ mutual)), 2)

  # the correct answer, by the documented definition:
  expect_equal(as.numeric(summary(g ~ gwespmatch(0, type = "RTP", match = "grp"))), 1)

  # ergm agrees once its shared-partner cache is switched off ...
  expect_equal(
    as.numeric(summary(g ~ dgwesp(0, fixed = TRUE, type = "RTP"),
                       term.options = list(cache.sp = FALSE))),
    1)
  # ... and disagrees with itself when the cache is on (documents the ergm bug).
  expect_equal(
    as.numeric(summary(g ~ dgwesp(0, fixed = TRUE, type = "RTP"),
                       term.options = list(cache.sp = TRUE))),
    0)
})
