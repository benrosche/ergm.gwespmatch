library(ergm)
library(network)

# ---------------------------------------------------------------------------
# Helpers: build small networks with a "grp" vertex attribute.
# ---------------------------------------------------------------------------
make_net <- function(n, edges, grp) {
  net <- network.initialize(n, directed = TRUE)
  for (i in seq_len(nrow(edges))) add.edge(net, edges[i, 1], edges[i, 2])
  set.vertex.attribute(net, "grp", grp)
  net
}

make_net_u <- function(n, edges, grp) {
  net <- network.initialize(n, directed = FALSE)
  for (i in seq_len(nrow(edges))) add.edge(net, edges[i, 1], edges[i, 2])
  set.vertex.attribute(net, "grp", grp)
  net
}

# ---------------------------------------------------------------------------
# 1. Smoke test: term runs on faux.magnolia.high with Race (undirected)
# ---------------------------------------------------------------------------
test_that("gwespmatch runs on faux.magnolia.high with Race", {
  data(faux.magnolia.high, package = "ergm")
  s <- summary(faux.magnolia.high ~ gwespmatch(1.5, match = "Race"))
  expect_length(s, 1)
  expect_false(is.na(s))
  expect_true(is.finite(s))
  expect_gt(s, 0)
})

# ---------------------------------------------------------------------------
# 2. Stat is 0 when matched edges have no shared partners
#    Directed:   1->2 (same group), no two-path
#    Undirected: {1,2} (same group), no common neighbour
# ---------------------------------------------------------------------------
test_that("gwespmatch is 0 when no matched edge has a shared partner", {
  # Directed
  net_d <- make_net(3,
                    edges = matrix(c(1, 2), nrow = 1),
                    grp   = c(1, 1, 1))
  expect_equal(as.numeric(summary(net_d ~ gwespmatch(1.5, match = "grp"))), 0)

  # Undirected
  net_u <- make_net_u(3,
                      edges = matrix(c(1, 2), nrow = 1),
                      grp   = c(1, 1, 1))
  expect_equal(as.numeric(summary(net_u ~ gwespmatch(1.5, match = "grp"))), 0)
})

# ---------------------------------------------------------------------------
# 3. Known analytical value: minimal OTP triangle (directed), all same group
#
#    Edges: 1->2, 1->3, 3->2   (all nodes group=1)
#    Edge (1,2): 1 OTP shared partner (k=3: 1->3 AND 3->2)
#      g(1) = exp(a) * (1 - q^1) = 1.0  for any a
#    Edges (1,3) and (3,2): 0 OTP partners -> contribute 0
#    => gwespmatch = 1.0 regardless of decay value
# ---------------------------------------------------------------------------
test_that("gwespmatch equals 1 for minimal directed OTP triangle (any decay)", {
  net <- make_net(3,
                  edges = matrix(c(1, 2,
                                   1, 3,
                                   3, 2), nrow = 3, byrow = TRUE),
                  grp = c(1, 1, 1))

  for (decay in c(0.5, 1.0, 1.5, 2.5)) {
    s <- as.numeric(summary(net ~ gwespmatch(decay, match = "grp")))
    expect_equal(s, 1.0, tolerance = 1e-10,
                 label = paste0("decay = ", decay))
  }
})

# ---------------------------------------------------------------------------
# 4. Non-matching edges do not contribute their own shared-partner counts
#
#    Triangle with groups (2,1,1).
#
#    Directed (1->2, 1->3, 3->2):
#      (1,2) and (1,3) are cross-group -> skip.
#      (3,2) is matched but L2=0 (no k with 3->k AND k->2) -> 0.
#      => gwespmatch = 0
#
#    Undirected ({1,2}, {1,3}, {2,3}):
#      {1,2} and {1,3} are cross-group -> skip (no contribution from them).
#      {2,3} is matched; common neighbour = 1, L2=1, g(1)=1.0.
#      => gwespmatch = 1.0  (not 3.0, which would be the all-same-group value)
# ---------------------------------------------------------------------------
test_that("gwespmatch ignores edges between nodes of different groups", {
  edges_tri <- matrix(c(1, 2,
                        1, 3,
                        3, 2), nrow = 3, byrow = TRUE)
  grp_mixed <- c(2, 1, 1)

  # Directed: only matched edge has no OTP partners -> 0
  net_d <- make_net(3, edges = edges_tri, grp = grp_mixed)
  expect_equal(as.numeric(summary(net_d ~ gwespmatch(1.5, match = "grp"))),
               0, tolerance = 1e-10)

  # Undirected: only {2,3} is matched; node 1 is a shared partner -> 1.0
  net_u <- make_net_u(3,
                      edges = matrix(c(1, 2, 1, 3, 2, 3), nrow = 3, byrow = TRUE),
                      grp   = grp_mixed)
  expect_equal(as.numeric(summary(net_u ~ gwespmatch(1.5, match = "grp"))),
               1.0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Two separate matched triangles
#
#    Nodes 1-3 (group=1) and 4-6 (group=2); no cross-group edges.
#    Each triangle contributes g(1) = 1.0 per OTP edge / per triangle edge.
#
#    Directed (1->2, 1->3, 3->2, 4->5, 4->6, 6->5):
#      One OTP edge per triangle => 2 * 1.0 = 2.0
#
#    Undirected ({1,2},{1,3},{2,3}, {4,5},{4,6},{5,6}):
#      Each of the 6 edges has 1 shared partner => 6 * 1.0 = 6.0
# ---------------------------------------------------------------------------
test_that("gwespmatch sums contributions from separate matched triangles", {
  edges_d <- matrix(c(1, 2,
                      1, 3,
                      3, 2,
                      4, 5,
                      4, 6,
                      6, 5), nrow = 6, byrow = TRUE)
  edges_u <- matrix(c(1, 2,
                      1, 3,
                      2, 3,
                      4, 5,
                      4, 6,
                      5, 6), nrow = 6, byrow = TRUE)
  grp <- c(1, 1, 1, 2, 2, 2)

  net_d <- make_net(6,   edges = edges_d, grp = grp)
  net_u <- make_net_u(6, edges = edges_u, grp = grp)

  expect_equal(as.numeric(summary(net_d ~ gwespmatch(1.5, match = "grp"))),
               2.0, tolerance = 1e-10, label = "directed")
  expect_equal(as.numeric(summary(net_u ~ gwespmatch(1.5, match = "grp"))),
               6.0, tolerance = 1e-10, label = "undirected")
})

# ---------------------------------------------------------------------------
# 6. Equivalence with gwesp when all nodes share the same attribute
#
#    Directed:   6-node two-triangle network; gwespmatch should equal
#                gwesp(type="OTP") when every node has the same group.
#
#    Undirected: faux.magnolia.high with all-same Race; gwespmatch should
#                equal gwesp(type="OTP") (= shared-neighbour triangles).
# ---------------------------------------------------------------------------
test_that("gwespmatch equals gwesp(OTP) when all nodes share one attribute", {
  # Directed
  net_d <- make_net(6,
                    edges = matrix(c(1, 2,
                                     1, 3,
                                     3, 2,
                                     4, 5,
                                     4, 6,
                                     6, 5), nrow = 6, byrow = TRUE),
                    grp = rep(1L, 6))
  s_d_match <- as.numeric(summary(net_d ~ gwespmatch(1.5, match = "grp")))
  s_d_gwesp  <- as.numeric(summary(net_d ~ gwesp(1.5, fixed = TRUE, type = "OTP")))
  expect_equal(s_d_match, s_d_gwesp, tolerance = 1e-8, label = "directed")

  # Undirected
  data(faux.magnolia.high, package = "ergm")
  net_u <- faux.magnolia.high
  set.vertex.attribute(net_u, "Race", rep(1L, network.size(net_u)))
  s_u_match <- as.numeric(summary(net_u ~ gwespmatch(1.5, match = "Race")))
  s_u_gwesp  <- as.numeric(summary(net_u ~ gwesp(1.5, fixed = TRUE, type = "OTP")))
  expect_equal(s_u_match, s_u_gwesp, tolerance = 1e-8, label = "undirected")
})

# ---------------------------------------------------------------------------
# 7. Change-stat consistency: C_CHANGESTAT agrees with S_CHANGESTAT
#
#    ergm.godfather() applies a single edge toggle using the C_ change-stat.
#    The running total is compared with a fresh summary() call to verify that
#    c_gwespmatch and s_gwespmatch are mutually consistent.
#
#    Directed:   6-node two-triangle network; remove edge (1->2).
#                Initial = 2.0, after = 1.0 (one OTP triangle lost).
#
#    Undirected: faux.magnolia.high; remove the first same-Race edge that
#                changes the statistic (found dynamically).
# ---------------------------------------------------------------------------
test_that("c_gwespmatch change stat is consistent with s_gwespmatch", {

  # --- Directed -------------------------------------------------------
  net_d <- make_net(6,
                    edges = matrix(c(1, 2,
                                     1, 3,
                                     3, 2,
                                     4, 5,
                                     4, 6,
                                     6, 5), nrow = 6, byrow = TRUE),
                    grp = c(1, 1, 1, 2, 2, 2))
  s_d_init <- as.numeric(summary(net_d ~ gwespmatch(1.5, match = "grp")))
  net_d_after        <- net_d; net_d_after[1, 2] <- FALSE
  s_d_after          <- as.numeric(summary(net_d_after ~ gwespmatch(1.5, match = "grp")))

  gf_d     <- ergm.godfather(net_d ~ gwespmatch(1.5, match = "grp"),
                              changes      = matrix(c(1L, 2L, 0L), nrow = 1),
                              stats.start  = TRUE,
                              changes.only = FALSE)
  gf_d_mat <- as.matrix(gf_d)

  expect_equal(as.numeric(gf_d_mat[1, 1]), s_d_init,  tolerance = 1e-8,
               label = "directed: initial stat")
  expect_equal(as.numeric(gf_d_mat[2, 1]), s_d_after, tolerance = 1e-8,
               label = "directed: post-toggle stat")

  # --- Undirected (faux.magnolia.high) ---------------------------------
  data(faux.magnolia.high, package = "ergm")
  net_u  <- faux.magnolia.high
  race   <- net_u %v% "Race"
  el     <- as.edgelist(net_u)

  s_u_init <- as.numeric(summary(net_u ~ gwespmatch(1.5, match = "Race")))

  same_idx <- which(race[el[, 1]] == race[el[, 2]])
  e_test   <- NULL
  s_u_after <- NULL
  for (k in same_idx) {
    net_tmp <- net_u; net_tmp[el[k, 1], el[k, 2]] <- FALSE
    s_tmp <- as.numeric(summary(net_tmp ~ gwespmatch(1.5, match = "Race")))
    if (abs(s_tmp - s_u_init) > 1e-10) {
      e_test    <- el[k, ]
      s_u_after <- s_tmp
      break
    }
  }
  stopifnot(!is.null(e_test))

  gf_u     <- ergm.godfather(net_u ~ gwespmatch(1.5, match = "Race"),
                              changes      = matrix(c(e_test[1], e_test[2], 0L), nrow = 1),
                              stats.start  = TRUE,
                              changes.only = FALSE)
  gf_u_mat <- as.matrix(gf_u)

  expect_equal(as.numeric(gf_u_mat[1, 1]), s_u_init,  tolerance = 1e-8,
               label = "undirected: initial stat")
  expect_equal(as.numeric(gf_u_mat[2, 1]), s_u_after, tolerance = 1e-8,
               label = "undirected: post-toggle stat")
})

# ---------------------------------------------------------------------------
# 8. Undirected triangle = 3.0 (analytical check, any decay)
#
#    Three-node undirected triangle {1,2}, {1,3}, {2,3}, all same group.
#    Each of the 3 edges has 1 shared partner; g(1) = 1.0 for any decay.
#    => gwespmatch = 3.0
# ---------------------------------------------------------------------------
test_that("gwespmatch equals 3 for undirected triangle (any decay)", {
  net_u <- network.initialize(3, directed = FALSE)
  add.edge(net_u, 1, 2)
  add.edge(net_u, 1, 3)
  add.edge(net_u, 2, 3)
  set.vertex.attribute(net_u, "grp", c(1, 1, 1))

  for (decay in c(0.5, 1.0, 1.5, 2.5)) {
    s <- as.numeric(summary(net_u ~ gwespmatch(decay, match = "grp")))
    expect_equal(s, 3.0, tolerance = 1e-10,
                 label = paste0("decay = ", decay))
  }
})

# ---------------------------------------------------------------------------
# 9. ISP type: directed star 1<-2, 1<-3, 2->3, all same group
#
#    ISP(a,b) = # k with k->a AND k->b.
#    Edges in the network: 2->1, 3->1, 2->3.
#
#    Edge (2->1): ISP partners = k with k->2 AND k->1. No such k -> 0.
#    Edge (3->1): ISP partners = k with k->3 AND k->1. k=2: 2->3 AND 2->1. => 1.
#    Edge (2->3): ISP partners = k with k->2 AND k->3. No such k -> 0.
#
#    Matched edge with ISP>0: only (3->1). g(1) = 1.0.
#    => gwespmatch(type="ISP") = 1.0 for any decay.
# ---------------------------------------------------------------------------
test_that("gwespmatch ISP type gives correct value", {
  net <- make_net(3,
                  edges = matrix(c(2, 1,
                                   3, 1,
                                   2, 3), nrow = 3, byrow = TRUE),
                  grp   = c(1, 1, 1))

  for (decay in c(0.5, 1.0, 1.5)) {
    s <- as.numeric(summary(net ~ gwespmatch(decay, type = "ISP", match = "grp")))
    expect_equal(s, 1.0, tolerance = 1e-10,
                 label = paste0("ISP decay=", decay))
  }
})

# ---------------------------------------------------------------------------
# 10. OSP type: directed network, verify OSP count
#
#     Edges: 1->3, 2->3, 1->2  (nodes 1,2,3 all group=1)
#
#     OSP(a,b) = # k with a->k AND b->k.
#     Edge (1->2): OSP partners = k with 1->k AND 2->k. k=3: 1->3 AND 2->3. => 1.
#     Edge (1->3): OSP partners = k with 1->k AND 3->k. k=2? 1->2 but 3->2? No. => 0.
#     Edge (2->3): OSP partners = k with 2->k AND 3->k. None. => 0.
#
#     => gwespmatch(type="OSP") = g(1) = 1.0 for any decay.
# ---------------------------------------------------------------------------
test_that("gwespmatch OSP type gives correct value", {
  net <- make_net(3,
                  edges = matrix(c(1, 3,
                                   2, 3,
                                   1, 2), nrow = 3, byrow = TRUE),
                  grp   = c(1, 1, 1))

  for (decay in c(0.5, 1.0, 1.5)) {
    s <- as.numeric(summary(net ~ gwespmatch(decay, type = "OSP", match = "grp")))
    expect_equal(s, 1.0, tolerance = 1e-10,
                 label = paste0("OSP decay=", decay))
  }
})

# ---------------------------------------------------------------------------
# 11. fixed = FALSE: ESP distribution (directed OTP triangle)
#
#     Directed OTP triangle (nodes 1-3, all group=1):
#       1->2, 1->3, 3->2
#     Matched edges with OTP>0:  only (1->2) with 1 OTP partner.
#     => esp#1 = 1, esp#2 = esp#3 = ... = 0
#
#     The weighted sum with theta=1:
#       theta * g(1; alpha) = 1 * 1.0 = 1.0  (for any alpha)
#     which equals the fixed=TRUE result.
# ---------------------------------------------------------------------------
test_that("gwespmatch fixed=FALSE returns correct ESP distribution (directed)", {
  net <- make_net(3,
                  edges = matrix(c(1, 2,
                                   1, 3,
                                   3, 2), nrow = 3, byrow = TRUE),
                  grp = c(1, 1, 1))

  dist <- as.numeric(summary(net ~ gwespmatch(1.5, fixed = FALSE, cutoff = 5,
                                              match = "grp")))
  expect_length(dist, 5)
  expect_equal(dist[1], 1.0, tolerance = 1e-10, label = "esp#1 = 1")
  expect_equal(dist[2], 0.0, tolerance = 1e-10, label = "esp#2 = 0")
  expect_equal(dist[3], 0.0, tolerance = 1e-10, label = "esp#3 = 0")
})

# ---------------------------------------------------------------------------
# 12. fixed = FALSE: ESP distribution (undirected triangle)
#
#     Undirected triangle {1,2}, {1,3}, {2,3}, all group=1.
#     Each edge has 1 common neighbour.
#     => esp#1 = 3, esp#2 = ... = 0
# ---------------------------------------------------------------------------
test_that("gwespmatch fixed=FALSE returns correct ESP distribution (undirected)", {
  net_u <- network.initialize(3, directed = FALSE)
  add.edge(net_u, 1, 2); add.edge(net_u, 1, 3); add.edge(net_u, 2, 3)
  set.vertex.attribute(net_u, "grp", c(1, 1, 1))

  dist <- as.numeric(summary(net_u ~ gwespmatch(1.0, fixed = FALSE, cutoff = 5,
                                                match = "grp")))
  expect_length(dist, 5)
  expect_equal(dist[1], 3.0, tolerance = 1e-10, label = "esp#1 = 3")
  expect_equal(dist[2], 0.0, tolerance = 1e-10, label = "esp#2 = 0")
})

# ---------------------------------------------------------------------------
# 13. fixed=FALSE weighted sum equals fixed=TRUE result
#
#     For two matched OTP triangles (test 5, directed), fixed=TRUE gives 2.0.
#     With fixed=FALSE, applying the GW map with theta=1 and the right alpha
#     should recover the same value:
#       sum_l [ theta * g(l; alpha) * dist_l ] = theta * gwespmatch(alpha)
# ---------------------------------------------------------------------------
test_that("gwespmatch fixed=FALSE weighted sum matches fixed=TRUE", {
  edges_d <- matrix(c(1, 2, 1, 3, 3, 2, 4, 5, 4, 6, 6, 5), nrow = 6, byrow = TRUE)
  grp     <- c(1, 1, 1, 2, 2, 2)
  net_d   <- make_net(6, edges = edges_d, grp = grp)

  alpha  <- 1.5
  cutoff <- 10
  dist   <- as.numeric(summary(net_d ~
               gwespmatch(alpha, fixed = FALSE, cutoff = cutoff, match = "grp")))

  # GW weights for each level l: g(l; alpha) = exp(alpha)*(1-(1-exp(-alpha))^l)
  l   <- seq_len(cutoff)
  q   <- 1 - exp(-alpha)
  gw  <- exp(alpha) * (1 - q^l)

  weighted_sum <- sum(gw * dist)   # theta = 1

  s_fixed <- as.numeric(summary(net_d ~ gwespmatch(alpha, match = "grp")))
  expect_equal(weighted_sum, s_fixed, tolerance = 1e-8,
               label = "weighted sum (fixed=FALSE) == fixed=TRUE")
})
