rm(list = ls())
library(ergm)
library(ergm.gwespmatch)

# =============================================================================
# gwespmatch on a real directed network: Sampson's monastery data
#
# samplike: 18 monks, directed "liking" network, 88 edges.
# Vertex attribute "group": Loyal (7), Outcasts (4), Turks (7).
#
# Sociological question: beyond general tie formation and homophily, do monks
# within the same faction form closed triangular structures?  The gwespmatch
# term lets us test this and choose which type of triangle is most relevant.
# =============================================================================

data(sampson)   # loads samplike

cat("Nodes:", network.size(samplike), "\n")
cat("Directed:", is.directed(samplike), "\n")
cat("Edges:", network.edgecount(samplike), "\n")
cat("Groups:\n"); print(table(samplike %v% "group")); cat("\n")

# =============================================================================
# 1.  How much within-faction triangulation is there, and what shape?
#
# OTP (k: tail->k->head): classic "shared liking path" — two monks who both
#   like a third faction-mate then also like each other.
# ITP (k: k->tail, head->k): reverse paths
# OSP (k: tail->k, head->k): both endpoints send ties to the same k
# ISP (k: k->tail, k->head): both endpoints receive ties from the same k
# RTP (k: tail<->k, head<->k): mutual ties to a shared faction-mate
# =============================================================================

cat("gwespmatch at decay=1.0 for each type:\n")
for (type in c("OTP", "ITP", "OSP", "ISP", "RTP")) {
  s <- as.numeric(summary(samplike ~ gwespmatch(1.0, type = type, match = "group")))
  cat(sprintf("  %-4s  %.3f\n", type, s))
}

# Compare with gwesp (all edges) and nodematch (dyadic homophily):
cat(sprintf("\n  gwesp(OTP, all edges):  %.3f\n",
            as.numeric(summary(samplike ~ gwesp(1.0, fixed = TRUE, type = "OTP")))))
cat(sprintf("  nodematch('group'):     %d\n\n",
            as.numeric(summary(samplike ~ nodematch("group")))))

# =============================================================================
# 2.  How does the OTP statistic vary with decay?
# =============================================================================

cat("OTP gwespmatch across decay values:\n")
for (decay in c(0.2, 0.5, 1.0, 1.5, 2.5)) {
  s <- as.numeric(summary(samplike ~ gwespmatch(decay, type = "OTP", match = "group")))
  cat(sprintf("  decay=%.1f:  %.3f\n", decay, s))
}

# =============================================================================
# 3.  ERGM fit
#
# Model: edges + nodematch + gwespmatch(OTP)
#
#   edges         -- baseline density
#   nodematch     -- dyadic same-faction homophily
#   gwespmatch    -- within-faction triangular closure (OTP)
#
# The gwespmatch coefficient tells us whether same-faction triads close
# *above and beyond* the homophily already captured by nodematch.
# =============================================================================

cat("\nFitting ERGM: edges + nodematch('group') + gwespmatch(OTP, decay=0.5) ...\n\n")

data(sampson)
data(faux.magnolia.high)

fit <- ergm(
  faux.magnolia.high ~ edges +
    nodematch("Race") +
    gwesp(0.1, fixed = T, type = "OTP") +
    gwespmatch(0.1, fixed = T, type = "OTP", match = "Race"),
  estimate = "MLE",
  control = control.ergm(seed = 1)
)

summary(fit)
