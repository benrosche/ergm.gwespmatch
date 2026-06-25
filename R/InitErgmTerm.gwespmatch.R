InitErgmTerm.gwespmatch <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw, arglist,
                      varnames      = c("decay", "fixed", "cutoff", "type", "match"),
                      vartypes      = c("numeric", "logical", "numeric", "character", "character"),
                      defaultvalues = list(0, TRUE, 30, "OTP", NULL),
                      required      = c(TRUE, FALSE, FALSE, FALSE, TRUE))

  match_attr <- a$match
  node_attr  <- get.vertex.attribute(nw, match_attr)
  if (!is.numeric(node_attr))
    node_attr <- as.numeric(as.factor(node_attr))

  # Directed type code (ergm convention: OTP=1 ITP=2 RTP=3 OSP=4 ISP=5).
  # Undirected networks always use common neighbours; type_code = 0.
  if (is.directed(nw)) {
    type_map  <- c(OTP = 1, ITP = 2, RTP = 3, OSP = 4, ISP = 5)
    type_code <- type_map[toupper(a$type)]
    if (is.na(type_code))
      stop(paste0("gwespmatch: invalid type '", a$type,
                  "'. Choose from OTP, ITP, RTP, OSP, ISP."))
  } else {
    type_code <- 0L
  }

  # INPUT_PARAM layout shared by both C functions:
  #   [decay_or_cutoff, type_code, attr_1, ..., attr_N]
  base_inputs <- c(type_code, node_attr)

  if (a$fixed) {
    # -----------------------------------------------------------------------
    # fixed = TRUE: single GW-weighted statistic.
    # C function: s_gwespmatch / c_gwespmatch
    # -----------------------------------------------------------------------
    inputs    <- c(a$decay, base_inputs)
    coef.name <- paste("gwespmatch", match_attr, a$decay,
                       if (is.directed(nw)) toupper(a$type) else "UTP",
                       sep = ".")

    list(
      name       = "gwespmatch",
      coef.names = coef.name,
      inputs     = inputs,
      pkgname    = "ergm.gwespmatch",
      soname     = "ergm.gwespmatch",
      dependence = TRUE
    )

  } else {
    # -----------------------------------------------------------------------
    # fixed = FALSE: curved ERGM — estimate decay jointly with the coefficient.
    # C function: s_gwespmatch_dist / c_gwespmatch_dist
    # Returns `cutoff` statistics: esp_l = # matched edges with exactly l
    # shared partners (l = 1..cutoff).
    # -----------------------------------------------------------------------
    cutoff  <- as.integer(a$cutoff)
    inputs  <- c(cutoff, base_inputs)
    sp_type <- if (is.directed(nw)) toupper(a$type) else "UTP"

    coef.names <- paste0("esp.", sp_type, "#", seq_len(cutoff))

    # Curved ERGM map:  theta_l = theta * exp(alpha) * (1 - (1-exp(-alpha))^l)
    # Parameters: x[1] = theta (freely estimated), x[2] = alpha = decay (>= 0)
    gw_map <- function(x, n, ...) {
      i     <- seq_len(n)
      alpha <- x[2]
      x[1] * exp(alpha) * (1 - (1 - exp(-alpha))^i)
    }

    # Jacobian rows: d(theta_l)/d(theta), d(theta_l)/d(alpha)
    gw_gradient <- function(x, n, ...) {
      i     <- seq_len(n)
      alpha <- x[2]
      q     <- 1 - exp(-alpha)
      w     <- exp(alpha) * (1 - q^i)         # d/d(theta)
      dw    <- w - i * q^(i - 1)              # exp(alpha)*(1-q^i) - i*q^(i-1)
      rbind(w, x[1] * dw)
    }

    # params: NULL = freely estimated; numeric = initial value for decay
    term_name <- paste("gwespmatch", match_attr, sep = ".")
    params    <- setNames(list(NULL,    a$decay),
                          c(term_name, paste0(term_name, ".decay")))

    list(
      name       = "gwespmatch_dist",
      coef.names = coef.names,
      inputs     = inputs,
      params     = params,
      map        = gw_map,
      gradient   = gw_gradient,
      minpar     = c(-Inf, 0),
      pkgname    = "ergm.gwespmatch",
      soname     = "ergm.gwespmatch",
      dependence = TRUE
    )
  }
}
