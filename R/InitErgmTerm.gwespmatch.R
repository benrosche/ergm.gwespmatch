#' Shared-partner cache auxiliary for gwespmatch (internal).
#'
#' Maintains shared-partner counts on the within-group subgraph when
#' homophily = "triad", or on the full network when homophily = "dyad".
#' This is what lets c_gwespmatch look a count up in O(1) rather than
#' recomputing it in O(degree).  Not intended to be used directly.
InitErgmTerm..gwespmatch.spcache <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw, arglist,
                      varnames      = c("type", "match", "homophily"),
                      vartypes      = c("character", "character", "character"),
                      defaultvalues = list("OTP", NULL, "triad"),
                      required      = c(FALSE, TRUE, FALSE))

  node_attr <- get.vertex.attribute(nw, a$match)
  if (!is.numeric(node_attr)) node_attr <- as.numeric(as.factor(node_attr))

  match_k <- if (match.arg(tolower(a$homophily), c("triad", "dyad")) == "triad") 1 else 0

  type_code <- if (is.directed(nw))
                 c(OTP = 1, ITP = 2, RTP = 3, OSP = 4, ISP = 5)[toupper(a$type)]
               else 0L

  list(
    name       = "_gwespmatch_spcache",
    coef.names = character(0),
    inputs     = c(type_code, match_k, node_attr),
    dependence = TRUE,
    pkgname    = "ergm.gwespmatch",
    soname     = "ergm.gwespmatch"
  )
}

InitErgmTerm.gwespmatch <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw, arglist,
                      varnames      = c("decay", "fixed", "cutoff", "type", "match", "homophily"),
                      vartypes      = c("numeric", "logical", "numeric", "character", "character", "character"),
                      defaultvalues = list(0, TRUE, 30, "OTP", NULL, "triad"),
                      required      = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))

  match_attr <- a$match
  node_attr  <- get.vertex.attribute(nw, match_attr)
  if (!is.numeric(node_attr))
    node_attr <- as.numeric(as.factor(node_attr))

  # homophily: does the shared partner k also have to match?
  #   "triad" (default): all three nodes share the attribute value.
  #                      Grund & Densley (2015); Hong et al. (2024).
  #   "dyad"           : only the focal edge (tail,head) matches; k is free.
  homophily <- match.arg(tolower(a$homophily), c("triad", "dyad"))
  match_k   <- if (homophily == "triad") 1 else 0

  # Directed type code (ergm convention: OTP=1 ITP=2 RTP=3 OSP=4 ISP=5).
  # Undirected networks always use common neighbours; type_code = 0.
  if (is.directed(nw)) {
    type_map  <- c(OTP = 1, ITP = 2, RTP = 3, OSP = 4, ISP = 5)
    type_code <- type_map[toupper(a$type)]
    if (is.na(type_code))
      stop(paste0("gwespmatch: invalid type '", a$type,
                  "'. Choose from OTP, ITP, RTP, OSP, ISP."))
    sp_type <- toupper(a$type)
  } else {
    type_code <- 0L
    sp_type   <- "UTP"
  }

  # INPUT_PARAM layout shared by both C functions:
  #   [decay_or_cutoff, type_code, match_k, attr_1, ..., attr_N]
  base_inputs <- c(type_code, match_k, node_attr)

  # Shared-partner cache, so the change statistics are O(degree) per toggle
  # rather than O(degree^2).  Built on the within-group subgraph for "triad"
  # and on the full network for "dyad" -- see gwespmatch_spcache.c.
  aux <- substitute(~ .gwespmatch.spcache(type = TY, match = MA, homophily = HO),
                    list(TY = sp_type, MA = match_attr, HO = homophily))
  aux <- as.formula(aux)

  if (a$fixed) {
    # -----------------------------------------------------------------------
    # fixed = TRUE: single GW-weighted statistic.
    # C function: s_gwespmatch / c_gwespmatch
    # -----------------------------------------------------------------------
    inputs    <- c(a$decay, base_inputs)
    coef.name <- paste("gwespmatch", match_attr, a$decay, sp_type, homophily,
                       sep = ".")

    list(
      name        = "gwespmatch",
      coef.names  = coef.name,
      inputs      = inputs,
      auxiliaries = aux,
      pkgname     = "ergm.gwespmatch",
      soname      = "ergm.gwespmatch",
      dependence  = TRUE
    )

  } else {
    # -----------------------------------------------------------------------
    # fixed = FALSE: curved ERGM — estimate decay jointly with the coefficient.
    # C function: s_gwespmatch_dist / c_gwespmatch_dist
    # Returns `cutoff` statistics: esp_l = # matched edges with exactly l
    # shared partners (l = 1..cutoff).
    # -----------------------------------------------------------------------
    cutoff <- as.integer(a$cutoff)
    inputs <- c(cutoff, base_inputs)

    coef.names <- paste0("esp.", sp_type, ".", homophily, "#", seq_len(cutoff))

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
    term_name <- paste("gwespmatch", match_attr, homophily, sep = ".")
    params    <- setNames(list(NULL,    a$decay),
                          c(term_name, paste0(term_name, ".decay")))

    list(
      name        = "gwespmatch_dist",
      coef.names  = coef.names,
      inputs      = inputs,
      params      = params,
      map         = gw_map,
      gradient    = gw_gradient,
      minpar      = c(-Inf, 0),
      auxiliaries = aux,
      pkgname     = "ergm.gwespmatch",
      soname      = "ergm.gwespmatch",
      dependence  = TRUE
    )
  }
}
