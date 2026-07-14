\name{gwespmatch}
\alias{gwespmatch}
\alias{ergm.gwespmatch}
\title{Geometrically Weighted Edgewise Shared Partners with Attribute Matching}
\description{
An ERGM term that counts geometrically weighted edgewise shared partners
(GWESP) restricted to edges whose two endpoints share the same value of a
specified vertex attribute.  Supports all five directed shared-partner types
and a curved-ERGM mode for jointly estimating the decay parameter.
}
\usage{
gwespmatch(decay = 0, match, fixed = TRUE, cutoff = 30, type = "OTP",
           homophily = "triad")
}
\arguments{
  \item{decay}{Non-negative GWESP decay parameter \eqn{\alpha}.  Controls
    how rapidly the weight of additional shared partners diminishes.  When
    \code{fixed = FALSE} this is used as an initial value for estimation.}

  \item{match}{Character string naming the vertex attribute.  Only edges
    where \code{attr[tail] == attr[head]} contribute to the statistic.}

  \item{homophily}{Character string controlling whether the shared partner
    \eqn{k} must also carry the same attribute value.
    \describe{
      \item{\code{"triad"} (default)}{\emph{All three} nodes must match: the
        focal edge \eqn{(t,h)} \emph{and} the closing node \eqn{k}.  Only
        homogeneous triads are counted.  This is the statistic of Grund &
        Densley (2015) (\code{GWESP_MATCH}) and the "categorical
        transitivity" of Hong et al. (2024).}
      \item{\code{"dyad"}}{Only the focal edge \eqn{(t,h)} must match; the
        closing node \eqn{k} is unconstrained.  A within-group tie counts as
        closed even when the third party is out-group.}
    }
    With a single nominal attribute these are the only two possibilities: if
    \eqn{x_t = x_h}, then \dQuote{\eqn{k} matches \eqn{t}} and \dQuote{\eqn{k}
    matches \eqn{h}} are the same condition.}

  \item{fixed}{Logical (default \code{TRUE}).
    \itemize{
      \item \code{TRUE}: returns a single GW-weighted statistic with the
        supplied \code{decay} value fixed.
      \item \code{FALSE}: curved ERGM mode.  Returns \code{cutoff} raw
        matched-ESP counts and jointly estimates the decay parameter via
        the standard GWDECAY map
        \eqn{\theta \cdot e^\alpha (1-(1-e^{-\alpha})^l)}{theta * exp(alpha) * (1-(1-exp(-alpha))^l)}.
    }}

  \item{cutoff}{Integer (default \code{30}).  Used only when
    \code{fixed = FALSE}: the number of ESP-level statistics returned
    (i.e.\ the histogram is truncated at \code{cutoff} shared partners).}

  \item{type}{Character string specifying the directed shared-partner type
    (ignored for undirected networks, which always use common neighbours).
    \describe{
      \item{\code{"OTP"} (default)}{Out-two-path: \eqn{k} such that
        \code{tail -> k -> head}.}
      \item{\code{"ITP"}}{In-two-path: \eqn{k} such that
        \code{k -> tail} and \code{head -> k}.}
      \item{\code{"RTP"}}{Reciprocated two-path: \eqn{k} mutually tied to
        both \code{tail} and \code{head}.}
      \item{\code{"OSP"}}{Outgoing shared partner: \eqn{k} such that
        \code{tail -> k} and \code{head -> k}.}
      \item{\code{"ISP"}}{Incoming shared partner: \eqn{k} such that
        \code{k -> tail} and \code{k -> head}.}
    }}
}
\details{
This term behaves like \code{\link[ergm]{gwesp}} but restricts the sum to
\emph{matched} edges, i.e.\ edges \eqn{(t,h)} (tail--head) with
\eqn{x_t = x_h}, where \eqn{x} is the \code{match} attribute.  Let
\eqn{\alpha} denote the decay, \eqn{q = 1 - e^{-\alpha}}, and let
\eqn{ESP(t,h)} be the number of shared partners of edge \eqn{(t,h)}
under the chosen \code{type} and \code{homophily} rule.

With \code{homophily = "triad"} (the default) this is exactly Grund &
Densley's (2015) \code{GWESP_MATCH} statistic, and -- at \code{decay = 0} --
Hong et al.'s (2024) categorical transitivity.  It is equivalent to computing
\code{\link[ergm]{gwesp}} on the subgraph induced by within-group edges, i.e.
\code{F(~gwesp(decay, fixed = TRUE), ~nodematch(attr))}.

\subsection{Statistic (\code{fixed = TRUE})}{
Writing \eqn{D_l} for the number of matched edges with exactly \eqn{l}
shared partners, the term contributes the single statistic
\deqn{g(y;\,\alpha) = e^{\alpha}\!\!\sum_{(t,h):\,x_t=x_h}\!\!
      \bigl(1 - q^{\mathrm{ESP}(t,h)}\bigr)
      = e^{\alpha}\sum_{l\ge 1}\bigl(1 - q^{l}\bigr) D_l.}{%
      g(y; alpha) = exp(alpha) * sum_{(t,h): x_t=x_h} (1 - q^ESP(t,h))
                  = exp(alpha) * sum_l (1 - q^l) D_l.}
A single matched edge with \eqn{L} shared partners thus contributes
\eqn{e^{\alpha}(1-q^{L})}, which equals 1 for \eqn{L=1} regardless of
\eqn{\alpha} and increases toward \eqn{L} as \eqn{\alpha \to \infty}.  Edges
with \eqn{ESP=0} contribute zero; at \eqn{\alpha=0} the statistic
collapses to \eqn{sum_l D_l}, the number of matched edges lying in at
least one triangle of the chosen type.
}

\subsection{Change statistic}{
Toggling a matched edge \eqn{(t,h)} affects \eqn{g} through the edge's own
term plus every existing matched edge \eqn{(a,b)} for which \eqn{t} or \eqn{h}
becomes (or ceases to be) a shared partner; each such edge contributes
\eqn{q^{L_{ab}-edgestate}}, where \eqn{L_{ab}} is its current shared-
partner count.  The neighbour enumeration is specific to the chosen directed
type.
}

\subsection{Curved form (\code{fixed = FALSE})}{
The model carries the raw matched-ESP distribution
\eqn{D_1,\dots,D_{cutoff}} as sufficient statistics and estimates the
decay jointly with the coefficient via the standard GWDECAY map
\deqn{\theta_l(\theta,\alpha) = \theta\, e^{\alpha}\bigl(1 - q^{l}\bigr),}{%
      theta_l(theta, alpha) = theta * exp(alpha) * (1 - q^l),}
with gradient
\eqn{d\theta_l/d\theta = e^{\alpha}(1-q^{l})} and
\eqn{d\theta_l/d\alpha = \theta[e^{\alpha}(1-q^{l}) - l q^{l-1}]}.
Summing \eqn{sum_l \theta_l D_l} recovers the \code{fixed = TRUE}
statistic, so the fixed term is the collapse of the curved term at a fixed
\eqn{\alpha}.
}

\subsection{Undirected networks}{
All five type codes are accepted but the statistic always uses common
neighbours (equivalent to the UTP type in \code{\link[ergm]{gwesp}}).
}
}
\references{
Grund, T. U. and Densley, J. A. (2015). Ethnic Homophily and Triad Closure:
Mapping Internal Gang Structure Using Exponential Random Graph Models.
\emph{Journal of Contemporary Criminal Justice} 31(3), 354--370.

Hong, C.-S., Paik, A., Ballakrishnen, S., Silver, C. and Boutcher, S. (2024).
Categorical closure: Transitivity and identities in longitudinal networks.
\emph{Social Networks} 79, 76--92.

Hunter, D. R. (2007). Curved exponential family models for social networks.
\emph{Social Networks} 29(2), 216--230.
}
\note{
\strong{RTP and ergm's shared-partner cache.}  \code{gwespmatch} implements
ergm's documented RTP rule (\eqn{k} is an RTP partner of \eqn{(i,j)} iff
\eqn{i \leftrightarrow k \leftrightarrow j}).  ergm's own shared-partner cache
miscomputes \emph{edgewise} RTP; \code{ergm}'s \code{dgwesp(type = "RTP")}
agrees with \code{gwespmatch} only when called with
\code{term.options = list(cache.sp = FALSE)}.  Bear this in mind when
comparing the two.
}
\author{Benjamin Rosche \email{benrosche@nyu.edu}}
\seealso{\code{\link[ergm]{gwesp}}, \code{\link[ergm]{nodematch}}}
\examples{
library(ergm)
library(ergm.gwespmatch)

## Undirected network ---------------------------------------------------------
data(faux.magnolia.high)
summary(faux.magnolia.high ~ gwespmatch(1.5, match = "Race"))

## Directed network (Sampson monastery) ---------------------------------------
data(sampson)

# Summary statistics for each directed type
for (type in c("OTP", "ITP", "OSP", "ISP", "RTP"))
  cat(type, as.numeric(summary(samplike ~
        gwespmatch(1.0, type = type, match = "group"))), "\n")

# ERGM: do monks form within-faction OTP triangles beyond simple homophily?
\dontrun{
fit <- ergm(samplike ~ edges
                     + nodematch("group")
                     + gwespmatch(0.5, type = "OTP", match = "group"),
            control = control.ergm(seed = 42))
summary(fit)
}
}
\keyword{models}
