# ergm.gwespmatch

An [**ergm**](https://github.com/statnet/ergm) extension package providing the
`gwespmatch` term: **geometrically weighted edgewise shared partners restricted
to within-group edges**.

`gwespmatch` is Hunter's (2007) `gwesp` term with two modifications:

1. **Matching restriction.** Only edges whose two endpoints share the same value
   of a nominal vertex attribute contribute to the statistic. This isolates
   triadic closure that happens *within* groups (same school, same race, …) from
   closure that crosses group boundaries.
2. **Directed shared-partner typing.** For directed networks you choose which
   shared-partner shape (OTP/ITP/RTP/OSP/ISP) defines a shared partner, following
   the `ergm` convention.

The core change and summary statistics are implemented in C for speed.

## Installation

```r
# install.packages("remotes")
remotes::install_github("benrosche/ergm.gwespmatch")
```

Or from a local clone:

```
R CMD INSTALL ergm.gwespmatch
```

Requires `ergm` and `network`.

## Usage

```r
gwespmatch(decay, match, fixed = TRUE, cutoff = 30, type = "OTP")
```

| Argument | Default | Description |
|----------|---------|-------------|
| `decay`  | — (required) | Decay parameter α ≥ 0. Fixed value if `fixed = TRUE`; initial value if `fixed = FALSE`. |
| `match`  | — (required) | Name of the vertex attribute. Only edges with `attr[tail] == attr[head]` contribute. |
| `fixed`  | `TRUE` | `TRUE`: single GW-weighted statistic. `FALSE`: curved ERGM, decay estimated jointly. |
| `cutoff` | `30` | Number of ESP bins when `fixed = FALSE`. |
| `type`   | `"OTP"` | Directed shared-partner type: `OTP`, `ITP`, `RTP`, `OSP`, `ISP`. Ignored for undirected networks (common neighbours). |

### Shared-partner types

For an edge (tail → head), node *k* is a shared partner under:

| Type | Condition on *k* |
|------|------------------|
| **OTP** | tail → *k* and *k* → head (out two-path) |
| **ITP** | *k* → tail and head → *k* (in two-path) |
| **RTP** | tail ↔ *k* and head ↔ *k* (reciprocated) |
| **OSP** | tail → *k* and head → *k* (outgoing shared partner) |
| **ISP** | *k* → tail and *k* → head (incoming shared partner) |

Undirected networks always use common neighbours regardless of `type`.

## Examples

```r
library(ergm)
library(ergm.gwespmatch)

## Undirected network --------------------------------------------------
data(faux.magnolia.high)
summary(faux.magnolia.high ~ gwespmatch(1.5, match = "Race"))

## Directed network (Sampson monastery) --------------------------------
data(sampson)

# Summary statistic for each directed type
for (type in c("OTP", "ITP", "OSP", "ISP", "RTP"))
  cat(type, as.numeric(summary(samplike ~
        gwespmatch(1.0, type = type, match = "group"))), "\n")

# Do monks form within-faction OTP triangles beyond simple homophily?
fit <- ergm(samplike ~ edges
                     + nodematch("group")
                     + gwespmatch(0.5, type = "OTP", match = "group"),
            control = control.ergm(seed = 42))
summary(fit)
```

## The statistic

With decay α fixed (`fixed = TRUE`), writing *D*<sub>ℓ</sub> for the number of
matched edges with exactly ℓ shared partners and *q* = 1 − e<sup>−α</sup>:

```
g(y; α) = e^α · Σ_ℓ (1 − q^ℓ) · D_ℓ
```

A matched edge with one shared partner contributes 1 regardless of α; additional
shared partners contribute with geometrically diminishing weight. With
`fixed = FALSE`, the term becomes a curved ERGM that carries the raw matched-ESP
distribution and estimates α jointly via the standard GWDECAY map
θ<sub>ℓ</sub> = θ·e<sup>α</sup>(1 − q<sup>ℓ</sup>).

See `?gwespmatch` for the full mathematical details.

## Author

Benjamin Rosche (<benrosche@nyu.edu>)

## License

GPL-3
