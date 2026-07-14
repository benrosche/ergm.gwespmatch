# ergm.gwespmatch

An [**ergm**](https://github.com/statnet/ergm) extension providing `gwespmatch`:
**geometrically weighted edgewise shared partners, restricted to within-group triads.**

It lets you ask a question that plain ERGM terms cannot: *is triadic closure
stronger **within** groups than between them?*

---

## The problem it solves

Suppose you find that same-race people are disproportionately tied to each
other. The obvious explanation is **homophily** — a preference for similar
others. But there is a second, purely structural explanation: **triadic
closure**. Friends-of-friends become friends, and if your friends are already
mostly same-race, closure alone will manufacture same-race ties without anyone
preferring anything.

So homophily and closure are confounded. The standard fix is to put both
`nodematch()` and `gwesp()` in the model and let them compete. That's necessary
but not sufficient, because it assumes closure works *at the same rate*
everywhere. If closure is actually *stronger inside groups* than across them,
that interaction is a real mechanism — and leaving it out biases both the
homophily and the closure estimates.

`gwespmatch` is the missing interaction term. Grund & Densley (2015) call it
`GWESP_MATCH`; Hong et al. (2024) call it **categorical closure** and show that
adding it improves model fit and changes the homophily coefficients.

You cannot get this by multiplying `nodematch` by `gwesp`. As Grund & Densley put
it: "One cannot simply multiply the scores for GWESP and homophily. Instead, one
must incorporate a new configuration in the model."

## What it counts

For each edge whose two endpoints share the attribute, count that edge's shared
partners, then weight additional partners with geometrically diminishing returns
(exactly as `gwesp` does, so the term doesn't degenerate).

The one design choice is whether the *closing node* must match too:

| `homophily` | Who must share the attribute | Meaning |
|---|---|---|
| **`"triad"`** (default) | the edge's endpoints **and** the shared partner | All three nodes homogeneous. Grund & Densley's `GWESP_MATCH`; Hong et al.'s categorical transitivity. |
| `"dyad"` | only the edge's endpoints | A within-group tie counts as closed even when the third party is out-group. |

```
      homophily = "triad"              homophily = "dyad"
                                                                
            ● k                              ○ k          ● = same group
           / \                              / \           ○ = other group
          ●───●   counted                  ●───●   counted
          t   h                            t   h
                                     
          ○ k                              (both are counted under "dyad";
         / \                                only the top one under "triad")
        ●───●     NOT counted
        t   h
```

The default reproduces the published statistic. Use `"dyad"` only if you
deliberately want the looser version — and say so, because it is **not** what
either paper estimates.

## Installation

```r
# install.packages("remotes")
remotes::install_github("benrosche/ergm.gwespmatch")
```

Requires `ergm` and `network`.

## Usage

```r
gwespmatch(decay, match, fixed = TRUE, cutoff = 30, type = "OTP",
           homophily = "triad")
```

| Argument | Default | Description |
|----------|---------|-------------|
| `decay` | — (required) | Decay α ≥ 0. Higher α = more credit for *additional* shared partners. `decay = 0` counts each matched edge that closes at least one triangle. Fixed if `fixed = TRUE`; a starting value if `fixed = FALSE`. |
| `match` | — (required) | Name of the vertex attribute defining the groups. |
| `homophily` | `"triad"` | `"triad"` = all three nodes match (the published statistic). `"dyad"` = only the edge's endpoints match. |
| `fixed` | `TRUE` | `TRUE`: one statistic at the given `decay`. `FALSE`: curved ERGM, `decay` estimated jointly with the coefficient. |
| `cutoff` | `30` | Number of ESP bins when `fixed = FALSE`. |
| `type` | `"OTP"` | Directed shared-partner shape (below). Ignored for undirected networks. |

### Worked example

```r
library(ergm)
library(ergm.gwespmatch)

data(faux.magnolia.high)   # undirected friendship network, with a Race attribute

fit <- ergm(faux.magnolia.high ~ edges
                               + nodematch("Race")                    # dyadic homophily
                               + gwesp(0.5, fixed = TRUE)             # closure in general
                               + gwespmatch(0.5, match = "Race"))     # EXTRA closure within race
summary(fit)
```

Read the coefficients as a hierarchy:

- `nodematch("Race")` — do same-race ties form more often at the **dyad** level?
- `gwesp` — does closure happen at all, regardless of race?
- `gwespmatch` — **on top of both**, is closure *extra* likely when all three
  people are the same race?

A positive, significant `gwespmatch` is the categorical-closure finding. Expect
the `nodematch` coefficient to shrink once you add it — that shrinkage is the
share of apparent homophily that was really triadic closure in disguise.

### Directed networks

Choose which shared-partner shape counts, following ergm's convention. For an
edge *t* → *h*, node *k* is a shared partner under:

| `type` | Condition on *k* | Name |
|------|------------------|------|
| `"OTP"` (default) | *t* → *k* → *h* | outgoing two-path (transitivity) |
| `"ITP"` | *h* → *k* → *t* | incoming two-path (cyclicality) |
| `"RTP"` | *t* ↔ *k* ↔ *h* | reciprocated two-path |
| `"OSP"` | *t* → *k* ← *h* | outgoing shared partner |
| `"ISP"` | *t* ← *k* → *h* | incoming shared partner |

```r
data(sampson)   # directed liking network among monks, with a "group" attribute
summary(samplike ~ gwespmatch(0.5, type = "OTP", match = "group"))
```

### Estimating the decay (curved ERGM)

```r
fit <- ergm(faux.magnolia.high ~ edges + nodematch("Race")
                               + gwespmatch(0.5, fixed = FALSE, match = "Race"))
```

`decay` is then estimated jointly with the coefficient rather than fixed by you.

## The statistic

Let *q* = 1 − e<sup>−α</sup>, and let *D*<sub>ℓ</sub> be the number of matched
edges with exactly ℓ shared partners (under the chosen `type` and `homophily`
rule). Then

```
g(y; α) = e^α · Σ_ℓ (1 − q^ℓ) · D_ℓ
```

An edge with one shared partner contributes 1 regardless of α; further partners
add geometrically less. With `fixed = FALSE` the term becomes a curved ERGM
carrying *D*<sub>1</sub>…*D*<sub>cutoff</sub> and estimating α via the standard
GWDECAY map θ<sub>ℓ</sub> = θ·e<sup>α</sup>(1 − q<sup>ℓ</sup>).

See `?gwespmatch` for the change statistic and the curved-form gradient.

## Validation

The `triad` statistic is checked against:

- **Grund & Densley (2015), Table 1** — reproduces their published
  `GWESP_Match` values (5.00 / 5.22 / 5.39) exactly.
- **ergm itself** — `gwespmatch(d, match = "a")` equals
  `F(~gwesp(d, fixed = TRUE), ~nodematch("a"))` on random networks, for all five
  directed types and undirected.
- **Change statistics** — verified toggle-by-toggle against recomputed summary
  statistics via `ergm.godfather`, in both `homophily` modes.

Run them with:

```r
testthat::test_dir(system.file("tests/testthat", package = "ergm.gwespmatch"))
```

### Do you even need this package?

Depends which version you want.

**`homophily = "triad"` — often not.** ergm's operator syntax already gets you
there, and this is exactly how Hong et al. (2024) did it:

```r
gwespmatch(0.5, match = "Race", homophily = "triad")   # these two
F(~gwesp(0.5, fixed = TRUE), ~nodematch("Race"))       # are the same number
```

`F()` restricts the network to within-group edges and then runs `gwesp` on that
subgraph. Since *every* edge of a homogeneous triad is a within-group edge, the
subgraph keeps exactly the triads you want. If you just need the Grund & Densley
statistic once, `F()` is the zero-dependency route.

**`homophily = "dyad"` — there is no `F()` equivalent.** `F()` filters the
*network*, which is all-or-nothing: it strips cross-group edges from the
two-path legs (*t*–*k* and *k*–*h*) at the same time as the focal edge. But
`"dyad"` needs the focal edge filtered while the legs are still evaluated on the
**full** graph — an out-group *k* must remain visible as a shared partner. No
combination of ergm's operators expresses that, because the restriction applies
to the focal edge only, not to the configuration the operator sees.

So: use `F()` if you want the published statistic and nothing else. Use this
package for the `"dyad"` variant, for curved-ERGM support (`fixed = FALSE`), or
for the purpose-written C change statistics (faster in long MCMC runs).

### A note on RTP

`gwespmatch` implements ergm's *documented* RTP rule (*k* is an RTP partner of
(*i*,*j*) iff *i* ↔ *k* ↔ *j*). **ergm's shared-partner cache miscomputes
edgewise RTP**: `dgwesp(type = "RTP")` only agrees with the documentation — and
with this package — when called with `term.options = list(cache.sp = FALSE)`.
Keep that in mind if you compare the two.

## References

- Grund, T. U. & Densley, J. A. (2015). Ethnic Homophily and Triad Closure:
  Mapping Internal Gang Structure Using Exponential Random Graph Models.
  *Journal of Contemporary Criminal Justice*, 31(3), 354–370.
- Hong, C.-S., Paik, A., Ballakrishnen, S., Silver, C. & Boutcher, S. (2024).
  Categorical closure: Transitivity and identities in longitudinal networks.
  *Social Networks*, 79, 76–92.
- Hunter, D. R. (2007). Curved exponential family models for social networks.
  *Social Networks*, 29(2), 216–230.

## Author

Benjamin Rosche (<benrosche@nyu.edu>) · GPL-3
