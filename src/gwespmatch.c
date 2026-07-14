#include "ergm_userterms.h"
#include <math.h>

/*
 * gwespmatch / gwespmatch_dist
 *
 * INPUT_PARAM layout (both functions share the same offset scheme):
 *   [0]  decay  (gwespmatch)  OR  cutoff cast to double (gwespmatch_dist)
 *   [1]  type code: 1=OTP, 2=ITP, 3=RTP, 4=OSP, 5=ISP
 *          Ignored for undirected networks (common neighbours used).
 *   [2]  match_k: 1 = the shared partner k must ALSO carry the same attribute
 *                     value ("triad" homophily; Grund & Densley 2015,
 *                     Hong et al. 2024).
 *                 0 = only the focal edge (tail,head) must match; k is
 *                     unconstrained ("dyad" homophily).
 *   [3..N+2]  vertex attribute values for nodes 1..N_NODES
 *
 * Directed shared-partner types for edge (tail -> head):
 *   OTP (1): k s.t.  tail->k  AND  k->head
 *   ITP (2): k s.t.  k->tail  AND  head->k
 *   RTP (3): k s.t.  tail<->k AND  head<->k  (mutually adjacent to both)
 *   OSP (4): k s.t.  tail->k  AND  head->k
 *   ISP (5): k s.t.  k->tail  AND  k->head
 *
 * Undirected: common neighbours regardless of type code.
 *
 * Only edges where attr[tail] == attr[head] contribute.  Because the ESP
 * helpers are only ever called on such matched edges, "k matches a" is
 * equivalent to "k matches b", i.e. to "all three nodes match".
 *
 * Note on decay == 0: there is deliberately NO special case.  With
 * q = 1 - exp(-0) = 0, the general formulae reduce correctly, since C's
 * pow(0,0) == 1:
 *     exp(0)*(1 - pow(0, L))     = 1 iff L >= 1, else 0
 *     pow(0, L - edgestate)      = 1 iff L == edgestate
 * which is exactly the limit of the GW statistic as decay -> 0, i.e. every
 * matched edge with at least one shared partner counts 1.
 */

/* -------------------------------------------------------------------------
 * Helpers: count ESP partners for a single edge (a,b), which the caller has
 * already verified is matched (attrs[a-1] == attrs[b-1]).
 * When match_k is set, only partners k with attrs[k-1] == attrs[a-1] count.
 * ------------------------------------------------------------------------- */

#define K_MATCHES(K) (!match_k || attrs[(K) - 1] == ref)

static int count_directed_esp(Vertex a, Vertex b, int type,
                              const double *attrs, int match_k, Network *nwp)
{
    int count = 0;
    const double ref = attrs[a - 1];
    Edge e; Vertex k;
    switch (type) {
    case 1: /* OTP: a->k AND k->b */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && K_MATCHES(k) && IS_OUTEDGE(k, b)) count++;
        }
        break;
    case 2: /* ITP: k->a AND b->k */
        STEP_THROUGH_INEDGES(a, e, k) {
            if (k != b && K_MATCHES(k) && IS_OUTEDGE(b, k)) count++;
        }
        break;
    case 3: /* RTP: a<->k AND b<->k */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && K_MATCHES(k) &&
                IS_OUTEDGE(k, a) && IS_OUTEDGE(b, k) && IS_OUTEDGE(k, b))
                count++;
        }
        break;
    case 4: /* OSP: a->k AND b->k */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && K_MATCHES(k) && IS_OUTEDGE(b, k)) count++;
        }
        break;
    case 5: /* ISP: k->a AND k->b */
        STEP_THROUGH_INEDGES(a, e, k) {
            if (k != b && K_MATCHES(k) && IS_OUTEDGE(k, b)) count++;
        }
        break;
    default:
        break;
    }
    return count;
}

static int count_undirected_esp(Vertex a, Vertex b,
                                const double *attrs, int match_k, Network *nwp)
{
    /* Common neighbours of a and b. */
    int count = 0;
    const double ref = attrs[a - 1];
    Edge e; Vertex k;
    STEP_THROUGH_OUTEDGES(a, e, k) {   /* k > a */
        if (k != b && K_MATCHES(k) && IS_UNDIRECTED_EDGE(k, b)) count++;
    }
    STEP_THROUGH_INEDGES(a, e, k) {    /* k < a */
        if (k != b && K_MATCHES(k) && IS_UNDIRECTED_EDGE(k, b)) count++;
    }
    return count;
}

/* Count ESP for edge (a,b) dispatching on directedness. */
#define COUNT_ESP(A, B) (DIRECTED                                            \
    ? count_directed_esp((A), (B), type, attrs, match_k, nwp)                \
    : count_undirected_esp((A), (B), attrs, match_k, nwp))

/* =========================================================================
 * s_gwespmatch  –  summary statistic (fixed = TRUE, returns 1 value)
 * ========================================================================= */
S_CHANGESTAT_FN(s_gwespmatch) {
    const double  decay   = INPUT_PARAM[0];
    const int     type    = (int)INPUT_PARAM[1];
    const int     match_k = (int)INPUT_PARAM[2];
    const double *attrs   = INPUT_PARAM + 3;

    CHANGE_STAT[0] = 0.0;

    const double q         = 1.0 - exp(-decay);
    const double exp_decay = exp(decay);

    for (Vertex tail = 1; tail <= N_NODES; tail++) {
        Edge e; Vertex head;
        STEP_THROUGH_OUTEDGES(tail, e, head) {
            /* For undirected, STEP_THROUGH_OUTEDGES gives head > tail,
             * so each edge is visited exactly once. */
            if (attrs[tail - 1] != attrs[head - 1]) continue;

            int esp = COUNT_ESP(tail, head);
            if (esp == 0) continue;

            CHANGE_STAT[0] += exp_decay * (1.0 - pow(q, esp));
        }
    }
}

/* =========================================================================
 * c_gwespmatch  –  change statistic (fixed = TRUE, returns 1 value)
 *
 * When edge (tail,head) is toggled, three groups of edges are affected:
 *   1. The edge (tail,head) itself.
 *   2. Existing edges whose type-X ESP count gains/loses tail as a partner.
 *   3. Existing edges whose type-X ESP count gains/loses head as a partner.
 *
 * For each affected matched edge (a,b) with current count L:
 *   marginal weight = q^(L - edgestate)
 *   (edgestate correction: when removing, L already includes the partner)
 *
 * Note: the set of affected edges does NOT depend on match_k.  If the focal
 * edge (tail,head) is matched and an affected edge (a,b) is matched, then the
 * newly gained partner (tail or head) necessarily shares their attribute
 * value, so it counts under "triad" homophily automatically.  Only the
 * current count L is affected by match_k, and COUNT_ESP handles that.
 * ========================================================================= */
C_CHANGESTAT_FN(c_gwespmatch) {
    const double  decay   = INPUT_PARAM[0];
    const int     type    = (int)INPUT_PARAM[1];
    const int     match_k = (int)INPUT_PARAM[2];
    const double *attrs   = INPUT_PARAM + 3;

    CHANGE_STAT[0] = 0.0;

    const int matched_focal = (attrs[tail - 1] == attrs[head - 1]);

    /* Under "triad" homophily every contribution below is provably zero when
     * the focal edge is unmatched: the only nodes that can become new shared
     * partners are tail and head, and neither can match an affected matched
     * edge unless attrs[tail] == attrs[head].  Under "dyad" homophily this is
     * NOT true -- a heterophilous edge can still close a matched edge -- so we
     * must carry on. */
    if (match_k && !matched_focal) return;

    const int    echange = edgestate ? -1 : 1;
    const double q       = 1.0 - exp(-decay);
    double cumchange     = 0.0;

    /* Contribution of an existing matched edge (A,B) that gains/loses node P
     * as a shared partner.  Under "triad", P must itself match. */
#define ADD_CONTRIB(A, B, P)                                                 \
    do {                                                                     \
        if (attrs[(A) - 1] == attrs[(B) - 1] &&                              \
            (!match_k || attrs[(P) - 1] == attrs[(A) - 1])) {                \
            int _L = COUNT_ESP((A), (B));                                    \
            cumchange += pow(q, _L - edgestate);                             \
        }                                                                    \
    } while (0)

    /* Part 1: the toggled edge's own term (only if it is itself matched). */
    if (matched_focal)
        cumchange += exp(decay) * (1.0 - pow(q, COUNT_ESP(tail, head)));

    Edge e; Vertex v;

    if (!DIRECTED) {
        /* ----------------------------------------------------------------
         * Undirected: common-neighbour change statistic.
         * Edges {tail,v} gain head as a partner; {head,v} gain tail.
         * ---------------------------------------------------------------- */

        /* Part 2: edges {tail,v} for v in N(head), v != tail  -> partner head */
        STEP_THROUGH_OUTEDGES(head, e, v) {       /* v > head > tail */
            if (IS_UNDIRECTED_EDGE(tail, v)) ADD_CONTRIB(tail, v, head);
        }
        STEP_THROUGH_INEDGES(head, e, v) {        /* v < head */
            if (v != tail && IS_UNDIRECTED_EDGE(tail, v)) ADD_CONTRIB(tail, v, head);
        }

        /* Part 3: edges {head,v} for v in N(tail), v != head  -> partner tail */
        STEP_THROUGH_OUTEDGES(tail, e, v) {       /* v > tail */
            if (v != head && IS_UNDIRECTED_EDGE(head, v)) ADD_CONTRIB(head, v, tail);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {        /* v < tail < head */
            if (IS_UNDIRECTED_EDGE(head, v)) ADD_CONTRIB(head, v, tail);
        }

        CHANGE_STAT[0] = echange * cumchange;
        return;
    }

    /* ------------------------------------------------------------------
     * Directed: type-specific change statistic.
     * ------------------------------------------------------------------ */
    switch (type) {

    case 1: /* OTP: k s.t. a->k->b.  New edge tail->head can serve as either leg.
             * as k->b leg (k=tail, b=head): edges (v,head) with v->tail  -> partner tail
             * as a->k leg (a=tail, k=head): edges (tail,v) with head->v  -> partner head */
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(v, head)) ADD_CONTRIB(v, head, tail);
        }
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(tail, v)) ADD_CONTRIB(tail, v, head);
        }
        break;

    case 2: /* ITP: k s.t. k->a and b->k.
             * as k->a leg (k=tail, a=head): edges (head,v) with v->tail  -> partner tail
             * as b->k leg (b=tail, k=head): edges (v,tail) with head->v  -> partner head */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(v, tail)) ADD_CONTRIB(head, v, tail);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(head, v)) ADD_CONTRIB(v, tail, head);
        }
        break;

    case 3: /* RTP: k s.t. a<->k and k<->b.
             * Toggling tail->head only matters if head->tail already exists, in
             * which case the dyad {tail,head} gains/loses mutuality.  BOTH tail
             * and head can then become new partners:
             *   k=tail : edges (v,head) and (head,v) for v <-> tail
             *   k=head : edges (v,tail) and (tail,v) for v <-> head   <-- was missing
             */
        if (IS_OUTEDGE(head, tail)) {
            STEP_THROUGH_OUTEDGES(tail, e, v) {            /* tail->v */
                if (v != head && IS_OUTEDGE(v, tail)) {    /* v <-> tail */
                    if (IS_OUTEDGE(v, head)) ADD_CONTRIB(v, head, tail);
                    if (IS_OUTEDGE(head, v)) ADD_CONTRIB(head, v, tail);
                }
            }
            STEP_THROUGH_OUTEDGES(head, e, v) {            /* head->v */
                if (v != tail && IS_OUTEDGE(v, head)) {    /* v <-> head */
                    if (IS_OUTEDGE(v, tail)) ADD_CONTRIB(v, tail, head);
                    if (IS_OUTEDGE(tail, v)) ADD_CONTRIB(tail, v, head);
                }
            }
        }
        break;

    case 4: /* OSP: k s.t. a->k and b->k.
             * Adding tail->head makes head a new out-nbr of tail, so head can
             * become a new partner of edges (tail,v) / (v,tail) with v->head. */
        STEP_THROUGH_INEDGES(head, e, v) {     /* v->head */
            if (v != tail) {
                if (IS_OUTEDGE(tail, v)) ADD_CONTRIB(tail, v, head);
                if (IS_OUTEDGE(v, tail)) ADD_CONTRIB(v, tail, head);
            }
        }
        break;

    case 5: /* ISP: k s.t. k->a and k->b.
             * Adding tail->head makes tail point at head, so tail can become a
             * new partner of edges (v,head) / (head,v) with tail->v. */
        STEP_THROUGH_OUTEDGES(tail, e, v) {    /* tail->v */
            if (v != head) {
                if (IS_OUTEDGE(v, head)) ADD_CONTRIB(v, head, tail);
                if (IS_OUTEDGE(head, v)) ADD_CONTRIB(head, v, tail);
            }
        }
        break;

    default:
        break;
    }

    CHANGE_STAT[0] = echange * cumchange;

#undef ADD_CONTRIB
}

/* =========================================================================
 * s_gwespmatch_dist  –  summary statistic (fixed = FALSE)
 *
 * Returns cutoff = N_CHANGE_STATS values: CHANGE_STAT[l-1] = number of
 * matched edges with exactly l shared partners (l = 1..cutoff).
 * ========================================================================= */
S_CHANGESTAT_FN(s_gwespmatch_dist) {
    const int     type    = (int)INPUT_PARAM[1];
    const int     match_k = (int)INPUT_PARAM[2];
    const double *attrs   = INPUT_PARAM + 3;

    for (int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] = 0.0;

    for (Vertex tail = 1; tail <= N_NODES; tail++) {
        Edge e; Vertex head;
        STEP_THROUGH_OUTEDGES(tail, e, head) {
            if (attrs[tail - 1] != attrs[head - 1]) continue;

            int esp = COUNT_ESP(tail, head);
            if (esp >= 1 && esp <= N_CHANGE_STATS)
                CHANGE_STAT[esp - 1] += 1.0;
        }
    }
}

/* =========================================================================
 * c_gwespmatch_dist  –  change statistic (fixed = FALSE)
 *
 * Updates the matched ESP distribution bins when edge (tail,head) is toggled.
 *
 * For affected matched edge (a,b) with current count L (before the toggle):
 *   CHANGE_STAT[L-1]         -= 1   (it leaves bin L)
 *   CHANGE_STAT[L+echange-1] += 1   (it enters bin L+echange)
 * ========================================================================= */
C_CHANGESTAT_FN(c_gwespmatch_dist) {
    const int     type    = (int)INPUT_PARAM[1];
    const int     match_k = (int)INPUT_PARAM[2];
    const double *attrs   = INPUT_PARAM + 3;

    for (int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] = 0.0;

    const int matched_focal = (attrs[tail - 1] == attrs[head - 1]);

    /* See c_gwespmatch: this early exit is only sound under "triad". */
    if (match_k && !matched_focal) return;

    const int echange = edgestate ? -1 : 1;

    /* Matched edge (A,B) gains/loses node P as a shared partner: move it from
     * bin L to bin L+echange.  Under "triad", P must itself match. */
#define MOVE_BIN(A, B, P)                                                    \
    do {                                                                     \
        if (attrs[(A) - 1] == attrs[(B) - 1] &&                              \
            (!match_k || attrs[(P) - 1] == attrs[(A) - 1])) {                \
            int _L = COUNT_ESP((A), (B));                                    \
            if (_L >= 1 && _L <= N_CHANGE_STATS) CHANGE_STAT[_L - 1] -= 1.0; \
            int _Ln = _L + echange;                                          \
            if (_Ln >= 1 && _Ln <= N_CHANGE_STATS) CHANGE_STAT[_Ln - 1] += 1.0; \
        }                                                                    \
    } while (0)

    /* Part 1: the toggled edge's own bin (only if it is itself matched). */
    if (matched_focal) {
        int L_th = COUNT_ESP(tail, head);
        if (L_th >= 1 && L_th <= N_CHANGE_STATS) CHANGE_STAT[L_th - 1] += echange;
    }

    Edge e; Vertex v;

    if (!DIRECTED) {
        /* {tail,v} gain head as partner; {head,v} gain tail. */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (IS_UNDIRECTED_EDGE(tail, v)) MOVE_BIN(tail, v, head);
        }
        STEP_THROUGH_INEDGES(head, e, v) {
            if (v != tail && IS_UNDIRECTED_EDGE(tail, v)) MOVE_BIN(tail, v, head);
        }
        STEP_THROUGH_OUTEDGES(tail, e, v) {
            if (v != head && IS_UNDIRECTED_EDGE(head, v)) MOVE_BIN(head, v, tail);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (IS_UNDIRECTED_EDGE(head, v)) MOVE_BIN(head, v, tail);
        }
        return;
    }

    /* Directed */
    switch (type) {

    case 1: /* OTP */
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(v, head)) MOVE_BIN(v, head, tail);
        }
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(tail, v)) MOVE_BIN(tail, v, head);
        }
        break;

    case 2: /* ITP */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(v, tail)) MOVE_BIN(head, v, tail);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(head, v)) MOVE_BIN(v, tail, head);
        }
        break;

    case 3: /* RTP -- both tail and head can become new partners */
        if (IS_OUTEDGE(head, tail)) {
            STEP_THROUGH_OUTEDGES(tail, e, v) {            /* v <-> tail */
                if (v != head && IS_OUTEDGE(v, tail)) {
                    if (IS_OUTEDGE(v, head)) MOVE_BIN(v, head, tail);
                    if (IS_OUTEDGE(head, v)) MOVE_BIN(head, v, tail);
                }
            }
            STEP_THROUGH_OUTEDGES(head, e, v) {            /* v <-> head */
                if (v != tail && IS_OUTEDGE(v, head)) {
                    if (IS_OUTEDGE(v, tail)) MOVE_BIN(v, tail, head);
                    if (IS_OUTEDGE(tail, v)) MOVE_BIN(tail, v, head);
                }
            }
        }
        break;

    case 4: /* OSP */
        STEP_THROUGH_INEDGES(head, e, v) {
            if (v != tail) {
                if (IS_OUTEDGE(tail, v)) MOVE_BIN(tail, v, head);
                if (IS_OUTEDGE(v, tail)) MOVE_BIN(v, tail, head);
            }
        }
        break;

    case 5: /* ISP */
        STEP_THROUGH_OUTEDGES(tail, e, v) {
            if (v != head) {
                if (IS_OUTEDGE(v, head)) MOVE_BIN(v, head, tail);
                if (IS_OUTEDGE(head, v)) MOVE_BIN(head, v, tail);
            }
        }
        break;

    default:
        break;
    }

#undef MOVE_BIN
}
