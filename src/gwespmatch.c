#include "ergm_userterms.h"
#include <math.h>

/*
 * gwespmatch / gwespmatch_dist
 *
 * INPUT_PARAM layout (both functions share the same offset scheme):
 *   [0]  decay  (gwespmatch)  OR  cutoff cast to double (gwespmatch_dist)
 *   [1]  type code: 1=OTP, 2=ITP, 3=RTP, 4=OSP, 5=ISP
 *          Ignored for undirected networks (common neighbours used).
 *   [2..N+1]  vertex attribute values for nodes 1..N_NODES
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
 * Only edges where attr[tail] == attr[head] contribute.
 */

/* -------------------------------------------------------------------------
 * Helpers: count ESP partners for a single edge.
 * k == a and k == b are excluded automatically by the iteration or checks.
 * ------------------------------------------------------------------------- */

static int count_directed_esp(Vertex a, Vertex b, int type, Network *nwp)
{
    int count = 0;
    Edge e; Vertex k;
    switch (type) {
    case 1: /* OTP: a->k AND k->b */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && IS_OUTEDGE(k, b)) count++;
        }
        break;
    case 2: /* ITP: k->a AND b->k */
        STEP_THROUGH_INEDGES(a, e, k) {
            if (k != b && IS_OUTEDGE(b, k)) count++;
        }
        break;
    case 3: /* RTP: a<->k AND b<->k */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && IS_OUTEDGE(k, a) && IS_OUTEDGE(b, k) && IS_OUTEDGE(k, b))
                count++;
        }
        break;
    case 4: /* OSP: a->k AND b->k */
        STEP_THROUGH_OUTEDGES(a, e, k) {
            if (k != b && IS_OUTEDGE(b, k)) count++;
        }
        break;
    case 5: /* ISP: k->a AND k->b */
        STEP_THROUGH_INEDGES(a, e, k) {
            if (k != b && IS_OUTEDGE(k, b)) count++;
        }
        break;
    default:
        break;
    }
    return count;
}

static int count_undirected_esp(Vertex a, Vertex b, Network *nwp)
{
    /* Common neighbours of a and b. */
    int count = 0;
    Edge e; Vertex k;
    STEP_THROUGH_OUTEDGES(a, e, k) {   /* k > a */
        if (k != b && IS_UNDIRECTED_EDGE(k, b)) count++;
    }
    STEP_THROUGH_INEDGES(a, e, k) {    /* k < a */
        if (k != b && IS_UNDIRECTED_EDGE(k, b)) count++;
    }
    return count;
}

/* =========================================================================
 * s_gwespmatch  –  summary statistic (fixed = TRUE, returns 1 value)
 * ========================================================================= */
S_CHANGESTAT_FN(s_gwespmatch) {
    const double  decay  = INPUT_PARAM[0];
    const int     type   = (int)INPUT_PARAM[1];
    const double *attrs  = INPUT_PARAM + 2;

    CHANGE_STAT[0] = 0.0;

    const int    decay_is_zero = (decay == 0.0);
    const double q             = decay_is_zero ? 0.0 : 1.0 - exp(-decay);
    const double exp_decay     = decay_is_zero ? 0.0 : exp(decay);

    for (Vertex tail = 1; tail <= N_NODES; tail++) {
        Edge e; Vertex head;
        STEP_THROUGH_OUTEDGES(tail, e, head) {
            /* For undirected, STEP_THROUGH_OUTEDGES gives head > tail,
             * so each edge is visited exactly once. */
            if (attrs[tail - 1] != attrs[head - 1]) continue;

            int esp;
            if (DIRECTED)
                esp = count_directed_esp(tail, head, type, nwp);
            else
                esp = count_undirected_esp(tail, head, nwp);

            if (esp == 0) continue;

            CHANGE_STAT[0] += decay_is_zero
                ? (double)esp
                : exp_decay * (1.0 - pow(q, esp));
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
 * ========================================================================= */
C_CHANGESTAT_FN(c_gwespmatch) {
    const double  decay  = INPUT_PARAM[0];
    const int     type   = (int)INPUT_PARAM[1];
    const double *attrs  = INPUT_PARAM + 2;

    CHANGE_STAT[0] = 0.0;

    const double attr_t = attrs[tail - 1];
    const double attr_h = attrs[head - 1];
    if (attr_t != attr_h) return;

    const int    echange       = edgestate ? -1 : 1;
    const int    decay_is_zero = (decay == 0.0);
    const double q             = decay_is_zero ? 0.0 : 1.0 - exp(-decay);
    double cumchange = 0.0;

    /* Macro: accumulate contribution of an existing matched edge (A,B)
     * that gains or loses one directed ESP partner.  For decay == 0 the
     * marginal is always 1; otherwise it is q^(L2 - edgestate). */
#define ADD_CONTRIB_D(A, B)                                              \
    do {                                                                 \
        if (attrs[(A) - 1] == attrs[(B) - 1]) {                         \
            if (decay_is_zero) {                                         \
                cumchange += 1.0;                                        \
            } else {                                                     \
                int _L = count_directed_esp((A), (B), type, nwp);       \
                cumchange += pow(q, _L - edgestate);                    \
            }                                                            \
        }                                                                \
    } while (0)

#define ADD_CONTRIB_U(A, B)                                              \
    do {                                                                 \
        if (attrs[(A) - 1] == attrs[(B) - 1]) {                         \
            if (decay_is_zero) {                                         \
                cumchange += 1.0;                                        \
            } else {                                                     \
                int _L = count_undirected_esp((A), (B), nwp);           \
                cumchange += pow(q, _L - edgestate);                    \
            }                                                            \
        }                                                                \
    } while (0)

    if (!DIRECTED) {
        /* ----------------------------------------------------------------
         * Undirected: common-neighbour change statistic.
         * ---------------------------------------------------------------- */

        /* Part 1: edge {tail,head} itself */
        int L2_th = count_undirected_esp(tail, head, nwp);
        if (decay_is_zero)
            cumchange += (double)L2_th;
        else
            cumchange += exp(decay) * (1.0 - pow(q, L2_th));

        Edge e; Vertex v;

        /* Part 2: edges {tail,v} for v in N(head), v != tail */
        STEP_THROUGH_OUTEDGES(head, e, v) {       /* v > head > tail */
            if (IS_UNDIRECTED_EDGE(tail, v)) ADD_CONTRIB_U(tail, v);
        }
        STEP_THROUGH_INEDGES(head, e, v) {        /* v < head */
            if (v != tail && IS_UNDIRECTED_EDGE(tail, v)) ADD_CONTRIB_U(tail, v);
        }

        /* Part 3: edges {head,v} for v in N(tail), v != head */
        STEP_THROUGH_OUTEDGES(tail, e, v) {       /* v > tail */
            if (v != head && IS_UNDIRECTED_EDGE(head, v)) ADD_CONTRIB_U(head, v);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {        /* v < tail < head */
            if (IS_UNDIRECTED_EDGE(head, v)) ADD_CONTRIB_U(head, v);
        }

        CHANGE_STAT[0] = echange * cumchange;
        return;
    }

    /* ------------------------------------------------------------------
     * Directed: type-specific change statistic.
     * ------------------------------------------------------------------ */

    /* Part 1: the (tail->head) edge itself */
    {
        int L2_th = count_directed_esp(tail, head, type, nwp);
        if (decay_is_zero)
            cumchange += (double)L2_th;
        else
            cumchange += exp(decay) * (1.0 - pow(q, L2_th));
    }

    Edge e; Vertex v;
    switch (type) {

    case 1: /* OTP: k: tail->k AND k->head
             * Part 2: (u,head) where u->tail  [k=tail is new OTP partner]
             * Part 3: (tail,v) where head->v  [k=head is new OTP partner] */
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(v, head)) ADD_CONTRIB_D(v, head);
        }
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(tail, v)) ADD_CONTRIB_D(tail, v);
        }
        break;

    case 2: /* ITP: k: k->tail AND head->k
             * Adding tail->head makes tail a new in-nbr of head.
             * Part 2: (head,b) where b->tail  [k=tail is new ITP partner]
             * Part 3: (a,tail) where head->a  [k=head is new ITP partner] */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(v, tail)) ADD_CONTRIB_D(head, v);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(head, v)) ADD_CONTRIB_D(v, tail);
        }
        break;

    case 3: /* RTP: k: tail<->k AND head<->k
             * k=tail can become new RTP partner only if head->tail exists.
             * Part 2: (v,head) where v<->tail  [k=tail new RTP partner]
             * Part 3: (head,v) where v<->tail  [k=tail new RTP partner] */
        if (IS_OUTEDGE(head, tail)) {
            STEP_THROUGH_OUTEDGES(tail, e, v) {
                if (v != head && IS_OUTEDGE(v, tail) && IS_OUTEDGE(v, head))
                    ADD_CONTRIB_D(v, head);
            }
            STEP_THROUGH_OUTEDGES(head, e, v) {
                if (v != tail && IS_OUTEDGE(v, tail) && IS_OUTEDGE(tail, v))
                    ADD_CONTRIB_D(head, v);
            }
        }
        break;

    case 4: /* OSP: k: tail->k AND head->k
             * Adding tail->head makes head a new out-nbr of tail.
             * Part 2: (tail,v) where v->head  [k=head new OSP partner]
             * Part 3: (v,tail) where v->head  [k=head new OSP partner] */
        STEP_THROUGH_INEDGES(head, e, v) {     /* v->head */
            if (v != tail) {
                if (IS_OUTEDGE(tail, v)) ADD_CONTRIB_D(tail, v);
                if (IS_OUTEDGE(v, tail)) ADD_CONTRIB_D(v, tail);
            }
        }
        break;

    case 5: /* ISP: k: k->tail AND k->head
             * Adding tail->head makes tail a new in-nbr of head.
             * Part 2: (v,head) where tail->v  [k=tail new ISP partner]
             * Part 3: (head,v) where tail->v  [k=tail new ISP partner] */
        STEP_THROUGH_OUTEDGES(tail, e, v) {    /* tail->v */
            if (v != head) {
                if (IS_OUTEDGE(v, head)) ADD_CONTRIB_D(v, head);
                if (IS_OUTEDGE(head, v)) ADD_CONTRIB_D(head, v);
            }
        }
        break;

    default:
        break;
    }

    CHANGE_STAT[0] = echange * cumchange;

#undef ADD_CONTRIB_D
#undef ADD_CONTRIB_U
}

/* =========================================================================
 * s_gwespmatch_dist  –  summary statistic (fixed = FALSE)
 *
 * Returns cutoff = N_CHANGE_STATS values: CHANGE_STAT[l-1] = number of
 * matched edges with exactly l shared partners (l = 1..cutoff).
 * ========================================================================= */
S_CHANGESTAT_FN(s_gwespmatch_dist) {
    const int     cutoff = (int)INPUT_PARAM[0];   /* == N_CHANGE_STATS */
    const int     type   = (int)INPUT_PARAM[1];
    const double *attrs  = INPUT_PARAM + 2;

    for (int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] = 0.0;

    for (Vertex tail = 1; tail <= N_NODES; tail++) {
        Edge e; Vertex head;
        STEP_THROUGH_OUTEDGES(tail, e, head) {
            if (attrs[tail - 1] != attrs[head - 1]) continue;

            int esp;
            if (DIRECTED)
                esp = count_directed_esp(tail, head, type, nwp);
            else
                esp = count_undirected_esp(tail, head, nwp);

            if (esp >= 1 && esp <= cutoff)
                CHANGE_STAT[esp - 1] += 1.0;
        }
    }
    (void)cutoff;   /* suppress unused-variable warning */
}

/* =========================================================================
 * c_gwespmatch_dist  –  change statistic (fixed = FALSE)
 *
 * Updates the matched ESP distribution bins when edge (tail,head) is toggled.
 *
 * For affected matched edge (a,b) with current count L (= count before toggle):
 *   CHANGE_STAT[L-1]       -= 1  (it leaves bin L)
 *   CHANGE_STAT[L+echange-1]+= 1  (it enters bin L+echange)
 * ========================================================================= */
C_CHANGESTAT_FN(c_gwespmatch_dist) {
    const int     cutoff = (int)INPUT_PARAM[0];
    const int     type   = (int)INPUT_PARAM[1];
    const double *attrs  = INPUT_PARAM + 2;

    for (int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] = 0.0;

    const double attr_t = attrs[tail - 1];
    const double attr_h = attrs[head - 1];
    if (attr_t != attr_h) return;

    const int echange = edgestate ? -1 : 1;

    /* Move edge (A,B) from its current bin to the adjacent bin. */
#define MOVE_BIN_D(A, B)                                                      \
    do {                                                                       \
        if (attrs[(A) - 1] == attrs[(B) - 1]) {                               \
            int _L = count_directed_esp((A), (B), type, nwp);                 \
            if (_L >= 1 && _L <= N_CHANGE_STATS) CHANGE_STAT[_L - 1] -= 1.0; \
            int _Ln = _L + echange;                                            \
            if (_Ln >= 1 && _Ln <= N_CHANGE_STATS) CHANGE_STAT[_Ln-1] += 1.0;\
        }                                                                      \
    } while (0)

#define MOVE_BIN_U(A, B)                                                      \
    do {                                                                       \
        if (attrs[(A) - 1] == attrs[(B) - 1]) {                               \
            int _L = count_undirected_esp((A), (B), nwp);                     \
            if (_L >= 1 && _L <= N_CHANGE_STATS) CHANGE_STAT[_L - 1] -= 1.0; \
            int _Ln = _L + echange;                                            \
            if (_Ln >= 1 && _Ln <= N_CHANGE_STATS) CHANGE_STAT[_Ln-1] += 1.0;\
        }                                                                      \
    } while (0)

    if (!DIRECTED) {
        /* Part 1: toggled edge {tail,head} itself */
        int L_th = count_undirected_esp(tail, head, nwp);
        if (L_th >= 1 && L_th <= N_CHANGE_STATS) CHANGE_STAT[L_th - 1] += echange;

        Edge e; Vertex v;

        /* Part 2: {tail,v} for v in N(head) */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (IS_UNDIRECTED_EDGE(tail, v)) MOVE_BIN_U(tail, v);
        }
        STEP_THROUGH_INEDGES(head, e, v) {
            if (v != tail && IS_UNDIRECTED_EDGE(tail, v)) MOVE_BIN_U(tail, v);
        }

        /* Part 3: {head,v} for v in N(tail) */
        STEP_THROUGH_OUTEDGES(tail, e, v) {
            if (v != head && IS_UNDIRECTED_EDGE(head, v)) MOVE_BIN_U(head, v);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (IS_UNDIRECTED_EDGE(head, v)) MOVE_BIN_U(head, v);
        }

        (void)cutoff;
        return;
    }

    /* Directed */

    /* Part 1: (tail->head) edge itself */
    {
        int L_th = count_directed_esp(tail, head, type, nwp);
        if (L_th >= 1 && L_th <= N_CHANGE_STATS) CHANGE_STAT[L_th - 1] += echange;
    }

    Edge e; Vertex v;
    switch (type) {

    case 1: /* OTP */
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(v, head)) MOVE_BIN_D(v, head);
        }
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(tail, v)) MOVE_BIN_D(tail, v);
        }
        break;

    case 2: /* ITP */
        STEP_THROUGH_OUTEDGES(head, e, v) {
            if (v != tail && IS_OUTEDGE(v, tail)) MOVE_BIN_D(head, v);
        }
        STEP_THROUGH_INEDGES(tail, e, v) {
            if (v != head && IS_OUTEDGE(head, v)) MOVE_BIN_D(v, tail);
        }
        break;

    case 3: /* RTP */
        if (IS_OUTEDGE(head, tail)) {
            STEP_THROUGH_OUTEDGES(tail, e, v) {
                if (v != head && IS_OUTEDGE(v, tail) && IS_OUTEDGE(v, head))
                    MOVE_BIN_D(v, head);
            }
            STEP_THROUGH_OUTEDGES(head, e, v) {
                if (v != tail && IS_OUTEDGE(v, tail) && IS_OUTEDGE(tail, v))
                    MOVE_BIN_D(head, v);
            }
        }
        break;

    case 4: /* OSP */
        STEP_THROUGH_INEDGES(head, e, v) {
            if (v != tail) {
                if (IS_OUTEDGE(tail, v)) MOVE_BIN_D(tail, v);
                if (IS_OUTEDGE(v, tail)) MOVE_BIN_D(v, tail);
            }
        }
        break;

    case 5: /* ISP */
        STEP_THROUGH_OUTEDGES(tail, e, v) {
            if (v != head) {
                if (IS_OUTEDGE(v, head))  MOVE_BIN_D(v, head);
                if (IS_OUTEDGE(head, v))  MOVE_BIN_D(head, v);
            }
        }
        break;

    default:
        break;
    }

    (void)cutoff;

#undef MOVE_BIN_D
#undef MOVE_BIN_U
}
