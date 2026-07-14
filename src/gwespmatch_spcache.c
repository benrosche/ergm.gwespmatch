/*
 * Shared-partner cache for gwespmatch.
 *
 * This is a port of ergm's changestats_spcache.c with one addition: an optional
 * attribute filter.  When the filter is on ("triad" homophily) the cache is
 * maintained on the WITHIN-GROUP SUBGRAPH -- i.e. cross-group edges are treated
 * as if they did not exist.
 *
 * That is precisely what makes it correct for "triad".  A two-path a - k - b
 * survives a nodematch filter iff BOTH legs (a,k) and (k,b) are within-group,
 * which holds iff attr[a] == attr[k] == attr[b] -- the triad condition.  So a
 * shared-partner cache built on the filtered subgraph counts exactly the
 * matching shared partners.  This is the same thing ergm's
 * F(~gwesp, ~nodematch) does, and it is why that composition is fast.
 *
 * With the filter off ("dyad" homophily) the cache is the plain, unfiltered
 * shared-partner cache: every k counts, which is what "dyad" needs.
 *
 * Either way, c_gwespmatch() can then read a shared-partner count in O(1)
 * instead of recomputing it in O(degree), which takes a toggle from
 * O(degree^2) to O(degree).
 *
 * INPUT_PARAM layout for this auxiliary:
 *   [0]     type code: 0=UTP/undirected, 1=OTP, 2=ITP, 3=RTP, 4=OSP, 5=ISP
 *   [1]     match_k: 1 = filter to within-group subgraph ("triad"), 0 = no filter
 *   [2..]   vertex attribute values for nodes 1..N_NODES
 *
 * Cache kind by type (mirrors ergm):
 *   OTP/ITP : DIRECTED map;  cache(i,k) = # j with i->j->k.
 *             ITP shares the OTP cache, read with the arguments reversed.
 *   OSP     : UNDIRECTED map; cache{i,k} = # j with i->j and k->j.
 *   ISP     : UNDIRECTED map; cache{i,k} = # j with j->i and j->k.
 *   RTP     : UNDIRECTED map; cache{i,k} = # j with i<->j<->k.
 *   UTP     : UNDIRECTED map; cache{i,k} = # common neighbours.
 */

#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_dyad_hashmap.h"

#define AUX_TYPE    ((int)INPUT_PARAM[0])
#define AUX_MATCHK  ((int)INPUT_PARAM[1])
#define AUX_ATTRS   (INPUT_PARAM + 2)

/* Is edge (a,b) present in the (possibly filtered) network we cache over?
   The caller has already established that the edge exists in nwp. */
#define EOK(a, b) (!match_k || attrs[(a) - 1] == attrs[(b) - 1])

/* ---------------------------------------------------------------------------
 * Initialiser: build the cache from scratch.
 * ------------------------------------------------------------------------- */
I_CHANGESTAT_FN(i__gwespmatch_spcache) {
  StoreStrictDyadMapUInt *spcache = AUX_STORAGE = kh_init(StrictDyadMapUInt);
  const int     type    = AUX_TYPE;
  const int     match_k = AUX_MATCHK;
  const double *attrs   = AUX_ATTRS;

  EXEC_THROUGH_NET_EDGES(i, j, e1, {                 /* edge i->j (i<j if undirected) */
      if (!EOK(i, j)) continue;                      /* not in the filtered subgraph */

      switch (type) {

      case 0:                                        /* UTP: common neighbours */
        EXEC_THROUGH_EDGES(i, e2, k, {               /* k ~ i  =>  {j,k} share i */
            if (j < k && EOK(i, k)) IncUDyadMapUInt(j, k, 1, spcache);
          });
        EXEC_THROUGH_EDGES(j, e2, k, {               /* k ~ j  =>  {i,k} share j */
            if (i < k && EOK(j, k)) IncUDyadMapUInt(i, k, 1, spcache);
          });
        break;

      case 1:                                        /* OTP (and ITP): i->j->k */
      case 2:
        EXEC_THROUGH_FOUTEDGES(j, e2, k, {
            if (i != k && EOK(j, k)) IncDDyadMapUInt(i, k, 1, spcache);
          });
        break;

      case 3:                                        /* RTP: i<->j and i<->k */
        if (IS_OUTEDGE(j, i)) {                      /* i<->j mutual */
          EXEC_THROUGH_FOUTEDGES(i, e2, k, {
              if (j < k && IS_OUTEDGE(k, i) && EOK(i, k))
                IncUDyadMapUInt(j, k, 1, spcache);   /* {j,k} share mutual nbr i */
            });
        }
        break;

      case 4:                                        /* OSP: i->j and k->j */
        EXEC_THROUGH_FINEDGES(j, e2, k, {
            if (i < k && EOK(k, j)) IncUDyadMapUInt(i, k, 1, spcache);
          });
        break;

      case 5:                                        /* ISP: i->j and i->k */
        EXEC_THROUGH_FOUTEDGES(i, e2, k, {
            if (j < k && EOK(i, k)) IncUDyadMapUInt(j, k, 1, spcache);
          });
        break;

      default:
        break;
      }
    });
}

/* ---------------------------------------------------------------------------
 * Updater: incrementally maintain the cache when (tail,head) is toggled.
 * ------------------------------------------------------------------------- */
U_CHANGESTAT_FN(u__gwespmatch_spcache) {
  GET_AUX_STORAGE(StoreStrictDyadMapUInt, spcache);
  const int     type    = AUX_TYPE;
  const int     match_k = AUX_MATCHK;
  const double *attrs   = AUX_ATTRS;

  /* A cross-group edge is not part of the filtered subgraph, so toggling it
     changes nothing there. */
  if (!EOK(tail, head)) return;

  const int echange = edgestate ? -1 : 1;

  switch (type) {

  case 0:                                            /* UTP */
    EXEC_THROUGH_EDGES(tail, e, k, {                 /* head gains/loses nbr k of tail */
        if (head != k && EOK(tail, k)) IncUDyadMapUInt(head, k, echange, spcache);
      });
    EXEC_THROUGH_EDGES(head, e, k, {                 /* tail gains/loses nbr k of head */
        if (tail != k && EOK(head, k)) IncUDyadMapUInt(tail, k, echange, spcache);
      });
    break;

  case 1:                                            /* OTP (and ITP) */
  case 2:
    EXEC_THROUGH_FOUTEDGES(head, e, k, {             /* new t->h->k two-paths */
        if (tail != k && EOK(head, k)) IncDDyadMapUInt(tail, k, echange, spcache);
      });
    EXEC_THROUGH_FINEDGES(tail, e, k, {              /* new k->t->h two-paths */
        if (k != head && EOK(k, tail)) IncDDyadMapUInt(k, head, echange, spcache);
      });
    break;

  case 3:                                            /* RTP */
    /* Only matters if the reciprocating edge exists, i.e. if this toggle makes
       the {tail,head} dyad mutual (or un-mutual). */
    if (!IS_OUTEDGE(head, tail)) return;
    EXEC_THROUGH_FOUTEDGES(tail, e, k, {             /* k <-> tail  =>  tail joins {head,k} */
        if (head != k && IS_OUTEDGE(k, tail) && EOK(tail, k))
          IncUDyadMapUInt(head, k, echange, spcache);
      });
    EXEC_THROUGH_FOUTEDGES(head, e, k, {             /* k <-> head  =>  head joins {tail,k} */
        if (tail != k && IS_OUTEDGE(k, head) && EOK(head, k))
          IncUDyadMapUInt(tail, k, echange, spcache);
      });
    break;

  case 4:                                            /* OSP: t->h<-k */
    EXEC_THROUGH_FINEDGES(head, e, k, {
        if (tail != k && EOK(k, head)) IncUDyadMapUInt(tail, k, echange, spcache);
      });
    break;

  case 5:                                            /* ISP: h<-t->k */
    EXEC_THROUGH_FOUTEDGES(tail, e, k, {
        if (head != k && EOK(tail, k)) IncUDyadMapUInt(head, k, echange, spcache);
      });
    break;

  default:
    break;
  }
}

/* ---------------------------------------------------------------------------
 * Finaliser.
 * ------------------------------------------------------------------------- */
F_CHANGESTAT_FN(f__gwespmatch_spcache) {
  GET_AUX_STORAGE(StoreStrictDyadMapUInt, spcache);
  kh_destroy(StrictDyadMapUInt, spcache);
  AUX_STORAGE = NULL;
}
