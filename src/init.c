#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Change-statistic functions are discovered by name at run time. */
void R_init_ergm_gwespmatch(DllInfo *dll) {
  R_useDynamicSymbols(dll, TRUE);
}
