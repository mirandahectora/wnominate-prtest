#include "R.h"
/*
 *
 * Wrappers to call R random number stuff from Fortran
 *
 */
void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
float F77_SUB(random)() { return (float)unif_rand(); }


