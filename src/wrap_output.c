#include "R.h"
/*
 *
 *   Wrappers to call R win32 process/output control stuff from Fortran
 *
 */

void R_FlushConsole(void);
void R_ProcessEvents(void);

void F77_SUB(flushcon)(void) { R_FlushConsole(); }
void F77_SUB(procevent)(void) {
#ifdef __win32
    R_ProcessEvents();
#endif /* __win32 */
}

void F77_SUB(echoevent)(int *action){
    const int a = *action;
    switch(a) {
        case 0:
            Rprintf(".");
            break;
        case 1:
            Rprintf("\t\tGetting legislator coordinates...\n");
            break;
        case 2:
            Rprintf("\t\tGetting bill parameters...\n");
            break;
        case 3:
            Rprintf("\t\tStarting bootstrap iterations...\n");
            break;
        case 4:
            Rprintf("\t\tStarting estimation of Beta...\n");
            break;
        case 5:
            Rprintf("\t\tEstimating weights...\n");
            break;
        case 6:
            Rprintf("\t\tComputing standard errors...\n");
            break;
        case 7:
            Rprintf("\n\n");
            break;
        case 8:
            Rf_error("Data set too small to recover estimates: GRID2() failed in wnom9707().\n");
            break;           
    }       
}
