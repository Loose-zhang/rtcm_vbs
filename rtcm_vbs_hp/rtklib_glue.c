/*------------------------------------------------------------------------------
 * rtklib_glue.c : provide storage for RTKLIB globals that are only declared
 *                 (extern) in rtklib.h but never defined in any of the .c
 *                 files we compile. Needed so satposs() / ephpos() from
 *                 ephemeris.c can link.
 *----------------------------------------------------------------------------*/
#include "rtklib.h"

/* All-zero is fine for VBS: RT_flag=0 (post mode keeps the "normal" branch),
 * BDS_CNV1_flag=0 (use the standard BDS path). */
PPPGlobal_t PPP_Glo = {0};
