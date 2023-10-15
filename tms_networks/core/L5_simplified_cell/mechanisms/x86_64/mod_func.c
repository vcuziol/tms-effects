#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _fdsexp2s_reg(void);
extern void _hh3_reg(void);
extern void _HH_traub_reg(void);
extern void _IM_cortex_reg(void);
extern void _kdr2_reg(void);
extern void _tms2_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," fdsexp2s.mod");
    fprintf(stderr," hh3.mod");
    fprintf(stderr," HH_traub.mod");
    fprintf(stderr," IM_cortex.mod");
    fprintf(stderr," kdr2.mod");
    fprintf(stderr," tms2.mod");
    fprintf(stderr, "\n");
  }
  _fdsexp2s_reg();
  _hh3_reg();
  _HH_traub_reg();
  _IM_cortex_reg();
  _kdr2_reg();
  _tms2_reg();
}
