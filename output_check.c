#include <stdio.h>
#include "ex_func.h"
#include "global_ex.h"

double vr_p[2]={0.0};
void check_disk(double r){
	printf("Properties at 1 AU:mdot=%e alpha=%f\n",mdot,alpha);
	printf("Surf_dens=%g\nTemperature=%f\nMidplane-density=%g\nSound_speed=%g\n\
H/R=%g\nMean free path=%g\nv_K=%g\nvt_gas=%g\nvr_peb=%g\nzeta=%g\n\
Omega_K=%g\nvr_gas=%g\n",\
Sigma(r),temperature(r),density(r)\
,sound_sp(r),height(r)/r/LUNIT,mean_path(r),v_K(r),\
vt_gas(r),vr_estimate(1.0,1.0,vr_p),yeta(r),w_K(r),vr_gas(r));
}
