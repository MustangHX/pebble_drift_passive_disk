#include <stdio.h>
#include "ex_func.h"
#include "global_ex.h"

double vr_p[2]={0.0};
void check_disk(double r){
	FILE *fp1,*fp2,*fp3,*fp4;
	int i;
	ITER=0;
	opa=10.0;
	printf("Properties at 1 AU:mdot=%e alpha=%f\n",mdot,alpha);
	printf("Surf_dens=%g\nTemperature=%f\nMidplane-density=%g\nSound_speed=%g\n\
H/R=%g\nMean free path=%g\nv_K=%g\nvt_gas=%g\nvr_peb=%g\nzeta=%g\n\
Omega_K=%g\nvr_gas=%g\n",\
Sigma(r),temperature(r),density(r)\
,sound_sp(r),height(r)/r/LUNIT,mean_path(r),v_K(r),\
vt_gas(r),vr_estimate(1.0,1.0,vr_p),yeta(r),w_K(r),vr_gas(r));

	fp1=fopen("opacity.txt","w");
	fp2=fopen("sigma.txt","w");
	fp3=fopen("temperature.txt","w");
	fp4=fopen("dust_mass.txt","w");
	for(i=0;i<ring_num;i++){
	fprintf(fp1,"%e\t%e\n",dust_budget[i].rad,func_spline3(dust_budget[i].rad,p_opa_line));
	fprintf(fp2,"%e\t%e\n",dust_budget[i].rad,Sigma(dust_budget[i].rad));
	fprintf(fp3,"%e\t%e\n",dust_budget[i].rad,temperature(dust_budget[i].rad));
	fprintf(fp4,"%e\t%e\n",dust_budget[i].rad,dust_budget[i].mass_out);

	}
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);}

