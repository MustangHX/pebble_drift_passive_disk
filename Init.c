#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double pp_vr_tau1[2]={0.0};

void Init(){
	int i,j;
	double AREA;
/*        for(i=0;i<100;i++){
                size=drag_group(i*1.0+0.5,0.01);
                for(j=0;j<100;j++){		
			if(j==(int)(size*10)){
                               pebble[i][j]=i*1.0*2*M_PI*dr;
		                             }
		                  }
		          }
	fp=fopen("velocity.dat","r");
	if(fp=None){
		for(j=0;j<100;j++){
			size=j*1.0;
			for(i=0;i<100;i++){
				velocity[i][j]
			
}
}
}*/
	
	for(i=0;i<ring_num;i++){
		dust_budget[i].rad=i*size_ring+1.0-size_ring;
		AREA=M_PI*((dust_budget[i].rad+size_ring/2.0)*(dust_budget[i].rad+size_ring/2.0)-(dust_budget[i].rad-size_ring/2.0)*(dust_budget[i].rad-size_ring/2.0))*LUNIT*LUNIT;
		dust_budget[i].AREA=AREA;
		dust_budget[i].dr=size_ring;
		dust_budget[i].mass_in=0.0;
		for(j=0;j<peb_size_num;j++){
                dust_budget[i].surf_dens=Sigma(dust_budget[i].rad)*dust_gas;
		dust_budget[i].rho=density(dust_budget[i].rad)*dust_gas;
		}
		dust_budget[i].mass_out=dust_budget[i].surf_dens*AREA;
	}
	i=0;
	for(i=0;i<ring_num;i++){
  //              printf("SIGMAAAAA%e\t%d\t",peb_map[i].surf_dens,i);
//		printf("%d\n",i);
		peb_map[i].dr=size_ring;
		peb_map[i].rad=i*size_ring+1.0-size_ring;
                peb_map[i].time=0.0;
		peb_map[i].AREA=M_PI*((peb_map[i].rad+size_ring/2.0)*(peb_map[i].rad+size_ring/2.0)-(peb_map[i].rad-size_ring/2.0)*(peb_map[i].rad-size_ring/2.0))*LUNIT*LUNIT;
		for(j=0;j<=peb_size_num;j++){
			peb_map[i].size[j]=0.1*pow(10,j*size_step);
		}
		for(j=0;j<peb_size_num;j++){
			peb_map[i].size_med[j]=0.5*(peb_map[i].size[j]+peb_map[i].size[j+1]);
			if(i==1) printf("SIZE=%fcm\n",peb_map[i].size[j]);
				//drag_group((i+1)*0.25,0.01*(j+1));
			AREA=M_PI*((peb_map[i].rad+size_ring/2.0)*(peb_map[i].rad+size_ring/2.0)-(peb_map[i].rad-size_ring/2.0)*(peb_map[i].rad-size_ring/2.0))*LUNIT*LUNIT;
			if ((j<1 && i==ring_num-1)|| 1) {
				//peb_map[i].mass_out[j]=0.1*AREA*0.01*(Sigma((i+4)*0.25)*exp(-1.0*peb_map[i].size[j]/0.1)+1e-10);
				peb_map[i].mass_out[j]=1.0*AREA*(0.01*Sigma(peb_map[i].rad_med)*exp(-1.0*peb_map[i].size_med[j]/0.1)+1e-10);
			}
			else peb_map[i].mass_out[j]=peb_low_lim*AREA;
			peb_map[i].mass_out[j]/=peb_size_num;//comes from dust_budget
			peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
			peb_map[i].mass_in[j]=0.0;
		}
  //              printf("SIGMAAAAA%e\t%d\n",peb_map[i].surf_dens,i);
	}
	for(i=0;i<ring_num;i++){
	for(j=0;j<peb_size_num;j++){
//		printf("%d%d %e\t",i,j,peb_map[i].surf_dens[j]);
	}
//	printf("\n");
	}
}

void Init2(){// disk with variable resolution
	int i,j,jj;
	double tau,AREA,size_ring1,size_ring2,rad1,rad2,mass_norm=0.0;//size1 is smaller ring
	double v1,v2,delta_v,Re,vr_temp;
	FILE *fp;
	size_ring1=size_ring/10.0;
        size_ring2=size_ring*10.0;
			
	for(i=0;i<ring_num;i++){
	if(i<i_lim1) {
		dust_budget[i].rad=i*size_ring1+0.1;
		dust_budget[i].dr=size_ring1;
	}
	if(i>=i_lim1 && i < i_lim2){
		dust_budget[i].rad=(i-i_lim1)*size_ring+i_lim1*size_ring1+0.1;
		dust_budget[i].dr=size_ring;
	}
	if(i>=i_lim2 && i < ring_num){
		dust_budget[i].rad=(i-i_lim2)*size_ring2+30.0;
		dust_budget[i].dr=size_ring2;
	}
	
	printf("ring#=%d\tRAD=%f\n",i,dust_budget[i].rad);
	}
	printf("ring_size:%f\t%f\t%f\n",size_ring1,size_ring,size_ring2);
	for(i=0;i<ring_num;i++){
	rad1=dust_budget[i].rad;
	if(i+1<ring_num) rad2=dust_budget[i+1].rad;
	else rad2=R_OUT;
	dust_budget[i].AREA=M_PI*(rad2*rad2-rad1*rad1)*LUNIT*LUNIT;
//	printf("ring#=%d AREA %g %g\n",i,M_PI*(rad2*rad2-rad1*rad1)*LUNIT*LUNIT,dust_budget[i].AREA);
	dust_budget[i].rad_med=0.5*(rad1+rad2);
//	printf("ring#=%d\tRAD_MED=%f\n",i,dust_budget[i].rad_med);
	dust_budget[i].mass_in=0.0;
	for(j=0;j<peb_size_num;j++){
	dust_budget[i].surf_dens=Sigma(dust_budget[i].rad_med)*dust_gas;
	dust_budget[i].rho=density(dust_budget[i].rad_med)*dust_gas;
	}
	dust_budget[i].mass_out=dust_budget[i].surf_dens*dust_budget[i].AREA;
	//printf("%d\t%g\n",i,dust_budget[i].surf_dens);
        printf("ring#=%d\tRAD_MED=%f %g %g\n",i,dust_budget[i].rad_med,M_PI*(rad2*rad2-rad1*rad1)*LUNIT*LUNIT,dust_budget[i].AREA);

	}
	for(i=0;i<ring_num;i++){
		mass_norm=0.0;
	peb_map[i].rad=dust_budget[i].rad; 
	peb_map[i].rad_med=dust_budget[i].rad_med;
	peb_map[i].dr=dust_budget[i].dr;
	peb_map[i].AREA=dust_budget[i].AREA;

	printf("ring#=%d\tPEB_RAD_MED=%f\n",i,peb_map[i].rad_med);

	for(j=0;j<=peb_size_num;j++){
		peb_map[i].size[j]=size_min*pow(10,j*size_step);
	}
	for(j=0;j<peb_size_num;j++){
                peb_map[i].size_med[j]=0.5*(peb_map[i].size[j]+peb_map[i].size[j+1]);
		//peb_map[i].vr[j]=vr_estimate(peb_map[i].rad_med,peb_map[i].size_med[j],pp_vr_tau1);
		peb_map[i].vr[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size_med[j],pp_vr_tau1);
		tau=pp_vr_tau1[1];
		peb_map[i].vt[j]=0.5*tau*peb_map[i].vr[j];
		peb_map[i].vr_med_s[j]=drift_vr(peb_map[i].rad,peb_map[i].size_med[j],pp_vr_tau1);
                //tau=pp_vr_tau1[1];
                //peb_map[i].vt_med_s[j]=0.5*tau*peb_map[i].vr_mid[j];
	//	if(i==1) printf("SIZE=%fcm\n",peb_map[i].size[j]);
        	peb_map[i].vr_med_r[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size[j],pp_vr_tau1);
                tau=pp_vr_tau1[1];
                peb_map[i].vt_med_r[j]=0.5*tau*peb_map[i].vr_med_r[j];
		peb_map[i].hei[j]=height(peb_map[i].rad_med)/sqrt(1+tau/alpha);
		//if(i<=151 && i > 145 && j < 10) printf("hei=%e tau= %e\n",peb_map[i].hei[j],tau);
        
	}
	j=peb_size_num;
	peb_map[i].vr_med_r[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size[j],pp_vr_tau1);
        tau=pp_vr_tau1[1];
        peb_map[i].vt_med_r[j]=0.5*tau*peb_map[i].vr_med_r[j];
	
	for(j=0;j<peb_size_num;j++) {
		if(peb_map[i].size_med[j]<peb_size_lim){
			mass_norm+=exp(-1.0*size_slope*peb_map[i].size_med[j]);
		}
	}
        for(j=0;j<peb_size_num;j++){
		AREA=peb_map[i].AREA;
		if (peb_map[i].size_med[j]<peb_size_lim && i > 4) {
			peb_map[i].mass_out[j]=peb_dust*AREA*(dust_budget[i].surf_dens*exp(-1.0*size_slope*peb_map[i].size_med[j])/mass_norm+1e-10);
			//peb_map[i].mass_out[j]=peb_dust*AREA*(0.1*Sigma(peb_map[i].rad_med)*exp(-1.0*peb_map[i].size_med[j])/mass_norm+1e-10);
//			peb_map[i].mass_out[j]=0.1*AREA*(dust_budget[i].surf_dens+1e-10);
		}
		else peb_map[i].mass_out[j]=peb_low_lim*AREA;
		//peb_map[i].mass_out[j]/=peb_size_num;//comes from dust_budget
		peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/AREA;
		
		peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];

		if(1 && 0){		
		vr_temp=drift_vr(peb_map[i].rad_med,peb_map[i].size[j],pp_vr_tau1);
		tau=pp_vr_tau1[1];				     
		peb_map[i].hei[j]=height(peb_map[i].rad_med)/sqrt(1+tau/alpha);
		peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];///2e9;

		}
		else	peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];

		//		printf("rho=%e\t%e\t%d\t%d\t",peb_map[i].rho[j],peb_map[i].rad_med,i,j);
		peb_map[i].mass_in[j]=0.0;
		if(i<=151 && i > 145 && j < 3) printf("INIT SURF_DENS=%e MID_DENS=%e HEI=%e SIZE=%d\n",peb_map[i].surf_dens[j],peb_map[i].rho[j],peb_map[i].hei[j],j);
	}
	}
	i=40;
	fp=fopen("relative_velocity.txt","w");
	for(j=0;j<peb_size_num;j++){
	for(jj=0;jj<peb_size_num;jj++){
		v1=fabs(peb_map[i].vr[j]-peb_map[i].vr[jj]);
		v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
		delta_v=sqrt(v1*v1+v2*v2);
		fprintf(fp,"%e\t",v1);
	}
	fprintf(fp,"\n");
	}
	fclose(fp);
	fp=fopen("Reynolds2d.txt","w");
	for(i=0;i<ring_num;i++){
	for(j=0;j<peb_size_num;j++){
	v1=peb_map[i].vr[j];
	v2=peb_map[i].vt[j];
	Re=2.0*peb_map[i].size_med[j]*sqrt(v1*v1+v2*v2)/viscosity(peb_map[i].rad_med);
	fprintf(fp,"%e\t",Re);
	}
	fprintf(fp,"\n");
	}
	fclose(fp);
}

void disk_evolve(){// evolving disk, update mdot, alpha, drift_velocity
	int i,j;
	double tau;//size1 is smaller ring
	
	
	for(i=0;i<ring_num;i++){

	for(j=0;j<peb_size_num;j++){
		peb_map[i].vr[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size_med[j],pp_vr_tau1);
		tau=pp_vr_tau1[1];
		peb_map[i].hei[j]=height(peb_map[i].rad_med)/sqrt(1+tau/alpha);

		peb_map[i].vt[j]=0.5*tau*peb_map[i].vr[j];
		peb_map[i].vr_med_s[j]=drift_vr(peb_map[i].rad,peb_map[i].size_med[j],pp_vr_tau1);
        	peb_map[i].vr_med_r[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size[j],pp_vr_tau1);
                tau=pp_vr_tau1[1];
                peb_map[i].vt_med_r[j]=0.5*tau*peb_map[i].vr_med_r[j];
	if(i<=151 && i > 145 && j < 3) printf("hei=%e tau= %e\n",peb_map[i].hei[j],tau);
        
	}
	j=peb_size_num;
	peb_map[i].vr_med_r[j]=drift_vr(peb_map[i].rad_med,peb_map[i].size[j],pp_vr_tau1);
        tau=pp_vr_tau1[1];
        peb_map[i].vt_med_r[j]=0.5*tau*peb_map[i].vr_med_r[j];
	
	}

}


void Restart(int rnum){
	double AREA,dens,dust,rad;
	int i,j;
	FILE *fp,*fp_dust;
	char name[256],name1[256];
	printf("RESTART=%d\n",rnum);
        sprintf(name,"out_sigma%d.txt",rnum);
	sprintf(name1,"dust_sigma%d.txt",rnum);
	printf("%s\n%s\n",name,name1);
	fp=fopen(name,"r");
	fp_dust=fopen(name1,"r");
        for(i=0;i<ring_num;i++){
//AREA=M_PI*((peb_map[i].rad+size_ring/2.0)*(peb_map[i].rad+size_ring/2.0)-(peb_map[i].rad-size_ring/2.0)*(peb_map[i].rad-size_ring/2.0))*LUNIT*LUNIT;
	AREA=dust_budget[i].AREA;
	fscanf(fp_dust,"%lf%lf",&rad,&dust);
	dust_budget[i].surf_dens=dust;
	dust_budget[i].rho=dust/height(dust_budget[i].rad_med)/sqrt(2.0*M_PI);
	dust_budget[i].mass_out=dust*AREA;
	dust_budget[i].mass_in=0.0;
        for(j=0;j<peb_size_num;j++){
		fscanf(fp,"%lf",&dens);
		peb_map[i].surf_dens[j]=dens;
		peb_map[i].mass_out[j]=AREA*peb_map[i].surf_dens[j];
		peb_map[i].mass_in[j]=0.0;
	}
	}
}
