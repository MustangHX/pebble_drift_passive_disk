#include "ex_func.h"
#include "global_ex.h"
#include "global_var.h"
#include<math.h>
#include<stdio.h>
void frag(){
int i,j,jl,js;
double size_lim,frac,AREA,dr,delta_mass;
dr=size_ring;
for(i=ring_num-1;i>=0;i--){
	dr=peb_map[i].dr;
	for(j=peb_size_num-1;j>=0;j--){
		size_lim=2.25*mean_path(peb_map[i].rad);		
		if(peb_map[i].size[j] > 0.8*size_lim){
			jl=j;
			js=floor(log10(peb_map[i].size[jl]/0.1/2.0)/size_step);
//			js=jl;
			if(js <0) js=0;
			frac=peb_map[i].size[jl]/size_lim;
			frac=frac*0.9;
			if (frac >1.0) frac=1.0;
			//frac=0.0;
			delta_mass=peb_map[i].mass_out[j]*frac;
			peb_map[i].mass_out[jl]-=delta_mass;
			peb_map[i].mass_out[js]+=delta_mass;
	}
}
}

for(i=0;i<ring_num;i++){
	AREA=M_PI*((peb_map[i].rad+dr/2.0)*(peb_map[i].rad+dr/2.0)-(peb_map[i].rad-dr/2.0)*(peb_map[i].rad-dr/2.0))*LUNIT*LUNIT;
	for(j=0;j<peb_size_num;j++){
		peb_map[i].surf_dens[j]=(peb_map[i].mass_out[j]+peb_map[i].mass_in[j])/AREA;
	}
}
}

double coag_kernel(double a_pb1,double a_pb2,double delta_v,double n_density1, double n_density2,double hei1, double hei2, int i, int j, int jj){
	double cross_section, path,n_density;
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*n_density1*n_density2/sqrt(2*M_PI*(hei1*hei1+hei2*hei2));
	if(delta_v>v_coag_max+v_tran_width) path=0.0;
	else if(delta_v>v_coag_max && delta_v<=v_tran_width) path\
		=(v_coag_max+v_tran_width-delta_v)/v_tran_width*path;
	if(path>1.0)	printf("path=%e N1=%e  N2=%e  sig1=%e hei=%e %d %d %d\n",\
			path,n_density1,n_density2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj);
	return 1.0*path;
}


double coag_kernel_sweep(double a_pb1,double a_pb2,double delta_v,double n_density1, double n_density2,double hei1, double hei2, int i, int j, int jj){
	double cross_section, path,n_density;
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*n_density1*n_density2/sqrt(2*M_PI*(hei1*hei1+hei2*hei2));
	if(path>1.0)	printf("path=%e N1=%e  N2=%e  sig1=%e hei=%e %d %d %d\n",\
			path,n_density1,n_density2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj);
	return 1.0*path;
}


double coag_kernel2(double a_pb1,double a_pb2,double delta_v,double rho1, double rho2, double dt0, double rad, int i, int j, int jj){
	double cross_section, path,n_density,n_density1,n_density2;
	n_density2=rho2/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0)*0.85;
	n_density1=rho1/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0)*0.85;//disk within -h to h
	if(a_pb1 >= a_pb2) n_density=n_density2;
	else n_density=n_density1;
	cross_section=M_PI*(a_pb1*a_pb1+a_pb2*a_pb2);
	cross_section=M_PI*(a_pb1+a_pb2)*(a_pb1+a_pb2);
	//delta_v=fabs(v1-v2);
	path=cross_section*delta_v*n_density;
	if(path>1.0)	printf("path=%e N1=%e rho1=%e N2=%e rho2=%e sig1=%e hei=%e %d %d %d\n",\
			path,n_density1,rho1,n_density2,rho2,peb_map[i].surf_dens[jj],peb_map[i].hei[jj],i,j,jj);
	return 1.0*path;
}

void coagulation(double dt0){
        int i,j,jj,jjj=0,j_new,N_diff=5;
	double a_pb1,a_pb2,a_pb3,v1,v2,delta_v,dr,dt1,sub_t;
	double frac_s,frac_s2,mass0,n1,n2,n0,n_delta,h1,h2,m_pb1,m_pb2,m_pb3,AREA;
	
	for (i=ring_num-1;i>=1;i--){
		dr=peb_map[i].dr;
		AREA=peb_map[i].AREA;
		frac_s=0.0;
		dt1=2.0*dt0;
		do{
		dt1=dt1/2.0;
		frac_s=0.0;
		for(j=peb_size_num-1;j>=1;j--){
			if(peb_map[i].surf_dens[j] < 1e-5) continue;
			for(jj=j-1;jj>=0;jj--){
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];
			n1=peb_map[i].surf_dens[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			n2=peb_map[i].surf_dens[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
                        v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
			if(0||jj>((j-N_diff)>0?(j-N_diff):0)){
			if(n1<n2) frac_s2=n_delta/n1;//(2*M_PI*h1*h2*n2/sqrt(2*M_PI*(h1*h1+h2*h2)));
			else frac_s2=n_delta/n2;//(2*M_PI*h1*h2*n1/sqrt(2*M_PI*(h1*h1+h2*h2)));
			}
			else if(1){
				frac_s2=n_delta/n2;
				if(frac_s2<pow(m_pb2*n_delta/n1/m_pb1+1,1.0/3.0)/a_pb1-1.0){
					frac_s2=pow(m_pb2*n_delta/n1/m_pb1+1,1.0/3.0)-1.0;
				}
	//			if(frac_s2>0.1) printf("AHAGAGA Frac_S=%e i=%d j=%d jj=%d\n",frac_s2,i,j,jj);
			}



			//frac_s2=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/n2*h2;
			//else frac_s2=n_delta/n1*h1;
			
			
			if(frac_s2>frac_s){
			frac_s=frac_s2;
			jjj=j;
			}
			
			}
			}
			
			}while(frac_s>0.02);
	if( 0 && dt1 < dt0) printf("ring=%d\tdt=%g %f max_frac_s=%f j=%d\n",i,dt1,dt0,frac_s,jjj);
	sub_t=0.0;
	while(sub_t<dt0){
		for(j=peb_size_num-1;j>=1;j--){
                        if(peb_map[i].surf_dens[j] < 1e-5) continue;
                        for(jj=j-1;jj>=0;jj--){
			
				
			a_pb1=peb_map[i].size_med[j];
			a_pb2=peb_map[i].size_med[jj];

			a_pb3=pow(a_pb1*a_pb1*a_pb1+a_pb2*a_pb2*a_pb2,1.0/3.0);
			j_new=floor(log10(a_pb3/size_min)/size_step);
			if(j_new>=peb_size_num) j_new=peb_size_num-1;
			
			m_pb1=rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0;
			m_pb2=rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0;
			m_pb3=rho_peb*4.0*M_PI*a_pb3*a_pb3*a_pb3/3.0;
			n1=peb_map[i].surf_dens[j]/m_pb1;
			n2=peb_map[i].surf_dens[jj]/m_pb2;

			//n1=peb_map[i].rho[j]/(rho_peb*4.0*M_PI*a_pb1*a_pb1*a_pb1/3.0);
			//n2=peb_map[i].rho[jj]/(rho_peb*4.0*M_PI*a_pb2*a_pb2*a_pb2/3.0);
			
			v1=peb_map[i].vr[j]-peb_map[i].vr[jj];
                        v2=peb_map[i].vt[j]-peb_map[i].vt[jj];
			h1=peb_map[i].hei[j];
			h2=peb_map[i].hei[jj];
			delta_v=sqrt(v1*v1+v2*v2);
			n_delta=dt1*TUNIT*coag_kernel(a_pb1,a_pb2,delta_v,n1,n2,h1,h2,i,j,jj);
			
			if(0||jj>((j-N_diff)>0?(j-N_diff):0)){
                        	peb_map[i].mass_in[j_new]+=m_pb3*n_delta*AREA;
                        	peb_map[i].mass_out[jj]-=m_pb2*n_delta*AREA;
                        	peb_map[i].mass_out[j]-=m_pb1*n_delta*AREA;			
			}
			else if(0){	
				j_new=j+1;
				if(j_new>peb_size_num-1) j_new=j;
				frac_s=pow(m_pb2*n_delta/n1/m_pb1+1,1.0/3.0)-1.0;
				peb_map[i].mass_in[jj]-=m_pb2*n_delta*AREA;
				peb_map[i].mass_in[j]-=frac_s*peb_map[i].mass_out[j];
				peb_map[i].mass_in[j_new]+=m_pb2*n_delta*AREA+frac_s*peb_map[i].mass_out[j];
				
			}

			//frac_s=dt1*TUNIT*coag_kernel2(a_pb1,a_pb2,delta_v,peb_map[i].rho[j],peb_map[i].rho[jj],dt1,peb_map[i].rad_med,i,j,jj);
			//if(n1<n2) frac_s2=n_delta/(2*M_PI*h1*h2*n2/sqrt(2*M_PI(h1*h1+h2*h2)));
			//else frac_s2=n_delta/(2*M_PI*h1*h2*n1/sqrt(2*M_PI(h1*h1+h2*h2)));
			
			
			//		if(frac_s> 0.8) printf("rho=%e\t%e\t%f\t%f\t%f\t%f\t%f\n",peb_map[i].rho[j],peb_map[i].rho[jj],frac_s,peb_map[i].rad_med,peb_map[i].size_med[j],peb_map[i].size_med[jj],delta_v);
			mass0=peb_map[i].mass_out[j];
//			peb_map[i].mass_in[j_new]+=mass0*frac_s*(pow(a_pb2/a_pb1,3)+1.0);
//			peb_map[i].mass_out[jj]-=mass0*frac_s*pow(a_pb2/a_pb1,3);
//			peb_map[i].mass_out[j]-=mass0*frac_s;
			
			peb_map[i].mass_in[j_new]+=m_pb3*n_delta*AREA;
			peb_map[i].mass_out[jj]-=m_pb2*n_delta*AREA;
			peb_map[i].mass_out[j]-=m_pb1*n_delta*AREA;
			}
			
		}
		for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        	}
		
	sub_t+=dt1;

	}
	}

        for(i=ring_num-1;i>=0;i--){
        for(j=0;j<peb_size_num;j++){
                peb_map[i].mass_out[j]+=peb_map[i].mass_in[j];
                peb_map[i].mass_in[j]=0.0;
                peb_map[i].surf_dens[j]=peb_map[i].mass_out[j]/peb_map[i].AREA;
                peb_map[i].rho[j]=peb_map[i].surf_dens[j]/sqrt(2.0*M_PI)/peb_map[i].hei[j];
        }
        }


}
