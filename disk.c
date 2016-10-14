#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#include "global_ex.h"

double v_K(double r){
    return 29.8*100000.0*sqrt(m_star)/sqrt(r);
}

double w_K(double r){
    return v_K(r)/(r*LUNIT);
}

double temperature (double r){
	r=r*LUNIT;
return pow(3.0,0.2)*pow(2.0,-1.4)*pow(M_PI,-0.4)\
	*pow(mu*m_p/gamma0/k_B,0.2)*pow(opa/sig_sb,0.2)*pow(alpha,-0.2)*pow(G*m_star*MUNIT,0.3)*pow((1-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.4)*pow(r,-0.9);
}

double sound_sp(double r) {
return sqrt(gamma0*k_B*temperature(r)/mu/m_p);
}

double density(double r) {
	r=r*LUNIT;
return pow(3.0,-1.3)*pow(2.0,1.6)*pow(M_PI,-0.9)\
	  *pow(mu*m_p/gamma0/k_B,1.2)*pow(opa/sig_sb,-0.3)*pow(alpha,-0.7)\
	  *pow(G*m_star*MUNIT,0.55)*pow((1-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.4)*pow(r,-1.65);
}
double mean_path(double r){

return mu*m_p/density(r)/2e-15;
}

double viscosity( double r){//http://www.ifu.ethz.ch/IE/education/AQAM/GasKinetics
	    //return 0.4991*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r)*100000.0;
	return 0.5*sqrt(8.0/gamma0/M_PI)*sound_sp(r)*mean_path(r);
}

double yeta(double r){
    return k_P*(sound_sp(r)/v_K(r))*(sound_sp(r)/v_K(r));
}

double vt_gas(double r){
        return v_K(r)*sqrt(1-yeta(r));
}

double Sigma (double r){
	r=r*LUNIT;
return pow(3.0,-1.2)*pow(2.0,1.4)*pow(M_PI,-0.6)\
	*pow(mu*m_p/gamma0/k_B,0.8)*pow(opa/sig_sb,-0.2)*pow(alpha,-0.8)\
	*pow(G*m_star*MUNIT,0.2)*pow((1.0-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.6)\
	*pow(r,-0.6); 
}

double height(double r){
	r=r*LUNIT;
return pow(3.0,0.1)*pow(2.0,-0.7)*pow(M_PI,-0.2)\
          *pow(mu*m_p/gamma0/k_B,-0.4)*pow(opa/sig_sb,0.1)*pow(alpha,-0.1)\
          *pow(G*m_star*MUNIT,-0.35)*pow((1-sqrt(r_star*LUNIT/r))*mdot*MUNIT/TUNIT,0.2)*pow(r,1.05); 
}

double vr_gas (double r){
        r=r*LUNIT;
return pow(3.0,1.2)*pow(2.0,-2.4)*pow(M_PI,-0.4)\
        *pow(mu*m_p/gamma0/k_B,-0.8)*pow(opa/sig_sb,0.2)*pow(alpha,0.8)\
        *pow(G*m_star*MUNIT,-0.2)*pow((1.0-sqrt(r_star*LUNIT/r)),-0.6)\
        *pow(mdot*MUNIT/TUNIT,0.4)*pow(r,-0.4);
}
