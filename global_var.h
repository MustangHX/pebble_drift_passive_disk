#include "global_ex.h"
#ifndef PEB_STRUCT
#define PEB_STRUCT
extern double dt_ring[ring_num];
typedef struct PEBBLE{
	int max_step;
	double mass[time_yr];
        double rad[time_yr];
        double size[time_yr];
        double time[time_yr];
        double vr[time_yr];
} PEBBLE;
extern PEBBLE peb_group[peb_num];

typedef struct PEBBLE_MAP{
        double dr;
	double rad;//inner edge radius
	double rad_med;//middle radius
	double AREA;//ring AREA
	double time;
        double size[peb_size_num+1];
	double size_med[peb_size_num];
	double mass_in[peb_size_num];
	double mass_out[peb_size_num];
        double surf_dens[peb_size_num];
	double rho[peb_size_num];
	double vr[peb_size_num];//vr when both r and size are in center of box
	double vt[peb_size_num];
	double vr_med_r[peb_size_num+1];//vr when r=r_med
	double vr_med_s[peb_size_num];//vr when a_pb=a_pb_med at inner edge radius
	double vt_med_r[peb_size_num+1];
	double hei[peb_size_num];
	double flux[peb_size_num];
} PEBBLE_MAP;

typedef struct DUST_MAP{
        double dr;
        double rad;//inner edge radius
	double rad_med;
	double AREA;
	double mass_in;
	double mass_out;
        double surf_dens;
	double rho;
} DUST_MAP;
extern PEBBLE peb_group[peb_num];
extern PEBBLE_MAP peb_map[ring_num];
extern DUST_MAP dust_budget[ring_num];
extern double alpha;
extern double mdot;
#endif
