#ifndef __FORWARD_H_INCLUDED
#define __FORWARD_H_INCLUDED

double _Complex impedance(double n0, double omega, double _Complex *n_arr,
                         double *depth_arr, int lay_num);

double _Complex u(double n0, double omega, double ver_dist, double *rho_arr,
                 double *depth_arr, int lay_num);

double _Complex spec_dens(double t, double *pars);

double _Complex H(double omega, double ver_dist, double hor_dist,
                 double *rho_arr, double *depth_arr, int lay_num);

double *forward_fun_fixed_net(double *freq_arr, int freq_num,
                              double ver_dist, double hor_dist,
                              double *rho_arr, int lay_num,
                              double first_thick, double step);

double *forward_fun_wrapper(double *x, double *pars, int x_dim, int par_dim);

#endif
