#ifndef __FORWARD_H_INCLUDED
#define __FORWARD_H_INCLUDED

long double _Complex impedance(long double n0, long double omega,
                               long double _Complex *n_arr,
                               long double *depth_arr, int lay_num);

long double _Complex u(long double n0, long double omega, long double ver_dist,
                       long double *rho_arr,
                       long double *depth_arr, int lay_num);

long double _Complex spec_dens(long double t, long double *pars);

long double _Complex H(long double omega, long double ver_dist,
                       long double hor_dist, long double *rho_arr,
                       long double *depth_arr, int lay_num);

long double *forward_fun_fixed_net(long double *freq_arr, int freq_num,
                              long double ver_dist, long double hor_dist,
                              long double *rho_arr, int lay_num,
                              long double first_thick, long double step);

long double *forward_fun_wrapper(long double *x, long double *pars);

#endif
