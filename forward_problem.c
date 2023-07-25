#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "math_utils.h"

#define mu0 (M_PI*.4e-6)

double complex impedance(double n0, double omega, double *n_arr, 
                         double *depth_arr, int lay_num)
{
/*
 * This function computes impedance of layered model
 *Params:
 *        n0: double, fourier parameter
 *        omega: double, frequency
 *        n_arr: array of double, wavenumbers for layers
 *        depth_arr: array of double, layer thicknesses
 */

        int i;
        double complex imp;
        double complex n1, n2;

        imp = 1;
        n1 = n_arr[lay_num - 1];

        for(i = lay_num - 1; i > 0; i --){
                n2 = n1;
                n1 = n_arr[i-1];
                imp = ctanh(n1*depth_arr[i-1] + catanh(n1/n2*imp));
        } 
        return imp;
}

double complex u(double n0, double omega, double ver_dist, double *rho_arr, 
                 double *depth_arr, int lay_num)
{
/*
 * This function computes fourier potential for magnetic field
 *Params:
 *        n0: base frequency
 *        omega: frequency, rad/sec
 *        ver_dist: sum of altitudes of receiver and emitter
 *        rho_arr: layer resistivities
 *        depth_arr: layer thicknesses
 *        lay_num: number of layers
 */


        double complex *n_arr;
        double complex imp, n1;
        int i;

        n_arr = (double complex *)malloc(lay_num*sizeof(double complex));

        for(i = 0; i < lay_num; i++){
                n_arr[i] = csqrt(n0*n0 - omega*mu0/rho_arr[i]);
        }

        imp = impedance(n0, omega, n_arr, depth_arr, lay_num);
        n1 = n_arr[1];

        return exp(-n0*ver_dist)*(n1 - n0*imp)/(2*(n1 + n0*imp));
}

double complex spec_dens(double t, double *pars)
{
/*
 * This function wraps fourier potential and Bessel function into integrable
 *Params:
 *        t: point of evaluation
 *        pars: parameters of fourier potential
 */
        return bessel_function(t*pars[0])*t*t*u(t, pars[1], pars[2], pars[3],
                               pars[4], (int)(pars[5] + 0.1));
}

double complex H(double omega, double ver_dist, double hor_dist, 
                 double *rho_arr, double *depth_arr, int lay_num)
{
/*
 * This function computes vertical component of magnetic field
 *Params:
 *        omega: frequency, Hz
 *        ver_dist: sum of altitudes of receiver and emitter
 *        hor_dist: horizontal distance between emitter and receiver
 *        rho_arr: layer resistivities
 *        depth_arr: layer thicknesses
 *        lay_num: number of layers
 */

        double complex int_result;
        double *rho_reduced, *depth_reduced;
        double pars[6];
        double result[2];
        int i, j, lay_num_reduced;

        rho_reduced = (double *)malloc(lay_num*sizeof(double));
        depth_reduced = (double *)malloc(lay_num*sizeof(double));

        rho_reduced[0] = rho_arr[0];
        depth_reduced[0] = depth_arr[0];

        j = 0;
        for(i = 1; i < lay_num; i ++){
                if(fabs(rho_reduced[j] - rho_arr[i]) < 0.01){
                        depth_reduced[j]+=depth_arr[i];
                } else {
                        j++;
                        rho_reduced[j] = rho_arr[i];
                        depth_reduced[j] = depth_arr[i];
                }
        }

        lay_num_reduced = j;
        omega = omega*2*M_PI;
        pars = {hor_dist, omega, ver_dist, rho_reduced,
                depth_reduced, (double)lay_num_reduced};

        int_results = integral(spec_dens, pars, 6, 0, 1, 0.001);

        result[0] = (0.5*M_PI)*int_result*10000*creal(int_result);
        result[1] = -(0.5*M_PI)*int_result*10000*cimag(int_result);

        return result[0] + result[1]*I;
}

double *forward_fun_fixed_net(double *freq_arr, int freq_num, 
                              double ver_dist, double hor_dist, 
                              double *rho_arr, int lay_num, 
                              double first_thick, double step)
{
/*
 * This function computes vertical component of magnetic field
 *Params:
 *        freq_arr: array of frequencies, Hz
 *        freq_num: number of frequencies
 *        ver_dist: sum of altitudes of receiver and emitter
 *        hor_dist: horizontal distance between emitter and receiver
 *        rho_arr: layer resistivities
 *        lay_num: number of layers
 *        first_thick: thickness of the first layer
 *        step: thickness of layer i = first_thick*(step^i)
 *Returns:
 *        vector of (all real component, all imaginary components)
 */
        double *res, *depth_arr;
        double thick;
        double complex res_fr;
        int i;

        res = (double *)malloc(2*freq_num*sizeof(double));
        depth_arr = (double *)malloc(lay_num*sizeof(double));
        thick = first_thick;
        for(i = 0; i < lay_num; i++){
                depth_arr[i] = thick;
                thick*=step;
        }

        Ampl = primal_field(28.8, 27.8);
        for(i = 0; i < freq_num; i++){
                res_fr = H(freq_arr[i], ver_dist, hor_dist,
                         rho_arr, depth_arr, lay_num);
                res_fr/=Ampl;
                res[i] = creal(res_fr);
                res[freq_num + i] = cimag(res_fr);
        }
        return res;
}

double *forward_fun_wrapper(double *x, double *pars, int x_dim, int par_dim)
{
/*
 * This function wraps forward fun into form convenient for gradient computation
 *Params:
 *        x: vector of log(resistivity)
 *        pars: frequencies and other known parameters
 *        *_dim: dimensionalities 
 */
        int freq_num, lay_num, i;
        double ver_dist, hor_dist, first_thick, step;
        double *freq_arr, *res_arr;

        freq_num = pars[0];
        lay_num = pars[1];

        ver_dist = pars[2];
        hor_dist = pars[3];
        first_thick = pars[4];
        step = pars[5];

        freq_arr = (double *)malloc(freq_num*sizeof(double));
        res_arr = (double *)malloc(lay_num*sizeof(double));

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = pars[5 + i];
        for(i = 0; i < lay_num; i++)
                res_arr[i] = exp(pars[5 + freq_num + i]);

        return forward_fun_fixed_net(freq_arr, freq_num, ver_dist, hor_dist,
                                     rho_arr, lay_num, first_thick, step);
}
