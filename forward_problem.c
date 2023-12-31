#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "math_utils.h"

#define mu0 (M_PI*.4e-6)

long double complex impedance(long double n0, long double omega, long double complex *n_arr, 
                         long double *depth_arr, int lay_num)
{
/*
 * This function computes impedance of layered model
 *Params:
 *        n0: long double, fourier parameter
 *        omega: long double, frequency
 *        n_arr: array of long double, wavenumbers for layers
 *        depth_arr: array of long double, layer thicknesses
 */

        int i;
        long double complex imp;
        long double complex n1, n2;

        imp = 1;
        n1 = n_arr[lay_num - 1];

        for(i = lay_num - 1; i > 0; i--){
                n2 = n1;
                n1 = n_arr[i-1];
                imp = ctanh(n1*depth_arr[i-1] + catanh(n1/n2*imp));
        } 
        return imp;
}

long double complex u(long double n0, long double omega, long double ver_dist, long double *rho_arr, 
                 long double *depth_arr, int lay_num)
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


        long double complex *n_arr;
        long double complex imp, n1;
        int i;

        n_arr = (long double complex *)malloc(lay_num*sizeof(long double complex));

        for(i = 0; i < lay_num; i++){
                n_arr[i] = csqrt(n0*n0 - I*omega*mu0/rho_arr[i]);
        }

        imp = impedance(n0, omega, n_arr, depth_arr, lay_num);
        n1 = n_arr[0];
        free(n_arr);

        return exp(-n0*ver_dist)*(n1 - n0*imp)/(2*(n1 + n0*imp));
}

long double complex spec_dens(long double t, long double *pars)
{
/*
 * This function wraps fourier potential and Bessel function into integrable
 *Params:
 *        t: point of evaluation
 *        pars: parameters of fourier potential
 */
        int i, lay_num;
        long double r, omega, ver_dist, bf;
        long double complex res;
        long double *rho_arr, *depth_arr;
        r = pars[0];
        omega = pars[1];
        ver_dist = pars[2];
        lay_num = (int)(pars[3]+0.1);
        rho_arr = (long double *)malloc(lay_num*sizeof(long double));
        depth_arr = (long double *)malloc((lay_num-1)*sizeof(long double));
        for(i = 0; i < lay_num; i++)
                rho_arr[i] = pars[4 + i];
        for(i = 1; i < lay_num; i++)
                depth_arr[i-1] = pars[3 + lay_num + i];

	bf = bessel_function(t*r);
	bf*=t;
	bf*=t;
      
        res = u(t, omega, ver_dist, rho_arr, depth_arr, lay_num);

        if(isnan(creal(res)))
		res = 0 + 0*I;

	res*=bf;
	//printf("%Lf", cimag(u(0.1, 1000, 100, rho_arr, depth_arr, lay_num)*1.E50 + 1*I));
	//res = cimag(csqrt(-1));
	//printf("%Lf", res);
        //printf("n0 = %f, sp_dens = %Lf + %Lfj \n", t, creal(res)*1.E10,
        //       cimag(res)*1.E10);
        //printf("n0 = %f, bf = %Lf\n", t, bessel_function(t*r));

        free(rho_arr);
        free(depth_arr);
        return res;
}

long double complex H(long double omega, long double ver_dist, long double hor_dist, 
                 long double *rho_arr, long double *depth_arr, int lay_num)
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

        long double complex int_result, a;
        long double *rho_reduced, *depth_reduced;
        long double *pars;
        long double result[2];
        long double n0;
        int i, j, lay_num_reduced;

        rho_reduced = (long double *)malloc(lay_num*sizeof(long double));
        depth_reduced = (long double *)malloc((lay_num-1)*sizeof(long double));

        rho_reduced[0] = rho_arr[0];
        depth_reduced[0] = depth_arr[0];

        j = 0;
        for(i = 1; i < lay_num; i++){
                if(fabs(rho_reduced[j] - rho_arr[i]) < 0.1){
                        if(i < lay_num - 1)
                                depth_reduced[j]+=depth_arr[i];
                } else {
                        j++;
                        rho_reduced[j] = rho_arr[i];
                        if(i < lay_num - 1)
                                depth_reduced[j] = depth_arr[i];
                }
        }

        lay_num_reduced = j+1;
        pars = (long double *)malloc((3 + 2*lay_num_reduced)*sizeof(long double));
        omega = omega*2*M_PI;
        pars[0] = hor_dist;
        pars[1] = omega;  
        pars[2] = ver_dist;
        pars[3] = lay_num_reduced;
        for(i = 0; i < lay_num_reduced; i++)
                pars[4 + i] = rho_reduced[i];
        for(i = 1; i < lay_num_reduced; i++)
                pars[3 + lay_num_reduced + i] = depth_reduced[i-1];

        //int_result = integral(spec_dens, pars, 3+2*lay_num_reduced,
        //                      0, 1, 0.001);
        int_result = 0;

        for(n0 = 0.001; n0 < 0.9; n0+=0.001){
		//printf("n0 = %Lf,    ", n0);
                a = spec_dens(n0, pars);
                //printf("a = %Lf \n", a);
                //printf("n0 = %f, %f, %f\n", n0, creal(a)*1.E10, cimag(a)*1.E10);
                int_result+=a*0.001;
        }
        result[0] = (0.5/M_PI)*10000*creal(int_result);
        result[1] = -(0.5/M_PI)*10000*cimag(int_result);
        free(rho_reduced);
        free(depth_reduced);
        return result[0] + result[1]*I;
}

long double *forward_fun_fixed_net(long double *freq_arr, int freq_num, 
                              long double ver_dist, long double hor_dist, 
                              long double *rho_arr, int lay_num, 
                              long double first_thick, long double step)
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
        long double *res, *depth_arr;
        long double thick, Ampl;
        long double complex res_fr;
        int i;

        res = (long double *)malloc(2*freq_num*sizeof(long double));
        depth_arr = (long double *)malloc((lay_num-1)*sizeof(long double));
        thick = first_thick;
        for(i = 1; i < lay_num; i++){
                depth_arr[i-1] = thick;
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
        free(depth_arr);
        return res;
}

long double *forward_fun_wrapper(long double *x, long double *pars)
{
/*
 * This function wraps forward fun into form convenient for gradient computation
 *Params:
 *        x: vector of log(resistivity)
 *        pars: frequencies and other known parameters
 */
        int freq_num, lay_num, i;
        long double ver_dist, hor_dist, first_thick, step;
        long double *freq_arr, *rho_arr, *result;

        freq_num = pars[0];
        lay_num = pars[1];

        ver_dist = pars[2];
        hor_dist = pars[3];
        first_thick = pars[4];
        step = pars[5];

        freq_arr = (long double *)malloc(freq_num*sizeof(long double));
        rho_arr = (long double *)malloc(lay_num*sizeof(long double));

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = pars[6 + i];
        for(i = 0; i < lay_num; i++)
                rho_arr[i] = exp(x[i]);
        
        result = forward_fun_fixed_net(freq_arr, freq_num, ver_dist, hor_dist,
                                     rho_arr, lay_num, first_thick, step);
        free(freq_arr);
        free(rho_arr);

        return result;
}

long double *forward_fun_inphase_diff(long double *x, long double *pars)
{
/*
 * This function computes the response for the case
 * when only DIFFERENCE between consecutive in-phase components
 * is available
 *Params:
 *        x: vector of log(resistivity)
 *        pars: frequencies and other known parameters
 */
        int freq_num, lay_num, i;
        long double ver_dist, hor_dist, first_thick, step;
        long double *freq_arr, *rho_arr, *result;

        freq_num = pars[0];
        lay_num = pars[1];

        ver_dist = pars[2];
        hor_dist = pars[3];
        first_thick = pars[4];
        step = pars[5];

        freq_arr = (long double *)malloc(freq_num*sizeof(long double));
        rho_arr = (long double *)malloc(lay_num*sizeof(long double));

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = pars[6 + i];
        for(i = 0; i < lay_num; i++)
                rho_arr[i] = exp(x[i]);
        
        result = forward_fun_fixed_net(freq_arr, freq_num, ver_dist, hor_dist,
                                     rho_arr, lay_num, first_thick, step);

	for(i = freq_num - 1; i > 0; i--){
                result[i] -= result[i-1];
	}
        free(freq_arr);
        free(rho_arr);

        return result;
}
