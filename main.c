#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "math_utils.h"
#include "forward_problem.h"
#include "kalman.h"
#include "data_processing.h"

int main()
{

        double args[19] = {15, 25, 4, 1.1085, 77, 231, 385, 540, 694, 848,
                           1003, 1466, 1620, 1775, 2855, 4244, 5324, 9799,
                           15046};
        data_inversion("input", "output", args);
        return 0;



        /* 
        double *rho_arr, *freq_arr, *resp, *res, *observed, *x0;
        double *lower_bounds, *upper_bounds;
        double **P0, **Q, **R;
        double complex imp;
        int i, j, lay_num, freq_num;
        lay_num = 2;
        freq_num = 2;

        rho_arr = (double *)malloc(lay_num*sizeof(double));
        freq_arr = (double *)malloc(freq_num*sizeof(double));
        lower_bounds = (double *)malloc(lay_num*sizeof(double));
        upper_bounds = (double *)malloc(lay_num*sizeof(double));
        observed = (double *)malloc((6 + freq_num)*sizeof(double));
        x0 = (double *)malloc(lay_num*sizeof(double));
        P0 = (double **)malloc(lay_num*sizeof(double *));
        Q = (double **)malloc(lay_num*sizeof(double *));
        R = (double **)malloc(2*freq_num*sizeof(double *));

        for(i = 0; i < lay_num; i++){
                rho_arr[i] = 100.+50*i;
                x0[i] = log(200. - 50*i);
                P0[i] = (double *)malloc(lay_num*sizeof(double));
                Q[i] = (double *)malloc(lay_num*sizeof(double));
                for(j = 0; j < lay_num; j++){
                        if(i==j){
                               P0[i][j] = 0.09;
                        } else{
                               P0[i][j] = 0.;
                        }
                        Q[i][j] = 0.;
                }
        }

        for(i = 0; i < 2*freq_num; i++){
                R[i] = (double *)malloc(2*freq_num*sizeof(double));
                for(j = 0; j < 2*freq_num; j++){
                        if(i==j){
                               R[i][j] = 0.5;
                        } else{
                               R[i][j] = 0.;
                        }
                }
        }

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = 100. + 100*i;
        for(i = 0; i < lay_num; i++){
                lower_bounds[i] = log(0.1);
                upper_bounds[i] = log(3000.);
        }
        //imp = H(100, 100, 50, rho_arr, depth_arr, 3)*1.E7;
        //imp = u(0.01, 100, 100, rho_arr, depth_arr, 3)*1.E7;
        //imp = spec_dens(0.001, pars);
        resp = forward_fun_fixed_net(freq_arr, freq_num, 100., 50.,
                                    rho_arr, lay_num, 30, 1.1);
        
        for(i = 0; i < 2*freq_num; i++)
                printf("%f\n", resp[i]);

        observed[0] = freq_num;
        observed[1] = lay_num;
        observed[2] = 100.;  
        observed[3] = 50.;
        observed[4] = 30;
        observed[5] = 1.1;
        for(i = 0; i < freq_num; i++)
                observed[6+i] = freq_arr[i];

        res = kalman_step(resp, observed, x0, 2*freq_num, 6+freq_num, lay_num,
                          forward_fun_wrapper,
                          P0, Q, R, lower_bounds, upper_bounds, 8, 1.);

        for(i = 0; i < lay_num; i++)
                printf("result: %f\n", exp(res[i]));
        for(i = lay_num; i < 2*lay_num+1; i++)
                printf("result: %f\n", res[i]);


        free(lower_bounds);
        free(upper_bounds);
        free(P0);
        free(Q);
        free(R);
        free(x0);
        free(freq_arr);
        free(rho_arr);
        free(res);
        free(resp);
        free(observed);
        return 0;
        */
        /*double *rho_arr, *depth_arr;
        double complex res;

        rho_arr = (double *)malloc(2*sizeof(double));
        depth_arr = (double *)malloc(1*sizeof(double));
        rho_arr[0] = 1150;
        rho_arr[1] = 100;
        depth_arr[0] = 10;

        res = H(100, 100, 50, rho_arr, depth_arr, 2)*1000000;
        printf("%f, %f", creal(res), cimag(res));

        free(rho_arr);
        free(depth_arr);
        */
        /*double *res, *freq_arr, *rho_arr;
        int i;
        rho_arr = (double *)malloc(2*sizeof(double));
        freq_arr = (double *)malloc(3*sizeof(double));
        rho_arr[0] = 1150;
        rho_arr[1] = 100;
        freq_arr[0] = 100;
        freq_arr[1] = 1000;
        freq_arr[2] = 10000;
        
        res = forward_fun_fixed_net(freq_arr, 3, 100, 50, rho_arr, 
                                    2, 4, 1.10851);
        for(i = 0; i < 6; i++)
                printf("%f\n", res[i]);
        free(rho_arr);
        free(freq_arr);
        free(res);
        */
}

