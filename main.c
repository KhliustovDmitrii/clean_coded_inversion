#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include "math_utils.h"
#include "forward_problem.h"
#include "kalman.h"
#include "data_processing.h"

int main()
{       

        FILE *conf = fopen("configuration.txt", "rt");
        char buf[2000];
        int i, freq_num, pos;
        size_t n = 0;
        long double *values, *args;

        memset(buf,0,sizeof(buf));

        //Read freq num
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        freq_num = atoi(buf);
        values = (long double *)malloc((2*freq_num)*sizeof(long double));
        args = (long double *)malloc((7 + 3*freq_num)*sizeof(long double));
        args[0] = freq_num;
        
        //Read lay_num, first_thick, step
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        char *p = buf;
        for(pos = 0; n < 3 &&
                sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                        ++n;
                }
        n = 0;
        args[1] = values[0];
        args[2] = values[1];
        args[3] = values[2];

        //Read freqs
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        p = buf;
        for(pos = 0; n < freq_num &&
                sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                        ++n;
                }
        for(i = 0; i < freq_num; i++)
                args[4+i] = values[i];
        n = 0;
        //Read R_i_i
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        p = buf;
        for(pos = 0; n < 2*freq_num &&
                sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                        ++n;
                }
        for(i = 0; i < 2*freq_num; i++)
                args[4+freq_num+i] = values[i];
        n = 0;
        //Read initial half-space resistivity
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        p = buf;
        for(pos = 0; n < freq_num &&
                sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                        ++n;
                }
        args[3*freq_num + 4] = values[0];
        n = 0;
        //Read P0[0][0], P0[0][1]
        fgets(buf, 2000, conf);
        fgets(buf, 2000, conf);
        p = buf;
        for(pos = 0; n < 2 &&
                sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                        ++n;
                }

        args[3*freq_num + 5] = values[0];
        args[3*freq_num + 6] = values[1];
        n = 0;
        i = data_inversion("data.txt", "inversion_result.txt", args);
        free(values);
        free(args);
        return 0;
        


        /* 
        long double *rho_arr, *freq_arr, *resp, *res, *observed, *x0;
        long double *lower_bounds, *upper_bounds;
        long double **P0, **Q, **R;
        long double complex imp;
        int i, j, lay_num, freq_num;
        lay_num = 2;
        freq_num = 5;

        rho_arr = (long double *)malloc(lay_num*sizeof(long double));
        freq_arr = (long double *)malloc(freq_num*sizeof(long double));
        lower_bounds = (long double *)malloc(lay_num*sizeof(long double));
        upper_bounds = (long double *)malloc(lay_num*sizeof(long double));
        observed = (long double *)malloc((6 + freq_num)*sizeof(long double));
        x0 = (long double *)malloc(lay_num*sizeof(long double));
        P0 = (long double **)malloc(lay_num*sizeof(long double *));
        Q = (long double **)malloc(lay_num*sizeof(long double *));
        R = (long double **)malloc(2*freq_num*sizeof(long double *));

        for(i = 0; i < lay_num; i++){
                rho_arr[i] = 100.+50*i;
                x0[i] = log(500.);
                P0[i] = (long double *)malloc(lay_num*sizeof(long double));
                Q[i] = (long double *)malloc(lay_num*sizeof(long double));
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
                R[i] = (long double *)malloc(2*freq_num*sizeof(long double));
                for(j = 0; j < 2*freq_num; j++){
                        if(i==j){
                               R[i][j] = 1;
                        } else{
                               R[i][j] = 0.;
                        }
                }
        }

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = 100. + 300*i;
        for(i = 0; i < lay_num; i++){
                lower_bounds[i] = log(0.1);
                upper_bounds[i] = log(3000.);
        }
        //imp = H(100, 100, 50, rho_arr, depth_arr, 3)*1.E7;
        //imp = u(0.01, 100, 100, rho_arr, depth_arr, 3)*1.E7;
        //imp = spec_dens(0.001, pars);
        resp = forward_fun_fixed_net(freq_arr, freq_num, 100., 50.,
                                    rho_arr, lay_num, 4, 1.1);
        
        for(i = 0; i < 2*freq_num; i++)
                printf("%Lf\n", resp[i]);

        observed[0] = freq_num;
        observed[1] = lay_num;
        observed[2] = 100.;  
        observed[3] = 50.;
        observed[4] = 4;
        observed[5] = 1.1;
        for(i = 0; i < freq_num; i++)
                observed[6+i] = freq_arr[i];

        res = kalman_sequential(resp, observed, x0, 2*freq_num, 6+freq_num,
                          lay_num, forward_fun_wrapper,
                          P0, Q, R, lower_bounds, upper_bounds, 10, 1.);

        for(i = 0; i < lay_num; i++)
                printf("result: %lf\n", exp(res[i]));
        for(i = lay_num; i < 2*lay_num+1; i++)
                printf("result: %Lf\n", res[i]);


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
        /*long double *rho_arr, *depth_arr;
        long double complex res;

        rho_arr = (long double *)malloc(2*sizeof(long double));
        depth_arr = (long double *)malloc(1*sizeof(long double));
        rho_arr[0] = 1150;
        rho_arr[1] = 100;
        depth_arr[0] = 10;

        res = H(100, 100, 50, rho_arr, depth_arr, 2)*1000000;
        printf("%f, %f", creal(res), cimag(res));

        free(rho_arr);
        free(depth_arr);
        */
        /*long double *res, *freq_arr, *rho_arr;
        int i;
        rho_arr = (long double *)malloc(2*sizeof(long double));
        freq_arr = (long double *)malloc(3*sizeof(long double));
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

