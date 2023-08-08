#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "math_utils.h"
#include "forward_problem.h"
#include "kalman.h"

int main()
{
        double *rho_arr, *freq_arr, *res;
        double complex imp;
        int i;

        rho_arr = (double *)malloc(10*sizeof(double));
        freq_arr = (double *)malloc(1000*sizeof(double));

        for(i = 0; i < 10; i++)
                rho_arr[i] = 100+50*i;
        for(i = 0; i < 1000; i++)
                freq_arr[i] = 100 + 100*i;

        //imp = H(100, 100, 50, rho_arr, depth_arr, 3)*1.E7;
        //imp = u(0.01, 100, 100, rho_arr, depth_arr, 3)*1.E7;
        //imp = spec_dens(0.001, pars);
        res = forward_fun_fixed_net(freq_arr, 1000, 100, 50,
                                    rho_arr, 10, 5, 1.1);
        for(i = 0; i < 500; i++)
                printf("%f, ", res[i]);
        printf("\n\n\n");
        for(i = 500; i < 1000; i++)
                printf("%f, ", res[i]);
        free(freq_arr);
        free(rho_arr);
        free(res);
        return 0;
}
