#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include "math_utils.h"
#include "forward_problem.h"
#include "kalman.h"

int main()
{
        double *depth_arr, *rho_arr;
        double complex imp;
        int lay_num;

        rho_arr = (double *)malloc(3*sizeof(double));
        depth_arr = (double *)malloc(2*sizeof(double));

        rho_arr[0] = 100;
        rho_arr[1] = 200;
        rho_arr[2] = 300;

        depth_arr[0] = 10;
        depth_arr[1] = 20;

        imp = H(100, 100, 50, rho_arr, depth_arr, 3)*1.E7;

        printf("%.10f + %.10fi", creal(imp), cimag(imp));
        free(depth_arr);
        free(rho_arr);
        return 0;
}
