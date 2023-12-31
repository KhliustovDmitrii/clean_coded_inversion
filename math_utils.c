#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "math_utils.h"

/* 
 *This file contains useful mathematical functions, which do not depend
 *on actual geological parameters
 */

long double bessel_function(long double x)
{
/*
 * This function computes Bessel function of first kind order zero
 *Params:
 *        x: long double, point in which to evaluate
 *Returns:
 *        res: long double, the value of the function
 */
        long double abs_x, x_squared, x_corrected, z, z_squared, term_1, term_2;
        long double res;

        abs_x = fabs(x);
        if(abs_x < 8.0){
                x_squared = x*x;
                term_1 = 57568490574.0 + x_squared*(-13362590354.0 + x_squared*
                         (651619640.7 + x_squared*(-11214424.18 + x_squared*(
                         77392.33017 + x_squared*(-184.9052456)))));
                term_2 = 57568490411.0 + x_squared*(1029532985.0 + x_squared*(
                         9494680.718 + x_squared*(59272.64853 + x_squared*(
                         267.8532712 + x_squared))));
                res = term_1/term_2;
        } else {
                z = 8.0/abs_x;
                z_squared = z*z;
                x_corrected = abs_x - 0.785398164;
                term_1 = 1.0 + z_squared*(-0.1098628627e-2 + z_squared*(
                         0.2734510407e-4 + z_squared*(-0.2073370639e-5 + 
                         z_squared * 0.2093887211e-6)));
                term_2 = -0.1562499995e-1 + z_squared*(0.1430488765e-3 + 
                         z_squared*(-0.6911147651e-5 + z_squared*(
                         0.7621095161e-6 - z_squared*0.934935152e-7)));
                res = sqrt(0.636619772/abs_x)*(cos(x_corrected)*
                      term_1 - z*sin(x_corrected)*term_2);
        }
        return res;
}

long double scalar_product(long double *u, long double *v, int vec_len)
{
/*
 * This function computes scalar product of two vectors
 *Params:
 *        u, v: array of long double, vectors
 *        vec_len: int, dimension of vectors
 *Returns:
 *        res: long double, scalar product
 */
        long double res;
        int i; 

        res = 0;
        for(i = 0; i < vec_len; i++)
                res = res + u[i]*v[i];
        return res;
}

long double **matrix_product(long double **M, long double **N, int m1, int m2, int n1, int n2)
{
/*
 * This function computes product of two matrices
 *Params:
 *        M, N: two dimensional array of long double, matrices
 *        m1, m2: int, number of rows and columns for M
 *        n1, n2: int, number of rows and columns for N
 *Returns:
 *        res: two dimensional array of long double, matrix product
 */

        long double **res;
	int i, j, k;

        res = (long double **)malloc(m1*sizeof(long double *));
        for(i = 0; i < m1; i++)
                res[i] = (long double *)malloc(n2*sizeof(long double));

        if(m2!=n1){
                printf("Matrix dimensions incompatible, %d != %d", m2, n1);
                return res;
        }

        for(i = 0; i < m1; i++){
                for(j = 0; j < n2; j++){
                        res[i][j] = 0;
                        for(k = 0; k < m2; k++)
                                res[i][j] = res[i][j] + M[i][k]*N[k][j];
                }
        }
        return res;
}

long double **pseudo_inv(long double **M, int m1, int m2)
{
/*
 * This function computes pseudoinverse of a matrix
 *Params:
 *        M: matrix
 *        m1, m2: row and column numbers of M
 */
        long double **res, **MT;
        MT = transpose(M, m1, m2);
        res = matrix_product(MT, M, m2, m1, m1, m2);
        res = invert(res, m2);
        res = matrix_product(res, MT, m2, m2, m2, m1);
        free(MT);
        return(res);
}

long double **invert(long double **mat, int n)
{
/*
 * This function inverts matrix
 *Params:
 *        M: matrix
 *        m: dimensionality
 */
        int i, j, k, ind_row;
        long double temp;
        long double **res, **M;

        res = (long double **)malloc(n*sizeof(long double *));
        M = (long double **)malloc(n*sizeof(long double *));
        for(i = 0; i < n; i++){
                res[i] = (long double *)malloc(n*sizeof(long double));
                M[i] = (long double *)malloc(n*sizeof(long double));
                for(j = 0; j < n; j++){
                        res[i][j] = 0;
                        M[i][j] = mat[i][j];
                }
        }
        for(i = 0; i < n; i++)
                res[i][i] = 1;

        for(i = 0; i < n; i++){
                ind_row = i;
                while(fabs(M[ind_row][i])<0.0001&&ind_row<n)
                        ind_row++;
                if(ind_row > i){
                        for(j = 0; j < n; j++){
                                res[i][j]+=res[ind_row][j];
                                res[ind_row][j] = res[i][j] - res[ind_row][j];
                                res[i][j]-=res[ind_row][j];
                                M[i][j]+=M[ind_row][j];
                                M[ind_row][j] = M[i][j] - M[ind_row][j];
                                M[i][j]-=M[ind_row][j];
                        }
                }
                temp = 1/M[i][i];
                for(j = 0; j < n; j++){
                        M[i][j]*=temp;
                        res[i][j]*=temp;
                }
                for(j = i+1; j < n; j++){
                        temp = M[j][i];
                        for(k = i; k < n; k++){
                                res[j][k] = res[j][k] - temp*res[i][k];
                                M[j][k] = M[j][k] - temp*M[i][k];
                        }
                }
        }
        for(i = n-1; i >= 0; i--){
                for(j = i-1; j >= 0; j--){
                        temp = -M[j][i]/M[i][i];
                        for(k = i; k > j-1; k--){
                                M[j][k] = M[j][k] + temp*M[i][k];
                                res[j][k] = res[j][k] + temp*res[i][k];
                        }
                }
        }
        for(i = 0; i <n; i++){
                temp = 1/M[i][i];
                for(j = 0; j < n; j++){
                         res[i][j]*=temp;
                }
        }

return res;
} 

long double **transpose(long double **M, int m1, int m2)
{
/*
 * This function transposes matrix
 *Params:
 *        M: matrix
 *        m1, m2: int, number of rows and columns for M
 */

        long double **res;
        int i, j;

        res = (long double **)malloc(m2*sizeof(long double *));
        for(i = 0; i < m2; i++){
                res[i] = (long double *)malloc(m1*sizeof(long double));
                for(j = 0; j < m1; j++)
                        res[i][j] = M[j][i];
        }
        
        return res;
}

long double **cholesky(long double **M, int m)
{
/*
 * This function computes Cholesky decomposition of matrix M
 *Params:
 *        M: matrix
 *        m: dimensionality of M
 */

        long double **res;
        long double temp;
        int i, j, k;
        
        res = (long double **)malloc(m*sizeof(long double *));
        for(i = 0; i < m; i++){
                res[i] = (long double *)malloc(m*sizeof(long double));
                for(j = 0; j < m; j++)
                        res[i][j] = 0;
        }
        for(i = 0; i < m; i++){
                temp = 0;
                for(j = 0; j < i; j++) 
                        temp+=res[i][j]*res[i][j];
                res[i][i] = sqrt(M[i][i] - temp);
                for(j = i+1; j < m; j++){
                        temp = 0;
                        for(k = 0; k < i; k++){
                               temp+=res[j][k]*res[i][k];
                        }
                        res[j][i] = 1/res[i][i]*(M[j][i] - temp);
                }
        }
        return res;
}



long double primal_field(long double hor_dist, long double ver_dist)
{
/*
 * This function computes primal field amplitude induced by emitter in receiver
 *Params:
 *        hor_dist: long double, horizontal distance between E and R
 *        ver_dist: long double, vertical distance between E and R
 *Returns:
 *        res: long double, primal field amplitude
 */

        long double k, M, MR, res;
        long double Hp[3];
        long double R[3] = {hor_dist, 0, ver_dist};
       
        MR = scalar_product(R, R, 3);
        k = 1 / MR / sqrt(MR) / 4 / M_PI;
        Hp[0] = (3*R[0]*R[2] / MR)*k;
        Hp[1] = (3*R[1]*R[2] / MR)*k;
        Hp[2] = (3*R[2]*R[2] / MR - 1)*k;
        
        res = Hp[0]*Hp[0] + Hp[1]*Hp[1] + Hp[2]*Hp[2];

        res = sqrt(res);
        return res;
}

long double complex integral(long double complex (*f)(long double, long double*), long double *f_pars, 
                        int par_num, long double t0, long double t1, long double delta)
{
/*
 * This function integrates function f from t0 to t1 with step delta
 *Params:
 *        *f: pointer to function
 *        f_pars: array of long double, fixed parameters of f
 *        t0, t1: long double, integration limits
 *        delta: long double, discretization parameter
 */

        long double complex result;
        long double step;

        result = 0;
        step = t0;

        while (step<t1){
                result = result + f(step, f_pars)*delta;
                step+=delta;
        } 

        return result;
}

long double **jacobian(long double* (*f)(long double*, long double*), long double *f_pars, long double *point,
                  int in_dim, int out_dim, int par_dim, long double delta)
{
/*
 * This function computes jacobian of f(x, par) in point x w.r.t. dx
 *Params:
 *        *f: pointer to function
 *        f_pars: array of long double, fixed parameters of f
 *        *_dim: int, dimension of x, f(x), and parameters
 *        point: array of long double, x
 *        delta: long double, discretization parameter
 */
        int i, j;
        long double *val, *val_shifted;
        long double **jac;

        jac = (long double **)malloc(out_dim*sizeof(long double *));
        for(i = 0; i < out_dim; i++)
                jac[i] = (long double *)malloc(in_dim*sizeof(long double));
        
        val = f(point, f_pars);
        for(i = 0; i < in_dim; i++){
                point[i]+=delta;
                val_shifted = f(point, f_pars);
                for(j = 0; j < out_dim; j++)
                        jac[j][i] = (val_shifted[j] - val[j])/delta;
                point[i]-=delta;        
        }
        return jac;
}

long double discrepancy(long double *diff, long double **R, int dim)
{
/*
 * This function computes weighted MSE used as discrepancy measure
 *Params:
 *        diff: vector observed - predicted
 *        R: covariance matrix
 *        dim: dimensionality of vector
 */
        long double res;
        int i;
        res = 0;
    
        for(i = 0; i < dim; i++)
                res+=(diff[i]*diff[i])/R[i][i];

        printf("diff = %Lf\n", res); 
        return sqrt(res);
}
