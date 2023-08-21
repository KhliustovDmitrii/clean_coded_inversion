#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "kalman.h"
#include "forward_problem.h"

int data_inversion(char *input_file, char *output_file, long double *args)
{
/* This function reads data from file and wrights
 * inversion results
 * Params:
 * input_file, output_file: names of files
 * args: parameters for Kalman filter
 */
        FILE *fin = fopen(input_file, "rt");
        FILE *fout = fopen(output_file, "wt");
        long double *data, *observed, *x0, *result;
        long double *freq_arr, *lower_bounds, *upper_bounds;
        long double **P0, **Q, **R;
        char buf[2000];
        long double first_thick, step, temp, hor_dist, alt_T, alt_B;
        int freq_num, lay_num, i, j, k; 

        if(!fin || !fout){
                printf("Error opening files");
                return -1;
        }

        freq_num = (int)(args[0] + 0.1);
        lay_num = (int)(args[1] + 0.1);
        first_thick = args[2];
        step = args[3];

        freq_arr = (long double *)malloc(freq_num*sizeof(long double));
        lower_bounds = (long double *)malloc(lay_num*sizeof(long double));
        upper_bounds = (long double *)malloc(lay_num*sizeof(long double));
        observed = (long double *)malloc((6 + freq_num)*sizeof(long double));
        x0 = (long double *)malloc(lay_num*sizeof(long double));
        P0 = (long double **)malloc(lay_num*sizeof(long double *));
        Q = (long double **)malloc(lay_num*sizeof(long double *));
        R = (long double **)malloc(2*freq_num*sizeof(long double *));

        for(i = 0; i < lay_num; i++){
                x0[i] = log(300);
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
                               R[i][j] = 0.5;
                        } else{
                               R[i][j] = 0.;
                        }
                }
        }

        for(i = 0; i < freq_num; i++)
                freq_arr[i] = args[4 + i];
        for(i = 0; i < lay_num; i++){
                lower_bounds[i] = log(0.1);
                upper_bounds[i] = log(5000.);
        }
        

        memset(buf,0,sizeof(buf));

        for(i=0;fgets(buf,2000,fin);i++) {

                if(buf[0]=='/') {  // reading comments
                        memset(buf,0,sizeof(buf));
                        continue;
                }

                if(buf[0]=='L' || buf[0]=='B' || buf[0]=='T') { // reading Lines
                        fputs(buf,fout);
                        printf("%s",buf);
                        memset(buf,0,sizeof(buf));
                        continue;
                }

                if(strstr(buf,"*")) { // all the nonmeasurements are skipped
                        memset(buf,0,sizeof(buf));
                        continue;
                } 
                hor_dist = 0;
                alt_T = 0;
                alt_B = 0;
                for(j = 0; j < 2*freq_num; j++)
                        data[i] = 0;
                for(k = 0; k < 4; k++){
                        fgets(buf,2000,fin);
                        sscanf(buf, "%Lf", &temp);
                        hor_dist+=temp;
                        sscanf(buf, "%Lf", &temp);
                        alt_T+=temp;
                        sscanf(buf, "%Lf", &temp);
                        alt_B+=temp;
                        for(j = 0; j < 2*freq_num; j++){
                                sscanf(buf, "%Lf", &temp);
                                data[i]+=temp;
                        }
                }
                hor_dist/=4;
                alt_T/=4;
                alt_B/=4;
                for(j = 0; j < 2*freq_num; j++)
                        data[i]/=4;

                observed[0] = freq_num;
                observed[1] = lay_num;
                observed[2] = alt_T + alt_B;
                observed[3] = hor_dist;
                observed[4] = first_thick;
                observed[5] = step;
                for(j = 0; j < freq_num; j++)
                        observed[6+j] = args[4+j];

                result = kalman_sequential(data, observed, x0, 2*freq_num,
                             6+freq_num,
                             lay_num, forward_fun_wrapper,
                             P0, Q, R, lower_bounds, upper_bounds, 10, 1.);

                for(j = 0; j < freq_num; j++)
                        fprintf(fout, "%lf ", exp(result[i]));

                for(j = freq_num; j < 2*freq_num+1; j++)
                        fprintf(fout, "%Lf ", result[i]);
                }
        free(lower_bounds);
        free(upper_bounds);
        free(P0);
        free(Q);
        free(R);
        free(x0);
        free(freq_arr);
        free(result);
        free(data);
        free(observed);
        fclose(fin);
        fclose(fout);
        return 0;
}
