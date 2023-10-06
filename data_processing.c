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
        int freq_num, lay_num, i, j, k, pos;
        int resp_start, resp_fin, hd_pos, altt_pos, altb_pos;	
        size_t n = 0;
        long double *values;


        if(!fin || !fout){
                printf("Error opening files");
                return -1;
        }

        freq_num = (int)(args[0] + 0.1);
        lay_num = (int)(args[1] + 0.1);
        first_thick = args[2];
        step = args[3];
	resp_start = (int)(args[4] + 0.1);
	resp_fin = (int)(args[5] + 0.1);
	hd_pos = (int)(args[6] + 0.1);
	altt_pos = (int)(args[7] + 0.1);
	altb_pos = (int)(args[8] + 0.1);

        freq_arr = (long double *)malloc(freq_num*sizeof(long double));
        lower_bounds = (long double *)malloc(lay_num*sizeof(long double));
        upper_bounds = (long double *)malloc(lay_num*sizeof(long double));
        observed = (long double *)malloc((6 + freq_num)*sizeof(long double));
        x0 = (long double *)malloc(lay_num*sizeof(long double));
        P0 = (long double **)malloc(lay_num*sizeof(long double *));
        Q = (long double **)malloc(lay_num*sizeof(long double *));
        R = (long double **)malloc(2*freq_num*sizeof(long double *));
        data = (long double *)malloc(2*freq_num*sizeof(long double));
        values = (long double *)malloc((3 + 2*freq_num + 100)*sizeof(long double));
        
        for(i = 0; i < freq_num; i++)
                freq_arr[i] = args[9 + i];
        for(i = 0; i < lay_num; i++){
                lower_bounds[i] = log(0.1);
                upper_bounds[i] = log(5000.);
        }

        for(i = 0; i < lay_num; i++){
                x0[i] = log(args[3*freq_num + 9]);
                P0[i] = (long double *)malloc(lay_num*sizeof(long double));
                Q[i] = (long double *)malloc(lay_num*sizeof(long double));
                for(j = 0; j < lay_num; j++){
                        if(i==j){
                               P0[i][j] = args[3*freq_num + 10];
                        } else{
                               P0[i][j] = 0.;
                        }
                        Q[i][j] = 0.;
                }
                if(i < lay_num - 1)
                        P0[i][i+1] = args[3*freq_num + 11];
                if(i > 0)
                        P0[i][i-1] = args[3*freq_num + 11];
        }

         for(i = 0; i < 2*freq_num; i++){
                R[i] = (long double *)malloc(2*freq_num*sizeof(long double));
                for(j = 0; j < 2*freq_num; j++){
                        if(i==j){
                               R[i][j] = args[9+freq_num+i];
                        } else{
                               R[i][j] = 0.;
                        }
                }
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
                        data[j] = 0;
                for(k = 0; k < 4; k++){
                        fgets(buf,2000,fin);
                        char *p = buf;
                        n = 0;
                        for(pos = 0; n < 3 + 2*freq_num + 100&&
                        sscanf(p, "%Lf%n", values + n, &pos)==1; p+=pos){
                                ++n;
                        }
                        hor_dist+=values[hd_pos];
                        alt_T+=values[altt_pos];
                        alt_B+=values[altb_pos];
                        for(j = resp_start; j < resp_fin + 1; j++){
                                data[j-resp_start]+=values[j];
                        }
                }
                hor_dist/=4;
                alt_T/=4;
                alt_B/=4;
                for(j = 0; j < 2*freq_num; j++)
                        data[j]/=4;

                observed[0] = freq_num;
                observed[1] = lay_num;
                observed[2] = alt_T + alt_B;
                observed[3] = hor_dist;
                observed[4] = first_thick;
                observed[5] = step;
                for(j = 0; j < freq_num; j++)
                        observed[6+j] = args[9+j];

                result = kalman_step(data, observed, x0, 2*freq_num,
                             6+freq_num,
                             lay_num, forward_fun_wrapper,
                             P0, Q, R, lower_bounds, upper_bounds, 10, 1.);

                for(j = 0; j < lay_num; j++){
                        fprintf(fout, "%lf, ", exp(result[j]));
                        x0[j] = result[j];
                }
                for(j = lay_num; j < 2*lay_num+1; j++)
                        fprintf(fout, "%Lf, ", result[j]);
                
                fprintf(fout, "\n");
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
