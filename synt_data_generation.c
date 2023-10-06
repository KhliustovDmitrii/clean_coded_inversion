#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include "math_utils.h"
#include "forward_problem.h"

int main()
{
        int i, j, freq_num, lay_num;
        long double *freqs, *vals, *rhos;
	FILE *fout = fopen("synt_data.txt", "wt");

        freq_num = 15;
	lay_num = 25;
        freqs = (long double *)malloc(freq_num*sizeof(long double));
        vals = (long double *)malloc(2*freq_num*sizeof(long double));
        rhos = (long double *)malloc(lay_num*sizeof(long double));

        for(i = 0; i < freq_num; i++)
                freqs[i] = 75*(i+1)*(i+1);

        for(j = 0; j < 10; j++){
                for(i = 0; i < lay_num; i++)
		        rhos[i] = (i-j)*(i-j)*2 + 100;
		vals = forward_fun_fixed_net(freqs, freq_num, 100, 50, rhos, 
				lay_num, 4, 1.1085);
		fprintf(fout, "50 50 50 ");
	        for(i = 0; i < 2*freq_num; i++)
			fprintf(fout, "%Lf ", vals[i]);
		fprintf(fout, "\n");
        }
	fclose(fout);
        free(freqs);
	free(vals);
	free(rhos);
}
