#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define mu0 (M_PI*.4e-6)
#define DELTA   .01

// Number of frequencies
// this is for modelled data:
//#define NFREQS 26

// this is for Croatia (Old Equator w/o 50Hz harm: 850 & 2700 Hz)
//#define NFREQS 21

//this is for Maar-Siene Yakutia
#define NFREQS 15
//#define NFREQS 23
//#define NFREQS 17
//#define NFREQS 10
//#define NFREQS 22
// Number of layers
//#define NLAYERS 10
//#define NPLAYERS 4
//#define DSTEP 5

//#define NLAYERS 1
// # of free thickness upper layers:
#define NFLAYERS 0
// # of fixed thickness layers under the free ones:
#define NLAYERS 25
// # of polarized layers:
#define NPLAYERS 0
// 4 meters is the best value for the salt water (according to Croatian data)

#define FIRST_THK 4.0
// The same as for Aarhus Geoph
#define DSTEP 1.10851
//#define DSTEP 1.0
#define MIN_RES 0.02
#define MAX_RES 20000.

//*********** Parameters to tune inversion *************//
// Channel number to calculate app. res. for initial halfspace
#define BASE_CH 2
// Initial halfspace resistivity
#define RES_INI 150
#define DEP_INI 25
// Averaging interval
#define AVERAGE 4
// Maximum iterations number
#define MAX_ITER 10
// Initial Error - RMS for ln(Resistivity)
#define ERR_INI 0.3
// Initial Cross Correlation - should not be greater then ERR_INI!!!
#define COR_INI 0.1
// Delta Altitude
#define DA -4
//#define DA 0
#define STOP_VAL 1.0

//#define CHARGEABILITY

typedef struct {
    double coord;
    double coord_y;
    double hor_dist;
    double ver_dist;
    double relief;
    double alt;
    double prim;
    double w[2*NFREQS];
} geometry;

double freqs[NFREQS] =
//0  1    2   3  4   5    6    7   8    9    10   11   12   13   14   15   16   17  18   19   20   21   22(3) 23(3) 24(3)   25(4)
//{77,231,385,540,694,848,1003,1466,1620,1774,1929,2083,2391,2546,2700,2855,3009,3163,3626,5324,6558,7484, 8101, 9645, 12679, 14004};
//{77,231,385,540,694,848,1003,1157,1311,1466,1620,1774,2855,3627,4243,4398,5169,5324,5478,9645,9799,14737,14891,15046,15200,15354};
//{77,231,386,540,694,849,1003,1466,1620,2083,2855,3009,3164,3781,4398,5324,5478,6096,6559,6867,7022,7639,8102,9645,9799,9954,10880,11034,11188,12423,12731,13194,14583,14738,14892,15046,15355};
//{77,231,386,540,694,1003,1466,1620,2083,3009,3781,4398,5324,6096,6559,6867,7022,7639,8102,9799,11034,12731,15046};
//{77,231,386,540,694,1003,1466,1620,2083,3009,3781,4398,5324,6559,8102,9799,15046};
//{77,231,386,3164,5324,6713,8102,9799,11188,14583};
//{77,231,385,540,694,1003,1157,1466,1620,1775,2083,2238,2701,3009,3318,3472,3627,3781,4244,6250,9799,12269};
//{77,231.5,386,540,694.5,1003,1157,1312,1466,1620,1775,1929,2083,2238,2855,3009,3164,3318,3472,6173,13503}; // from Croatia, 21 frqs
{77,231.5,386,540,694.5,848.77,1003,1466,1620,1775,2855,4244,5324,9799,15046}; // from Maar-Siene, 15 frqs
//{77,231,385,540,694,848,1003,1466,1620,1774,1929,2083,2391,2546,2700,2854,3009,3163,3626,5324,6558,7484,8101,9645,12731,14120}; // for modelling!!!

double upper[2*NPLAYERS+NLAYERS+2*NFLAYERS];
double lower[2*NPLAYERS+NLAYERS+2*NFLAYERS];

double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

double complex PartSum(double n0,double hh,double complex n1,double r,double complex Imp)
{
    double complex s;
    if(fabs(r)>.001)
        s = bessj0(n0*r);
    else {
        s  = 1;
    }
    double complex A = exp(-n0*hh)*(n1-n0*Imp)*n0*n0*.25/(n1+n0*Imp)/M_PI;
    s = A*s;
    return s;
}

double complex Impedance(int n, double n0, double om, double *rho, double *dep) {
    int i;
    double complex ni,nim1;
    double complex Imp;
    Imp = 1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[n-1]);
    for(i=n-1;i>0;i--) {
    ni = nim1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[i-1]);
    Imp = ctanh(nim1*dep[i-1]+catanh(nim1/ni*Imp));
    }
    return Imp;
}

double complex CImpedance(int n, double n0, double om, double complex *rho, double *dep) {
    int i;
    double complex ni,nim1;
    double complex Imp;
    Imp = 1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[n-1]);
    for(i=n-1;i>0;i--) {
    ni = nim1;
    nim1 = csqrt(n0*n0 - I*om*mu0/rho[i-1]);
    Imp = ctanh(nim1*dep[i-1]+catanh(nim1/ni*Imp));
    }
    return Imp;
}

double complex integral(int n,double hh,double r,double *rho,double *dep,double f)
{
    double complex PS;         // f  -  Hz!!!
    double complex intl = 0;
    double dn0;
    double n0=0;//in0;
    double complex n1,c;
    double sigma = 1./rho[0];
    double complex Imp;
    double om = f*2*M_PI;
    c = I*(-om*sigma*mu0);        //  4pi * 1.e-7 * 2pi (Hz -> 1/sec)

    #define VAL .001
    for(n0=VAL,dn0=VAL;n0<1;n0+=dn0) {
        n1 = csqrt(n0*n0+c);
        Imp = Impedance(n,n0,om,rho,dep);
        PS  = PartSum(n0,hh,n1,r,Imp);
        intl += dn0*PS;
	}
    return intl;
}

double complex Cintegral(int n,double hh,double r,double complex *rho,double *dep,double f)
{
    double complex PS;         // f  -  Hz!!!
    double complex intl = 0;
    double dn0;
    double n0=0;//in0;
    double complex n1,c;
    double complex sigma = 1./rho[0];
    double complex Imp;
    double om = f*2*M_PI;
    c = I*(-om*sigma*mu0);        //  4pi * 1.e-7 * 2pi (Hz -> 1/sec)

    #define VAL .001
    for(n0=VAL,dn0=VAL;n0<1;n0+=dn0) {
        n1 = csqrt(n0*n0+c);
        Imp = CImpedance(n,n0,om,rho,dep);
        PS  = PartSum(n0,hh,n1,r,Imp);
        intl += dn0*PS;
    }
    return intl;
}

double complex ImHz(int n, double r,double z,double f,double *rho, double *dep) {
    return integral(n,z,r,rho,dep,f);
}

double complex CImHz(int n, double r,double z,double f,double complex *rho, double *dep) {
    return Cintegral(n,z,r,rho,dep,f);
}

void fdfun(geometry geo,int nlay,int bfr,double *x,double *y) {
    int i;
    double complex refl;
    for(i=0;i<NFREQS;i++) {
        refl = ImHz(nlay,geo.hor_dist,2*geo.alt+geo.ver_dist,freqs[i],x,&(x[nlay]))*I/geo.prim ;
        y[2*i] = creal(refl);
        y[2*i+1] = cimag(refl);
        if(i>bfr) break;
    }
    if(1) {//bfr==NFREQS) {
        for(i=0;i<bfr;i++) {
            if(i<NFREQS-1)
                y[2*i+1] = y[2*(i+1)+1] - y[2*i+1];
            else
                y[2*i+1] = y[2*(i-1)+1];
        }
    }
}

void Cfdfun(geometry geo,int nlay,int bfr,double *dep,double *x,double *y) {
    int i,j;
    double complex refl;
    double complex xx[NLAYERS+NFLAYERS];
    for(i=0;i<NFREQS;i++) {
//        for(j=0;j<NLAYERS+NFLAYERS;j++) xx[j] = x[j]+I*x[NLAYERS+NFLAYERS+j]/freqs[i]/2/M_PI;
        for(j=0;j<NLAYERS+NFLAYERS;j++)
            // version 1: RealFrequencyIndependent and ComplexFrequencyDependent resistivities in series (Mine)
            xx[j] = x[j]+((j<NPLAYERS)?(x[NLAYERS+NFLAYERS+j]/(1-I*freqs[i]*2*M_PI*x[NLAYERS+NFLAYERS+NPLAYERS+j])):(0));
            //-------------------------------------------------------------------------------------------
            // version 2: Cole-Cole from Andrea's (inverted sign for i*w) - resistivities in
            // xx[j] = x[j]-((j<NPLAYERS)?(x[NLAYERS+NFLAYERS+j]*(1.-1./(1-I*freqs[i]*2*M_PI*x[NLAYERS+NFLAYERS+NPLAYERS+j]))):(0));
        refl = CImHz(nlay,geo.hor_dist,2*geo.alt+geo.ver_dist,freqs[i],xx,dep)*I/geo.prim ;
        y[2*i] = creal(refl);
        y[2*i+1] = cimag(refl);
    }
}

// FK 1 step
void proc(int n,double *x,double *S,double z,double *h,double sg2) {// Square Root Matrix method for LSM(Kalman)
    int i,j,k;
    double f[2*NPLAYERS+NLAYERS+2*NFLAYERS],e[2*NPLAYERS+NLAYERS+2*NFLAYERS],d[2],bk,ck,dz;
    d[0] = sg2;
    for(i=0;i<n;i++) {
        f[i] = 0;
        for(j=0;j<n;j++)
            f[i]+=S[j*n+i]*h[j]; // S^T !!!
    }
    memset(e,0,sizeof(e));
    for(k=0;k<n;k++) {
        d[1] = d[0] + f[k]*f[k];
        bk = sqrt(d[0]/d[1]);
        ck = f[k]/sqrt(d[0]*d[1]);
        for(i=0;i<n;i++) {
            double tmp = S[i*n+k]*f[k];
            S[i*n+k] = bk*S[i*n+k]-ck*e[i];
            e[i] += tmp;
        }
        d[0] = d[1];
    }
    dz = z;
    for(i=0;i<n;i++) dz -= h[i]*x[i];
    dz/=d[0];
    for(i=0;i<n;i++) x[i] += e[i]*dz;
}

// combined free and fixed layers inversion
// nlay shoud be the total number of layers
int flinversion(geometry geo,
                int bfr,// now it is the # of used frqss
                int nlay,
                double *x_ini,
                double *dpth,
                double *y_ini,
                double *y_mes,
                double *residual,
                int *up,
                double *S) {// both x_ and y_ are in log-axes

    double y1[2*NFREQS];
    double dx[2*NPLAYERS+2*NLAYERS+2*NFLAYERS],x0[2*NPLAYERS+2*NLAYERS+2*NFLAYERS],x1[2*NPLAYERS+2*NLAYERS+2*NFLAYERS],xini[2*NPLAYERS+2*NLAYERS+2*NFLAYERS];
    int charge = 1;


#ifdef CHARGEABILITY
//    double Jacob[2*NFREQS][2*NLAYERS+2*NFLAYERS];
//    charge =2;
    double Jacob[2*NFREQS][2*NPLAYERS+NLAYERS+2*NFLAYERS];
    charge =2;
#else
    double Jacob[2*NFREQS][NLAYERS+2*NFLAYERS];
#endif
    double res = 0;
    int i,j;

    memset(dx,0,sizeof(dx));
    memset(Jacob, 0,sizeof(Jacob));
#ifdef CHARGEABILITY
    for(i=0;i<charge*NPLAYERS+nlay;i++) { // we can use chargeability for fixed layers only
        xini[i] = x_ini[i]; // first nlay items are resistivities (real), the rest are "impedivities" (imaginary)
    }
#else
    for(i=0;i<nlay+NFLAYERS;i++) {
        xini[i] = x_ini[i];
        if(i<nlay-1)
            xini[i+nlay+NFLAYERS] = x0[i+nlay+NFLAYERS] = x1[i+nlay+NFLAYERS] = dpth[i];
    }
#endif

    // against nans!!
    for(i=0;i<2*bfr;i++)
        if(y_mes[i]<.001) y_mes[i] = .001;

    // first forward calculation for the model
    if(*residual<0) {
#ifdef CHARGEABILITY
        Cfdfun(geo,nlay,bfr,dpth,xini,y_ini);
#else
        fdfun(geo,nlay,bfr,xini,y_ini);
#endif
        *residual = 0;
        for(j=0;j<2*bfr;j++) {
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
            *residual += val;//(val>*residual)?val:*residual;
        }
    }

    // a good enough first approach
    if(*residual<0.01)
        return 0;

    // Jacobian matrix calculation
    int val = (bfr==1)?nlay:(charge*NPLAYERS+nlay+NFLAYERS);
    for(i=0;i<val;i++) {
        for(j=0;j<val;j++) {
            double k = 10;
            x1[j] = x_ini[j] * ( (i!=j)? 1. : (1.-k*DELTA) );
        }
#ifdef CHARGEABILITY
        Cfdfun(geo,nlay,bfr,dpth,x1,y1);
#else
        fdfun(geo,nlay,bfr,x1,y1);
#endif
        double div = log(x1[i]/x_ini[i]);
        for(j=0;j<2*bfr;j++) {
            Jacob[j][i] = log(fabs(y1[j]/y_ini[j]));
            Jacob[j][i]/=div;
        }
    }

    // Step (by step) a-la FK
    for(j=0;j<2*bfr;j++) {
        proc(val,dx,S,log(fabs(y_ini[j]/y_mes[j])),Jacob[j],1./(geo.w[j]*y_mes[j]*geo.w[j]*y_mes[j]));
    }

    for(i=0;i<val;i++) {
        x0[i] = xini[i]*exp(-dx[i]);
        if(isnan(x0[i])) x0[i] = xini[i];
        if(x0[i]>upper[i]) x0[i] = upper[i];
        if(x0[i]<lower[i]) x0[i] = lower[i];
    }
#ifndef CHARGEABILITY
    for(i=nlay;i<2*nlay-1;i++)
        x0[i] = x1[i];
#endif

#ifdef CHARGEABILITY
    Cfdfun(geo,nlay,bfr,dpth,x0,y_ini);
#else
    fdfun(geo,nlay,bfr,x0,y_ini);
#endif
    for(j=0;j<2*bfr;j++){
        double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
        res += val;//(val>res)?val:res;
    }

    int cntr = 0;
    while(res>*residual*1.01) {
        if(cntr++ > 3) break;
        for(i=0;i<val;i++) {
            dx[i]*=.5;
            x0[i] = x_ini[i]*exp(-dx[i]);
            if(isnan(x0[i])) x0[i] = xini[i];
            if(x0[i]>upper[i]) x0[i] = upper[i];
            if(x0[i]<lower[i]) x0[i] = lower[i];
        }

        res = 0;
#ifdef CHARGEABILITY
        Cfdfun(geo,nlay,bfr,dpth,x0,y_ini);
#else
        fdfun(geo,nlay,bfr,x0,y_ini);
#endif
        for(j=0;j<2*bfr;j++){
            double val = (y_ini[j]-y_mes[j])*(y_ini[j]-y_mes[j])*geo.w[j]*geo.w[j]/(2*bfr);
            res += val;//(val>res)?val:res;
        }
    }

    if(res>*residual) up[0]++;
    else up[0] = 0;
    *residual = res;
    for(i=0;i<charge*NPLAYERS+nlay;i++)
        x_ini[i] = x0[i];
    if(bfr>1) {
        for(i=0;i<NFLAYERS;i++) {
            dpth[i] = x0[i+charge*NPLAYERS+nlay];
            x_ini[i+nlay] = dpth[i];
        }
    }
    return 1;
}

double primField(double hd,
                 double vd)
{

    double E[10], R[4], RR[10], Hp[4], k, MR, Ampl, M;
    R[1]=hd; R[2]=0; R[3]=vd; //towed cable = 39.11
    //R[1]=33; R[2]=0; R[3]=21; //towed cable = 39.11
    //M[1]=0; M[2]=0; M[3]=1.e9;
    M = 1;
    E[1]=1; E[2]=0; E[3]=0;
    E[4]=0; E[5]=1; E[6]=0;
    E[7]=0; E[8]=0; E[9]=1;
    RR[1]=R[1]*R[1]; RR[2]=R[1]*R[2]; RR[3]=R[1]*R[3];
    RR[4]=R[2]*R[1]; RR[5]=R[2]*R[2]; RR[6]=R[2]*R[3];
    RR[7]=R[3]*R[1]; RR[8]=R[3]*R[2]; RR[9]=R[3]*R[3];
    MR=R[1]*R[1]+R[2]*R[2]+R[3]*R[3];
    k = M / MR / (sqrt(MR)) / 4 / M_PI;
    Hp[1] = (3*RR[3] / MR-E[3]) * k;
    Hp[2] = (3*RR[6] / MR-E[6]) * k;
    Hp[3] = (3*RR[9] / MR-E[9]) * k;

    Ampl = sqrt(Hp[1] * Hp[1] + Hp[2] * Hp[2] + Hp[3] * Hp[3]);
    return Ampl;
}

int main(int argc, char **argv)
{
    if(argc != 3) {
        printf("inversion file_in file_out\n");
        return 0;
    }

    geometry geo;
    char buf[2000];
    char *data;
    //char time[23];
    //for double spaces
    char time[13];
    double y_mes[2*NFREQS];
    int charge = 1;
#ifdef CHARGEABILITY
    charge = 2;
    double rho[charge*NPLAYERS+NLAYERS+NFLAYERS];
    double x_ini[charge*NPLAYERS+NLAYERS+2*NFLAYERS];
    double S_ini[(charge*NPLAYERS+NLAYERS+2*NFLAYERS)*(charge*NPLAYERS+NLAYERS+2*NFLAYERS)];
    double S_pre[(charge*NPLAYERS+NLAYERS+2*NFLAYERS)*(charge*NPLAYERS+NLAYERS+2*NFLAYERS)];
    double S0[(charge*NPLAYERS+NLAYERS+2*NFLAYERS)*(charge*NPLAYERS+NLAYERS+2*NFLAYERS)];
#else
    double rho[NLAYERS+NFLAYERS];
    double x_ini[NLAYERS+2*NFLAYERS];
    double S_ini[(NLAYERS+2*NFLAYERS)*(NLAYERS+2*NFLAYERS)];
    double S_pre[(NLAYERS+2*NFLAYERS)*(NLAYERS+2*NFLAYERS)];
    double S0[(NLAYERS+2*NFLAYERS)*(NLAYERS+2*NFLAYERS)];
#endif
    double y_ini[NFREQS*2];
    int up = 0;
    int s7 = AVERAGE, s7c = 0;
    double mesv[NFREQS*2];
    double alta = 0;
    double vda = 0;
    int ft =1;
    int data_cntr = 0;
    double dpth[NLAYERS+NFLAYERS],d=FIRST_THK;
    int nlay = NLAYERS+NFLAYERS;
    int i;
    double weight = 1000.;

    for(int i=0;i<charge*NPLAYERS+NLAYERS+NFLAYERS;i++)
        upper[i] = (i<NLAYERS+NFLAYERS+NPLAYERS) ? MAX_RES : 1.;
    for(int i=0;i<charge*NPLAYERS+NLAYERS+NFLAYERS;i++)
        lower[i] = (i<NLAYERS+NFLAYERS+NPLAYERS) ? MIN_RES : .000001;
    for(int i=charge*NPLAYERS+NLAYERS+NFLAYERS;i<charge*NPLAYERS+NLAYERS+2*NFLAYERS;i++)
        upper[i] = 100.;
    for(int i=charge*NPLAYERS+NLAYERS+NFLAYERS;i<charge*NPLAYERS+NLAYERS+2*NFLAYERS;i++)
        lower[i] = .1;

    memset(mesv,0,sizeof(mesv));

    memset(time, 0, sizeof(time));

    // Layers thicknesses are fixed here !!!
    for(int i=0;i<NFLAYERS;i++) dpth[i] = 10.;
    for(int i=NFLAYERS;i<NLAYERS+NFLAYERS-1;i++,d*=DSTEP) dpth[i] = d;

    data = buf + 13; // it's an offset to skip time in the string

    FILE *fin  = fopen(argv[1],"rt");
    FILE *fout = fopen(argv[2],"wt");

    //double primField(double hd, double vd)
    double mom = 10000/primField(28.8,27.8);
//    double mom = 10000/primField(28.28,28.33);
    geo.prim = 1./mom; //5099041845; // Equator's primary to get 10000 A/m in FD!!! - greater than 5.e9 Am2
                       //5563098240;//!!!!!! For 2200:5 in Croatia!!!! // primField(28.21,28.54)/.92;
                       //6489173050;//!!!!!! For 2110:2 in Croatia!!!! // primField(29.3,27.8)/.82;

    memset(geo.w,0,sizeof(geo.w));
    //21 Frequencies
    // weghts for the measurements are 1/RMS
    // Re
    geo.w[1] = 1./.27;    //
    geo.w[3] = 1./.24;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./.26;
    geo.w[7] = 1./.28;
    geo.w[9] = 1./.3;
    geo.w[11] = 1./.89;
    geo.w[13] = 1./.58;
    geo.w[15] = 1./0.47;
    geo.w[17] = 1./.55;
    geo.w[19] = 1./.35;
    geo.w[21] = 1./.38;
    geo.w[23] = 1./.91;
    geo.w[25] = 1./.78;
    geo.w[27] = 1./.86;
    geo.w[29] = 1./.92;

    for(int i=0;i<NFREQS-1;i++)
        geo.w[i*2+1] = sqrt(geo.w[i*2+1]*geo.w[i*2+1]+geo.w[(i+1)*2+1]*geo.w[(i+1)*2+1]);

    // Im
    geo.w[0] = 1./.21;
    geo.w[2] = 1./.1;
    geo.w[4] = 1./.07;
    geo.w[6] = 1./.12;
    geo.w[8] = 1./.35;
    geo.w[10] = 1/.82;
    geo.w[12] = 1./.41;
    geo.w[14] = 1./.53;
    geo.w[16] = 1./.28;
    geo.w[18] = 1./.42;
    geo.w[20] = 1./.15;
    geo.w[22] = 1./.27;
    geo.w[24] = 1./.82;
    geo.w[26] = 1./.89;
    geo.w[28] = 1./.94;

     #if 0
    //21 Frequencies
    // weghts for the measurements are 1/RMS
    // Re
    geo.w[1] = 1./.27;    //
    geo.w[3] = 1./.24;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./.26;
    geo.w[7] = 1./.28;
    geo.w[9] = 1./.3;
    geo.w[11] = 1./.49;
    geo.w[13] = 1./.58;
    geo.w[15] = 1./1.07;
    geo.w[17] = 1./.55;
    geo.w[19] = 1./.35;
    geo.w[21] = 1./.38;
    geo.w[23] = 1./.91;
    geo.w[25] = 1./.48;
    geo.w[27] = 1./.46;
    geo.w[29] = 1./.52;
    geo.w[31] = 1./.45;
    geo.w[33] = 1./.49;
    geo.w[35] = 1./.25;
    geo.w[37] = 1./.65;
    geo.w[39] = 1./.67;
    geo.w[41] = 1./.92;

    for(int i=0;i<NFREQS-1;i++)
        geo.w[i*2+1] = sqrt(geo.w[i*2+1]*geo.w[i*2+1]+geo.w[(i+1)*2+1]*geo.w[(i+1)*2+1]);
    geo.w[41] = 1./20.;//

    // Im
    geo.w[0] = 1./.21;
    geo.w[2] = 1./.1;
    geo.w[4] = 1./.07;
    geo.w[6] = 1./.12;
    geo.w[8] = 1./.35;
    geo.w[10] = 1/.42;
    geo.w[12] = 1./.41;
    geo.w[14] = 1./.53;
    geo.w[16] = 1./.28;
    geo.w[18] = 1./.42;
    geo.w[20] = 1./.15;
    geo.w[22] = 1./.27;
    geo.w[24] = 1./.82;
    geo.w[26] = 1./.39;
    geo.w[28] = 1./.44;
    geo.w[30] = 1./.59;
    geo.w[32] = 1./.38;
    geo.w[34] = 1./.17;
    geo.w[36] = 1./.44;
    geo.w[38] = 1./.6;
    geo.w[40] = 1./.9;
#endif
#if 0
    //10 frequencies
    // weghts for the measurements are 1/RMS
    // Re
    geo.w[1] = 1./2.;    //
    geo.w[3] = 1./2.;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./0.14;
    geo.w[7] = 1./0.11;
    geo.w[9] = 1./0.15;
    geo.w[11] = 1./0.29;
    geo.w[13] = 1./.24;
    geo.w[15] = 1./.77;
    geo.w[17] = 1./1.11;
    geo.w[19] = 1./1.11;

    //geo.w[11] = 1./0.23;
    //geo.w[13] = 1./0.21;
    //geo.w[15] = 1./0.49;
    //geo.w[17] = 1./0.91;
    //geo.w[19] = 1./0.91;

    // Im
    geo.w[0] = 1./0.09;
    geo.w[2] = 1./0.07;
    geo.w[4] = 1./0.07;
    geo.w[6] = 1./0.04;
    geo.w[8] = 1./0.09;
    geo.w[10] = 1/0.19;
    geo.w[12] = 1./0.16;
    geo.w[14] = 1./0.52;
    geo.w[16] = 1./0.72;
    //geo.w[18] = 1./2;
    geo.w[18] = 1./0.72;
#endif
#if 0
    //17 frequencies
    // weghts for the measurements are 1/RMS
    // Re
    geo.w[1] = 1./2.;    //
    geo.w[3] = 1./2.;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./0.14;
    geo.w[7] = 1./0.14;
    geo.w[9] = 1./0.16;
    geo.w[11] = 1./0.19;
    geo.w[13] = 1./0.15;
    geo.w[15] = 1./0.14;
    geo.w[17] = 1./0.16;
    geo.w[19] = 1./0.11;
    geo.w[21] = 1./0.14;
    geo.w[23] = 1./0.16;
    geo.w[25] = 1./0.15;
    geo.w[27] = 1./0.23;
    geo.w[29] = 1./0.21;
    geo.w[31] = 1./0.49;
    geo.w[33] = 1./0.91;

    // Im
    geo.w[0] = 1./0.09;
    geo.w[2] = 1./0.07;
    geo.w[4] = 1./0.07;
    geo.w[6] = 1./0.09;
    geo.w[8] = 1./0.13;
    geo.w[10] = 1/0.15;
    geo.w[12] = 1./0.06;
    geo.w[14] = 1./0.06;
    geo.w[16] = 1./0.09;
    geo.w[18] = 1./0.04;
    geo.w[20] = 1./0.12;
    geo.w[22] = 1./0.13;
    geo.w[24] = 1./0.09;
    geo.w[26] = 1./0.19;
    geo.w[28] = 1./0.16;
    geo.w[30] = 1./0.52;
    geo.w[32] = 1./0.72;
#endif
   #if 0
    // weghts for the measurements are 1/RMS
    // Re
    geo.w[1] = 1./2.;    //
    geo.w[3] = 1./2.;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./0.14;
    geo.w[7] = 1./0.14;
    geo.w[9] = 1./0.16;
    geo.w[11] = 1./0.19;
    geo.w[13] = 1./0.15;
    geo.w[15] = 1./0.14;
    geo.w[17] = 1./0.16;
    geo.w[19] = 1./0.11;
    geo.w[21] = 1./0.14;
    geo.w[23] = 1./0.16;
    geo.w[25] = 1./0.15;
    geo.w[27] = 1./0.26;
    geo.w[29] = 1./0.23;
    geo.w[31] = 1./0.28;
    geo.w[33] = 1./0.29;
    geo.w[35] = 1./0.32;
    geo.w[37] = 1./0.21;
    geo.w[39] = 1./0.49;
    geo.w[41] = 1./0.56;
    geo.w[43] = 1./0.46;
    geo.w[45] = 1./0.91;

    // Im
    geo.w[0] = 1./0.09;
    geo.w[2] = 1./0.07;
    geo.w[4] = 1./0.07;
    geo.w[6] = 1./0.09;
    geo.w[8] = 1./0.13;
    geo.w[10] = 1/0.15;
    geo.w[12] = 1./0.06;
    geo.w[14] = 1./0.06;
    geo.w[16] = 1./0.09;
    geo.w[18] = 1./0.04;
    geo.w[20] = 1./0.12;
    geo.w[22] = 1./0.13;
    geo.w[24] = 1./0.09;
    geo.w[26] = 1./0.24;
    geo.w[28] = 1./0.19;
    geo.w[30] = 1./0.25;
    geo.w[32] = 1./0.26;
    geo.w[34] = 1./0.27;
    geo.w[36] = 1./0.16;
    geo.w[38] = 1./0.52;
    geo.w[40] = 1./0.45;
    geo.w[42] = 1./0.37;
    geo.w[44] = 1./0.72;
#endif
#if 0
    geo.w[1] = 1./2.;    //
    geo.w[3] = 1./2.;    // base frequency for Re. All lower Re are fake.
    geo.w[5] = 1./0.10;
    geo.w[7] = 1./0.09;
    geo.w[9] = 1./0.11;
    geo.w[11] = 1./5.10; // 50Hz Noises
    geo.w[13] = 1./0.15;
    geo.w[15] = 1./0.10;
    geo.w[17] = 1./0.09;
    geo.w[19] = 1./0.12;
    geo.w[21] = 1./0.10;
    geo.w[23] = 1./0.10;
    geo.w[25] = 1./0.10;
    geo.w[27] = 1./0.12;
    geo.w[29] = 1./0.09;
    geo.w[31] = 1./0.08;
    geo.w[33] = 1./0.08;
    geo.w[35] = 1./0.09;
    geo.w[37] = 1./0.23;
    geo.w[39] = 1./0.10;
    geo.w[41] = 1./0.17;
    geo.w[43] = 1./0.20;
    geo.w[45] = 1./0.17;
    geo.w[47] = 1./0.44;
    geo.w[49] = 1./0.31;
    geo.w[51] = 1./0.35;
    // Im
    geo.w[0] = 1./0.10;
    geo.w[2] = 1./0.05;
    geo.w[4] = 1./0.05;
    geo.w[6] = 1./0.08;
    geo.w[8] = 1./0.08;
    geo.w[10] = 1./5.09; // 50 Hz
    geo.w[12] = 1./0.12;
    geo.w[14] = 1./0.05;
    geo.w[16] = 1./0.05;
    geo.w[18] = 1./0.08;
    geo.w[20] = 1./0.07;
    geo.w[22] = 1./0.06;
    geo.w[24] = 1./0.08;
    geo.w[26] = 1./0.09;
    geo.w[28] = 1./0.05;
    geo.w[30] = 1./0.06;
    geo.w[32] = 1./0.07;
    geo.w[34] = 1./0.07;
    geo.w[36] = 1./0.25;
    geo.w[38] = 1./0.09;
    geo.w[40] = 1./0.16;
    geo.w[42] = 1./0.17;
    geo.w[44] = 1./0.14;
    geo.w[46] = 1./0.31;
    geo.w[48] = 1./0.23;
    geo.w[50] = 1./0.20;
#endif
    for(int lr = 0;lr<charge*NPLAYERS+NLAYERS+NFLAYERS;lr++)
        rho[lr] = (lr<NPLAYERS+NLAYERS+NFLAYERS)?RES_INI:.0001;

    memset(x_ini,0,sizeof(x_ini));
    memset(S_ini,0,sizeof(S_ini));
    memset(buf,0,sizeof(buf));
    for(int i=0;i<charge*NPLAYERS+NLAYERS+NFLAYERS;i++)
        x_ini[i] = rho[i];

    for(i=0;fgets(buf,2000,fin);i++) {
        if(buf[0]=='/') {  // reading comments
            memset(buf,0,sizeof(buf));
            continue;
        }

        if(buf[0]=='L' || buf[0]=='B' || buf[0]=='T') { // reading Lines or Base or Trends
            fputs(buf,fout);
            printf("%s",buf);
            for(int lr = 0;lr<charge*NPLAYERS+NLAYERS+NFLAYERS;lr++)
                rho[lr] = (lr<NPLAYERS+NLAYERS+NFLAYERS)?RES_INI:.0001;
            memset(buf,0,sizeof(buf));
            ft = 1;
            continue;
        }

        memcpy(time,buf,sizeof(time)-1); // reading time

        if(strstr(data,"*")) { // all the nonmeasurements are skipped
            memset(buf,0,sizeof(buf));
            continue;
        }

#if 0
	// Reading data file after Time
        //           x   y   alt hd  vd  re0 re1 re2 re3 re4 re5 re6 re7
        //           re8 re9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19 r20
        //           r21 r22 r23 r24 r25 im0 im1 im2 im3 im4 im5 im6 im7
        //           im8 im9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20
        //           i21 i22 i23 i24 i25
        sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf",
               &geo.coord,&geo.coord_y,&geo.alt,&geo.hor_dist,&geo.ver_dist,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],               &y_mes[21],               &y_mes[23],
               &y_mes[25],               &y_mes[27],               &y_mes[29],
               &y_mes[31],               &y_mes[33],               &y_mes[35],
               &y_mes[37],               &y_mes[39],               &y_mes[41],
               &y_mes[43],               &y_mes[45],               &y_mes[47],
               &y_mes[49],               &y_mes[51],
               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18],               &y_mes[20],               &y_mes[22],
               &y_mes[24],               &y_mes[26],               &y_mes[28],
               &y_mes[30],               &y_mes[32],               &y_mes[34],
               &y_mes[36],               &y_mes[38],               &y_mes[40],
               &y_mes[42],               &y_mes[44],               &y_mes[46],
               &y_mes[48],               &y_mes[50]
               );
        // Reading data file after Time
            //           x   y   alt hd  vd  re0 re1 re2 re3 re4 re5 re6 re7
            //           re8 re9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19 r20
            //           r21 r22 im0 im1 im2 im3 im4 im5 im6 im7
            //           im8 im9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20
            //           i21 i22
        sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf",
               &geo.coord,&geo.coord_y,&geo.alt,&geo.hor_dist,&geo.ver_dist,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],               &y_mes[21],               &y_mes[23],
               &y_mes[25],               &y_mes[27],               &y_mes[29],
               &y_mes[31],               &y_mes[33],               &y_mes[35],
               &y_mes[37],               &y_mes[39],               &y_mes[41],
               &y_mes[43],               &y_mes[45],
               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18],               &y_mes[20],               &y_mes[22],
               &y_mes[24],               &y_mes[26],               &y_mes[28],
               &y_mes[30],               &y_mes[32],               &y_mes[34],
               &y_mes[36],               &y_mes[38],               &y_mes[40],
               &y_mes[42],               &y_mes[44]
               );
        memset(buf,0,sizeof(buf));
#endif
#if 0
        //for 17 frequencies
        // Reading data file after Time
            //           x   y   alt hd  vd
            //           re0 re1 re2 re3 re4 re5 re6 re7
            //           re8 re9 r10 r11 r12 r13 r14 r15 r16
            //           im0 im1 im2 im3 im4 im5 im6 im7
            //           im8 im9 i10 i11 i12 i13 i14 i15 i16
        sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &geo.coord,&geo.coord_y,&geo.alt,&geo.hor_dist,&geo.ver_dist,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],               &y_mes[21],               &y_mes[23],
               &y_mes[25],               &y_mes[27],               &y_mes[29],
               &y_mes[31],               &y_mes[33],
               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18],               &y_mes[20],               &y_mes[22],
               &y_mes[24],               &y_mes[26],               &y_mes[28],
               &y_mes[30],               &y_mes[32]
               );
        memset(buf,0,sizeof(buf));
#endif
#if 0
        //for 10 frequencies
        // Reading data file after Time
            //           x   y   alt hd  vd
            //           re0 re1 re2 re3 re4 re5 re6 re7
            //           re8 re9
            //           im0 im1 im2 im3 im4 im5 im6 im7
            //           im8 im9 i10
        sscanf(data,"%lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ",
               &geo.coord,&geo.coord_y,&geo.alt,&geo.hor_dist,&geo.ver_dist,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],
               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18]
               );
        memset(buf,0,sizeof(buf));
#endif
   #if 0
        //for 21 frequencies
        // Reading data file after Time
            //           hd  vd  alt re0 re1 re2 re3 re4 re5 re6 re7
            //           re8 re9 r10 r11 r12 r13 r14 r15 r16 r17 r18 r19 r20
            //           im0 im1 im2 im3 im4 im5 im6 im7
            //           im8 im9 i10 i11 i12 i13 i14 i15 i16 i17 i18 i19 i20

        sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf",
               &geo.hor_dist, &geo.ver_dist, &geo.alt,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],               &y_mes[21],               &y_mes[23],
               &y_mes[25],               &y_mes[27],               &y_mes[29],
               &y_mes[31],               &y_mes[33],               &y_mes[35],
               &y_mes[37],               &y_mes[39],               &y_mes[41],
               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18],               &y_mes[20],               &y_mes[22],
               &y_mes[24],               &y_mes[26],               &y_mes[28],
               &y_mes[30],               &y_mes[32],               &y_mes[34],
               &y_mes[36],               &y_mes[38],               &y_mes[40]
               );
        memset(buf,0,sizeof(buf));
     #endif
        //for 15 frequencies
        // Reading data file after Time
            //           hd  vd  alt re0 re1 re2 re3 re4 re5 re6 re7
            //           re8 re9 r10 r11 r12 r13 r14
            //           im0 im1 im2 im3 im4 im5 im6 im7
            //           im8 im9 i10 i11 i12 i13 i14

        sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \
                     %lf %lf %lf %lf %lf %lf %lf %lf %lf",
               &geo.hor_dist, &geo.ver_dist, &geo.alt,
               &y_mes[1],                &y_mes[3],                &y_mes[5],
               &y_mes[7],                &y_mes[9],                &y_mes[11],
               &y_mes[13],               &y_mes[15],               &y_mes[17],
               &y_mes[19],               &y_mes[21],               &y_mes[23],
               &y_mes[25],               &y_mes[27],               &y_mes[29],

               &y_mes[0],                &y_mes[2],                &y_mes[4],
               &y_mes[6],                &y_mes[8],                &y_mes[10],
               &y_mes[12],               &y_mes[14],               &y_mes[16],
               &y_mes[18],               &y_mes[20],               &y_mes[22],
               &y_mes[24],               &y_mes[26],               &y_mes[28]
               );
        memset(buf,0,sizeof(buf));
        data_cntr++;
        //if(data_cntr<1400) continue;

        geo.alt += DA;

    // averaging for AVERAGE samples
        for(int i=0;i<2*NFREQS;i++)
            mesv[i] += y_mes[i];
        alta += geo.alt;
        vda += geo.ver_dist;
        s7c++;
        if(--s7)
            continue;
        for(int i=0;i<2*NFREQS;i++)
            y_mes[i] = mesv[i]/s7c;
        geo.alt = alta/s7c;
        geo.ver_dist = vda/s7c;
        s7c = 0;
        s7 = AVERAGE;
        memset(mesv,0,sizeof(mesv));
        alta = 0;
        vda = 0;

	// excluding all negatives
        for(int i=0;i<2*NFREQS;i++)
           y_mes[i] = fabs(y_mes[i]);

    // It's a cricial point: we use differences of Inphase components instead of Inphase components!!!
        if(1) {//NFREQS == 21) {
            for(i=0;i<NFREQS;i++) {
                if(i<NFREQS-1)
                    y_mes[2*i+1] = y_mes[2*(i+1)+1] - y_mes[2*i+1];
                else
                    y_mes[2*i+1] = y_mes[2*(i-1)+1];
            }
        }

    // calculating of a half-space for the lowerest frequency!!!
    // itereative inversion for fixed layers
        double res = -1;
        double rho_ini = RES_INI;
        int itr;
        for (itr = 0;itr < MAX_ITER; itr++) {
            double d_ini = ERR_INI;
            int up = 0;
            flinversion(geo,1,1,&rho_ini,dpth,y_ini,y_mes,&res,&up,&d_ini);
            if(sqrt(res) <STOP_VAL) break;
            if(up) break;
        }

        res = sqrt(res);

        printf("%s %f %d %f ",time, res, itr, rho_ini);
        printf("\n");


	//inversion
        int iter = 0;
        //double res1;

#ifdef CHARGEABILITY
        double Sr[(2*NPLAYERS+NLAYERS+2*NFLAYERS)*(2*NPLAYERS+NLAYERS+2*NFLAYERS)];
#else
        double Sr[(NLAYERS+2*NFLAYERS)*(NLAYERS+2*NFLAYERS)];
#endif
        res = -1;
        up = 0;

        if(ft) { // for each new line restart KF
            memset(S_ini,0,sizeof(S_ini));
            // thickness
            for(int i=charge*NPLAYERS+NFLAYERS+nlay-1;i>=charge*NPLAYERS+nlay;i--) {
                S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = ERR_INI;
                x_ini[i] =DEP_INI;
            }
            // chargeability
            for(int i=charge*NPLAYERS+nlay-1;i>=NPLAYERS+nlay;i--) {
                if(i==charge*NPLAYERS+nlay-1)
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = ERR_INI;
                else {
                    S_ini[i+1+(charge*NPLAYERS+nlay)*i] =
                            COR_INI*ERR_INI/S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*(i+1)];
                    S_ini[i+(charge*NPLAYERS+nlay)*i] =
                     sqrt(ERR_INI*ERR_INI-S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]*
                            S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]);
                }
                x_ini[i] =.0001;
            }
            // ch.resistivity
            for(int i=NPLAYERS+nlay-1;i>=nlay;i--) {
                if(i==NPLAYERS+nlay-1)
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = ERR_INI;
                else {
                    S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i] =
                            COR_INI*ERR_INI/S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*(i+1)];
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] =
                     sqrt(ERR_INI*ERR_INI-S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]*
                            S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]);
                }
                x_ini[i] = rho_ini;//RES_INI;//*.5;
            }
            // ordinary resistivity
            for(int i=nlay-1;i>=0;i--) {
                if(i==nlay-1)
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = ERR_INI;
                else {
                    S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i] =
                            COR_INI*ERR_INI/S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*(i+1)];
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] =
                     sqrt(ERR_INI*ERR_INI-S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]*
                            S_ini[i+1+(charge*NPLAYERS+nlay+NFLAYERS)*i]);
                }
                x_ini[i] = rho_ini;//RES_INI;//*.5;
                //if(i<NFLAYERS) x_ini[i] *= .5;
            }

#if 0
            double PP[(charge*NPLAYERS+nlay)*(charge*NPLAYERS+nlay)];
            for(int i=0;i<(charge*NPLAYERS+nlay);i++) {
                for (int j=0;j<(charge*NPLAYERS+nlay);j++) {
                    PP[i*(charge*NPLAYERS+nlay)+j] = 0;
                    for (int k=0;k<(charge*NPLAYERS+nlay);k++) {
                        PP[i*(charge*NPLAYERS+nlay)+j]+=
                                S_ini[i*(charge*NPLAYERS+nlay)+k]*S_ini[j*(charge*NPLAYERS+nlay)+k];
                    }
                    printf("%.7f ",S_ini[i*(charge*NPLAYERS+nlay)+j]);
                }
                printf("    %f\n",x_ini[i]);
            }
#endif


#if 0
            for(int i=charge*nlay-1;i>=0;i--) {
                if(i==charge*nlay-1)
                    S_ini[i+charge*nlay*i] = ERR_INI;
                else {
                    S_ini[i+1+charge*nlay*i] = COR_INI*ERR_INI/S_ini[i+1+charge*nlay*(i+1)];
                    S_ini[i+charge*nlay*i] = sqrt(ERR_INI*ERR_INI-S_ini[i+1+charge*nlay*i]*S_ini[i+1+charge*nlay*i]);
                }
                x_ini[i] =(i<2*NLAYERS)?rho_ini/*RES_INI*/:.0001;
            }
#endif
            ft = !ft;
/*
            double P[NLAYERS*NLAYERS];
            for(int i=0;i<NLAYERS;i++) {
                for(int j=0;j<NLAYERS;j++) {
                    P[i+j*NLAYERS] = 0;
                    for(int k=0;k<NLAYERS;k++)
                        P[i+j*NLAYERS] += S_ini[k+i*NLAYERS]*S_ini[k+j*NLAYERS];
                }
            }
*/
            memcpy(S0,S_ini,sizeof(Sr));
        } else { // restart S for the next point. x_ini is tied to the initial resistivity
            memset(S_ini,0,sizeof(S_ini));
            for(int i=(charge*NPLAYERS+nlay+NFLAYERS)-1;i>=charge*NPLAYERS+nlay;i--) {
                S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = S0[i+(charge*NPLAYERS+nlay+NFLAYERS)*i];
                x_ini[i] = .97*x_ini[i]+.03*DEP_INI;
            }
            for(int i=(charge*NPLAYERS+nlay)-1;i>=0;i--) {
                S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i] = S0[i+(charge*NPLAYERS+nlay+NFLAYERS)*i];
                if(i<(charge*NPLAYERS+nlay+NFLAYERS)-1)
                    S_ini[i+(charge*NPLAYERS+nlay+NFLAYERS)*i+1] = S0[i+(charge*NPLAYERS+nlay+NFLAYERS)*i+1];
                if(NPLAYERS) {
                    x_ini[i] =
                        (i<NLAYERS+NPLAYERS+NFLAYERS)?
//                            sqrt(sqrt(rho_ini/*RES_INI*/*x_ini[i])*x_ini[i]):
                            (rho_ini*(1-1./weight)+x_ini[i]/weight):
                            sqrt(sqrt(.0001*x_ini[i])*x_ini[i]);
                } else {
                    x_ini[i] = (rho_ini*(1-1./weight)+x_ini[i]/weight);
                }
            }
        }
        memcpy(Sr,S_ini,sizeof(Sr));

	// itereative inversion for fixed layers
        for (iter = 0;iter < MAX_ITER; iter++) {
            memcpy(S_ini,Sr,sizeof(Sr));
            flinversion(geo,NFREQS,nlay,x_ini,dpth,y_ini,y_mes,&res,&up,S_ini);
            if(sqrt(res) <STOP_VAL) break;
            if(up) break;
        }
        
        res = sqrt(res);
        weight = (res>1)?res:1.;

        printf("%s %f %d ",time, res, iter);
        printf("\n");

        fprintf(fout, "%s %f %d ",time, res, iter);
        for(int i=0;i<(charge*NPLAYERS+nlay);i++) {
            if(i<NPLAYERS+nlay)
                fprintf(fout,"%.3f ",x_ini[i]);
            else
                fprintf(fout,"%.3f ",x_ini[i]*1000.);
            printf("%f ",x_ini[i]);
        }
        double depp = 0;
        for(int i=0;i<nlay-1;i++) {
            fprintf(fout,"%.3f ",depp+dpth[i]*.5);
            depp += dpth[i];
        }
        printf("%f ",dpth[0]);

        fprintf(fout,"%.3f %.3f ",depp+dpth[nlay-2]*.5,rho_ini);
        printf("\n");

        for(int i=0;i<(charge*NPLAYERS+nlay+NFLAYERS);i++) {
            double v = 0, v1 = 0;
            for (int j=i;j<(charge*NPLAYERS+nlay+NFLAYERS);j++) {
                v+= S_ini[i*(charge*NPLAYERS+nlay+NFLAYERS)+j]*S_ini[i*(charge*NPLAYERS+nlay+NFLAYERS)+j];
                v1+= Sr[i*(charge*NPLAYERS+nlay+NFLAYERS)+j]*Sr[i*(charge*NPLAYERS+nlay+NFLAYERS)+j];
            }
            fprintf(fout,"%.7f ",1 - sqrt(v/v1));
            printf("%f ",1 - sqrt(v/v1));
        }

        fprintf(fout,"\n");
        printf("\n");
        continue;
        // here ends the routine
    }
    printf("\n end");
    return 0;
}


