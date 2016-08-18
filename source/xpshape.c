/*

  simulate NMR peakshape for two-site chemical exchange

  written by Sung-Hun Bae
  Nov. 17, 2008

*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <gsl/gsl_fft_complex.h>

#define TINY  1.0e-12

using namespace std;

enum SITE_LABEL { A, B, N_SITE };
enum QUADRATURE { REAL, IMAG };
typedef struct {
  double Omegaa;// frequency of A site (rad/s)
  double Omegab;// frequency of B site (rad/s)
  double kab;	// reaction rate constant, A -> B (1/s)
  double kba;	// reaction rate constant, B -> A (1/s)
  double R20a;	// transverse relaxation rate of A site (1/s)
  double R20b;	// transverse relaxation rate of B site	(1/s)
  double M0a;	// initial magnetization
  double M0b;	// initial magnetization
  SITE_LABEL site;	// A or B
  QUADRATURE quad;	// REAL or IMAG
  } two_site_exchange;

//
// Magnetization (two-site exchange)
//

void N (double t, double &R, double &I, void * params) {
  two_site_exchange * x = (two_site_exchange *) params;
  complex<double> c1,c2,c3,lp,lm,a11,a12,a21,a22;
  complex<double> Omegaa = sqrt(complex<double>(-1.))*x->Omegaa;
  complex<double> Omegab = sqrt(complex<double>(-1.))*x->Omegab;
  double mag;
  c1 = (-Omegaa - Omegab + x->R20a + x->R20b + x->kab + x->kba);
  c2 = (-Omegaa + Omegab + x->R20a - x->R20b + x->kab - x->kba);
  c3 = sqrt(c2*c2+4.*x->kab*x->kba);
  lp = 0.5*(c1+c3);
  lm = 0.5*(c1-c3);
  if (x->site == A) {
    a11 = 0.5*((1.-c2/c3)*exp(-t*lm)+(1.+c2/c3)*exp(-t*lp)); 
    a12 = x->kba/c3*(exp(-t*lm)-exp(-t*lp));
    R += real(a11*x->M0a+a12*x->M0b);
    I += imag(a11*x->M0a+a12*x->M0b);
    } 
  if (x->site == B) {
    a22 = 0.5*((1.+c2/c3)*exp(-t*lm)+(1.-c2/c3)*exp(-t*lp)); 
    a21 = x->kab/c3*(exp(-t*lm)-exp(-t*lp));
    R += real(a21*x->M0a+a22*x->M0b);
    I += imag(a21*x->M0a+a22*x->M0b);
    }
}

int main (int argc, char *argv[])
{

  if (argc < 7 || argc > 17) {
    printf("\t\t\t\t\t\t\tSung-Hun Bae 2008.11\n");
    printf("   ");
    printf("Usage: xpshape [-model AB|ACB]\n");
    printf("   ");
    printf("               [-dw #(Hz)] [-kab #(1/s)] [-kba #(1/s)]\n");
    printf("   ");
    printf("               [-R20a #(1/s)] [-R20b #(1/s)]\n");
    printf("   ");
    printf("               [-a0 conc(M)] [-c0 conc(M)]\n");
    printf("   ");
    printf("               [-swh Hz(600)] [-fid pts(8192)]\n\n");
    printf("   ");
    printf("simulate NMR peakshape for two-site chemical exchange\n\n");
    printf("   ");
    printf("-model,AB  : A = B     (dw,kab,kba,R20a,R20b required)\n");
    printf("   ");
    printf("       ACB : A + C = B (dw,kab,kba,R20a,R20b,a0,c0 required)\n");
    printf("\n");
    printf("   ");
    printf("-dw,  (A.frequency= -dw/2, B.frequency= +dw/2)\n");
    printf("   ");
    printf("-R20a,(exchange-free R2 of A state)\n"); 
    printf("   ");
    printf("-R20b,(exchange-free R2 of B state)\n"); 
    printf("   ");
    printf("-kab, (rate constant of A->B)\n"); 
    printf("   ");
    printf("-kba, (rate constant of B->A)\n"); 
    printf("   ");
    printf("-a0,  (initial concentration of A)\n"); 
    printf("   ");
    printf("-c0,  (initial concentration of C)\n"); 
    printf("   ");
    printf("-swh, (spectrum = -swh/2 ... +swh/2)\n");
    printf("   ");
    printf("-fid, (must be ...,1024,2048,4096,8192,16384,32768,65536,..)\n");
    printf("\n");
    
    exit (1);
    }

  // default settings

  double SWH = 600;		  // (Hz)
  unsigned int FID_SIZE = 8192;	  // (complex data points)


  // receiving parameters

  unsigned int i = 1;

  string model;
  double dw,R20a,R20b,kab,kba,concA0,concC0;
  bool dw_,R20a_,R20b_,kab_,kba_,concA0_,concC0_;
  double pb,Kd,concA,concB,concC,C,u,v;

  dw_ = R20a_ = R20b_ = kab_ = kba_ = concA0_ = concC0_ = false;

  while ((argc - i) > 0) {
    if ((argc-i) > 0 && !strcmp(argv[i],"-model")) {
      model = string(argv[++i]); 
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-dw")) {
      dw = 2.*M_PI*atof(argv[++i]); // (Hz) -> (rad/s)
      dw_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-kab")) {
      kab = atof(argv[++i]);
      kab_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-kba")) {
      kba = atof(argv[++i]);
      kba_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-R20a")) {
      R20a = atof(argv[++i]);
      R20a_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-R20b")) {
      R20b = atof(argv[++i]);
      R20b_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-a0")) {
      concA0 = atof(argv[++i]);
      concA0_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-c0")) {
      concC0 = atof(argv[++i]);
      concC0_ = true;
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-swh")) {
      SWH = atof(argv[++i]);
      i++;
      }
    if ((argc-i) > 0 && !strcmp(argv[i],"-fid")) {
      FID_SIZE = atoi(argv[++i]);
      i++;
      }
    }

  double fid[2*FID_SIZE];	    // complex points
  double dwell = 1./(SWH);	    // (sec)
  double norm = 1./sqrt(FID_SIZE);  // scaling factor

  if (model == "AB" && dw_ && R20a_ && R20b_ && kab_ && kba_) {
    // dw, R20a, R20b, kab, kba
    Kd = kba/kab;
    pb = 1./(1.+Kd);
    C = 1.0;
    }
  else if (model == "ACB" && dw_ && R20a_ && R20b_ && kab_ && kba_ 
    && concA0_ && concC0_) {
    // dw, R20a, R20b, kab(=kon), kba(=koff), concA0, concC0
    kab = (concC0 > TINY ? kab : 0.0);
    kba = (concC0 > TINY ? kba : 0.0);
    Kd = ((kab > TINY) && (kba > TINY) ? kba/kab : 0.0);
    u = 0.5*(concA0+concC0+Kd);
    v = (concC0 > TINY ? (u-sqrt(u*u-concA0*concC0)) : 0.0);
    concA = (v > TINY ? concA0-v : concA0);
    concB = (v > TINY ? v : 0.0);
    concC = (v > TINY ? concC0 - v : concC0);
    pb = (concC0 > TINY ? kab*concC/(kab*concC+kba) : 0.0);
    C = concC;
    }
  else {
    printf("error: missing parameter(s)\n");
    exit(1);
    }

  two_site_exchange x = {-dw/2.,dw/2.,kab*C,kba,R20a,R20b,(1.-pb),pb,A,REAL};

  // header information
  printf("# peak shape simulation ");
  if (model == "AB")  printf("[A = B]\n");
  if (model == "ACB") printf("[A + C = B]\n");

  printf("# dw     = %g (rad/s)\n", dw);
  printf("# kab    = %g (1/s)\n",kab);
  printf("# kba    = %g (1/s)\n",kba);
  printf("# R20a   = %g (1/s)\n",R20a);
  printf("# R20b   = %g (1/s)\n",R20b);
  if (model == "AB") {
  printf("# pb     = %g\n", pb);
  }
  if (model == "ACB") {
  printf("# kon    = %g (1/s)\n",kab);
  printf("# koff   = %g (1/s)\n",kba);
  printf("# Kd     = %g (M)\n",Kd);
  printf("# pb     = %g\n",pb);
  printf("# [A]ini = %g (M)\n",concA0);
  printf("# [C]ini = %g (M)\n",concC0);
  printf("# [A]    = %g (M)\n",concA);
  printf("# [C]    = %g (M)\n",concC);
  printf("# [B]    = %g (M)\n",concB);
  }
  printf("# FID    = %d (complex points)\n",FID_SIZE);
  printf("# SWH    = %g (Hz)\n",SWH);
  printf("# X-axis unit: (Hz)\n");
  
  // creating complex data points
  for (i = 0; i < FID_SIZE; i++) {
    fid[2*i] = fid[2*i+1] = 0.0;
    x.site = A;
    N (i*dwell, fid[2*i], fid[2*i+1], &x);  // A
    x.site = B;
    N (i*dwell, fid[2*i], fid[2*i+1], &x);  // B
    }

  // Fast Fourier Transform (FFT)
  gsl_fft_complex_radix2_forward (fid, 1, FID_SIZE);

  // spectrum 
  for (i = FID_SIZE/2; i < FID_SIZE; i++)
    printf("%g %g\n",(-0.5*SWH + SWH*(i-FID_SIZE/2)/FID_SIZE),norm*fid[2*i]);
  for (i = 0; i <= FID_SIZE/2; i++)
    printf("%g %g\n",(-0.5*SWH + SWH*(i+FID_SIZE/2)/FID_SIZE),norm*fid[2*i]);
  printf("&\n");

  return 0;
}

/*

  For physical applications it is important to remember that the index 
  appearing in the DFT does not correspond directly to a physical frequency.
  If the time-step of the DFT is \Delta then the frequency-domain includes 
  both positive and negative frequencies, ranging from -1/(2\Delta) through 
  0 to +1/(2\Delta). 
  The positive frequencies are stored from the beginning of the array up to 
  the middle, and the negative frequencies are stored backwards from the 
  end of the array.
  
     index    z               x = FFT(z)
     
     0        z(t = 0)        x(f = 0)
     1        z(t = 1)        x(f = 1/(N Delta))
     2        z(t = 2)        x(f = 2/(N Delta))
     .        ........        ..................
     N/2      z(t = N/2)      x(f = +1/(2 Delta),
                                    -1/(2 Delta))
     .        ........        ..................
     N-3      z(t = N-3)      x(f = -3/(N Delta))
     N-2      z(t = N-2)      x(f = -2/(N Delta))
     N-1      z(t = N-1)      x(f = -1/(N Delta))

*/
