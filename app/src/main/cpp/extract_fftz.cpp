//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: extract_fftz.cpp
//
// MATLAB Coder version            : 3.3
// C/C++ source code generated on  : 23-May-2019 11:54:00
//

// Include Files
#include "rt_nonfinite.h"
#include "extract_fftz.h"

// Function Declarations
static void b_abs(const creal_T x[16000], double y[16000]);
static void b_r2br_r2dit_trig(const creal_T x[32768], const double costab[16385],
  const double sintab[16385], creal_T y[32768]);
static void bluestein_setup(creal_T wwc[31999]);
static void fft(const double x[16000], creal_T y[16000]);
static void filter(double b[4], double a[4], const double x[16018], const double
                   zi[3], double y[16018]);
static void filtfilt(const double x_in[16000], double y_out[16000]);
static void flipud(double x[16018]);
static void generate_twiddle_tables(double costab[16385], double sintab[16385],
  double sintabinv[16385]);
static void r2br_r2dit_trig(const creal_T x[31999], const double costab[16385],
  const double sintab[16385], creal_T y[32768]);
static void r2br_r2dit_trig_impl(const creal_T x[16000], int xoffInit, const
  double costab[16385], const double sintab[16385], creal_T y[32768]);
static double rt_hypotd_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : const creal_T x[16000]
//                double y[16000]
// Return Type  : void
//
static void b_abs(const creal_T x[16000], double y[16000])
{
  int k;
  for (k = 0; k < 16000; k++) {
    y[k] = rt_hypotd_snf(x[k].re, x[k].im);
  }
}

//
// Arguments    : const creal_T x[32768]
//                const double costab[16385]
//                const double sintab[16385]
//                creal_T y[32768]
// Return Type  : void
//
static void b_r2br_r2dit_trig(const creal_T x[32768], const double costab[16385],
  const double sintab[16385], creal_T y[32768])
{
  int ix;
  int ju;
  int iy;
  int i;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 32767; i++) {
    y[iy] = x[ix];
    iy = 32768;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 32767; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 8192;
  iheight = 32765;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 16384; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix += ix;
    iheight -= iy;
  }

  for (iy = 0; iy < 32768; iy++) {
    y[iy].re *= 3.0517578125E-5;
    y[iy].im *= 3.0517578125E-5;
  }
}

//
// Arguments    : creal_T wwc[31999]
// Return Type  : void
//
static void bluestein_setup(creal_T wwc[31999])
{
  int idx;
  int rt;
  int k;
  int y;
  double nt_im;
  double nt_re;
  idx = 15998;
  rt = 0;
  wwc[15999].re = 1.0;
  wwc[15999].im = 0.0;
  for (k = 0; k < 15999; k++) {
    y = ((k + 1) << 1) - 1;
    if (32000 - rt <= y) {
      rt = (y + rt) - 32000;
    } else {
      rt += y;
    }

    nt_im = -3.1415926535897931 * (double)rt / 16000.0;
    if (nt_im == 0.0) {
      nt_re = 1.0;
      nt_im = 0.0;
    } else {
      nt_re = std::cos(nt_im);
      nt_im = std::sin(nt_im);
    }

    wwc[idx].re = nt_re;
    wwc[idx].im = -nt_im;
    idx--;
  }

  idx = 0;
  for (k = 15998; k >= 0; k += -1) {
    wwc[k + 16000] = wwc[idx];
    idx++;
  }
}

//
// Arguments    : const double x[16000]
//                creal_T y[16000]
// Return Type  : void
//
static void fft(const double x[16000], creal_T y[16000])
{
  static double costab[16385];
  static double sintab[16385];
  static double sintabinv[16385];
  static creal_T wwc[31999];
  int xidx;
  int k;
  static creal_T fy[32768];
  static creal_T fv[32768];
  double fy_re;
  generate_twiddle_tables(costab, sintab, sintabinv);
  bluestein_setup(wwc);
  xidx = 0;
  for (k = 0; k < 16000; k++) {
    y[k].re = wwc[k + 15999].re * x[xidx];
    y[k].im = wwc[k + 15999].im * -x[xidx];
    xidx++;
  }

  r2br_r2dit_trig_impl(y, 0, costab, sintab, fy);
  r2br_r2dit_trig(wwc, costab, sintab, fv);
  for (xidx = 0; xidx < 32768; xidx++) {
    fy_re = fy[xidx].re;
    fy[xidx].re = fy[xidx].re * fv[xidx].re - fy[xidx].im * fv[xidx].im;
    fy[xidx].im = fy_re * fv[xidx].im + fy[xidx].im * fv[xidx].re;
  }

  b_r2br_r2dit_trig(fy, costab, sintabinv, fv);
  xidx = 0;
  for (k = 0; k < 16000; k++) {
    y[xidx].re = wwc[k + 15999].re * fv[k + 15999].re + wwc[k + 15999].im * fv[k
      + 15999].im;
    y[xidx].im = wwc[k + 15999].re * fv[k + 15999].im - wwc[k + 15999].im * fv[k
      + 15999].re;
    xidx++;
  }
}

//
// Arguments    : double b[4]
//                double a[4]
//                const double x[16018]
//                const double zi[3]
//                double y[16018]
// Return Type  : void
//
static void filter(double b[4], double a[4], const double x[16018], const double
                   zi[3], double y[16018])
{
  double a1;
  int k;
  int naxpy;
  int j;
  a1 = a[0];
  if ((!rtIsInf(a[0])) && (!rtIsNaN(a[0])) && (!(a[0] == 0.0)) && (a[0] != 1.0))
  {
    for (k = 0; k < 4; k++) {
      b[k] /= a1;
    }

    for (k = 0; k < 3; k++) {
      a[k + 1] /= a1;
    }

    a[0] = 1.0;
  }

  for (k = 0; k < 3; k++) {
    y[k] = zi[k];
  }

  memset(&y[3], 0, 16015U * sizeof(double));
  for (k = 0; k < 16018; k++) {
    naxpy = 16018 - k;
    if (!(naxpy < 4)) {
      naxpy = 4;
    }

    for (j = 0; j + 1 <= naxpy; j++) {
      y[k + j] += x[k] * b[j];
    }

    naxpy = 16017 - k;
    if (!(naxpy < 3)) {
      naxpy = 3;
    }

    a1 = -y[k];
    for (j = 1; j <= naxpy; j++) {
      y[k + j] += a1 * a[j];
    }
  }
}

//
// Arguments    : const double x_in[16000]
//                double y_out[16000]
// Return Type  : void
//
static void filtfilt(const double x_in[16000], double y_out[16000])
{
  double d0;
  double d1;
  int i;
  static double y[16018];
  double dv0[4];
  static const double dv1[4] = { 0.996863335697075, -2.990590007091225,
    2.990590007091225, -0.996863335697075 };

  double dv2[4];
  static const double dv3[4] = { 1.0, -2.9937168172766531, 2.9874533582428482,
    -0.993736510057099 };

  double b_y[16018];
  double a[3];
  static const double b_a[3] = { -0.996863335446352, 1.9937266708942794,
    -0.99686333544792249 };

  d0 = 2.0 * x_in[0];
  d1 = 2.0 * x_in[15999];
  for (i = 0; i < 9; i++) {
    y[i] = d0 - x_in[9 - i];
  }

  memcpy(&y[9], &x_in[0], 16000U * sizeof(double));
  for (i = 0; i < 9; i++) {
    y[i + 16009] = d1 - x_in[15998 - i];
  }

  for (i = 0; i < 4; i++) {
    dv0[i] = dv1[i];
    dv2[i] = dv3[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 16018U * sizeof(double));
  filter(dv0, dv2, b_y, a, y);
  flipud(y);
  for (i = 0; i < 4; i++) {
    dv0[i] = dv1[i];
    dv2[i] = dv3[i];
  }

  for (i = 0; i < 3; i++) {
    a[i] = b_a[i] * y[0];
  }

  memcpy(&b_y[0], &y[0], 16018U * sizeof(double));
  filter(dv0, dv2, b_y, a, y);
  flipud(y);
  memcpy(&y_out[0], &y[9], 16000U * sizeof(double));
}

//
// Arguments    : double x[16018]
// Return Type  : void
//
static void flipud(double x[16018])
{
  int i;
  double xtmp;
  for (i = 0; i < 8009; i++) {
    xtmp = x[i];
    x[i] = x[16017 - i];
    x[16017 - i] = xtmp;
  }
}

//
// Arguments    : double costab[16385]
//                double sintab[16385]
//                double sintabinv[16385]
// Return Type  : void
//
static void generate_twiddle_tables(double costab[16385], double sintab[16385],
  double sintabinv[16385])
{
  double costab1q[8193];
  int k;
  costab1q[0] = 1.0;
  for (k = 0; k < 4096; k++) {
    costab1q[k + 1] = std::cos(0.00019174759848570515 * ((double)k + 1.0));
  }

  for (k = 0; k < 4095; k++) {
    costab1q[k + 4097] = std::sin(0.00019174759848570515 * (8192.0 - ((double)k
      + 4097.0)));
  }

  costab1q[8192] = 0.0;
  costab[0] = 1.0;
  sintab[0] = 0.0;
  for (k = 0; k < 8192; k++) {
    sintabinv[k + 1] = costab1q[8191 - k];
    sintabinv[k + 8193] = costab1q[k + 1];
    costab[k + 1] = costab1q[k + 1];
    sintab[k + 1] = -costab1q[8191 - k];
    costab[k + 8193] = -costab1q[8191 - k];
    sintab[k + 8193] = -costab1q[k + 1];
  }
}

//
// Arguments    : const creal_T x[31999]
//                const double costab[16385]
//                const double sintab[16385]
//                creal_T y[32768]
// Return Type  : void
//
static void r2br_r2dit_trig(const creal_T x[31999], const double costab[16385],
  const double sintab[16385], creal_T y[32768])
{
  int i;
  int ix;
  int ju;
  int iy;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  for (i = 0; i < 32768; i++) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }

  ix = 0;
  ju = 0;
  iy = 0;
  for (i = 0; i < 31998; i++) {
    y[iy] = x[ix];
    iy = 32768;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 32767; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 8192;
  iheight = 32765;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 16384; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix += ix;
    iheight -= iy;
  }
}

//
// Arguments    : const creal_T x[16000]
//                int xoffInit
//                const double costab[16385]
//                const double sintab[16385]
//                creal_T y[32768]
// Return Type  : void
//
static void r2br_r2dit_trig_impl(const creal_T x[16000], int xoffInit, const
  double costab[16385], const double sintab[16385], creal_T y[32768])
{
  int i;
  int ix;
  int ju;
  int iy;
  boolean_T tst;
  double temp_re;
  double temp_im;
  int iheight;
  int istart;
  int j;
  double twid_re;
  double twid_im;
  int ihi;
  for (i = 0; i < 32768; i++) {
    y[i].re = 0.0;
    y[i].im = 0.0;
  }

  ix = xoffInit;
  ju = 0;
  iy = 0;
  for (i = 0; i < 15999; i++) {
    y[iy] = x[ix];
    iy = 32768;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }

    iy = ju;
    ix++;
  }

  y[iy] = x[ix];
  for (i = 0; i <= 32767; i += 2) {
    temp_re = y[i + 1].re;
    temp_im = y[i + 1].im;
    y[i + 1].re = y[i].re - y[i + 1].re;
    y[i + 1].im = y[i].im - y[i + 1].im;
    y[i].re += temp_re;
    y[i].im += temp_im;
  }

  iy = 2;
  ix = 4;
  ju = 8192;
  iheight = 32765;
  while (ju > 0) {
    for (i = 0; i < iheight; i += ix) {
      temp_re = y[i + iy].re;
      temp_im = y[i + iy].im;
      y[i + iy].re = y[i].re - temp_re;
      y[i + iy].im = y[i].im - temp_im;
      y[i].re += temp_re;
      y[i].im += temp_im;
    }

    istart = 1;
    for (j = ju; j < 16384; j += ju) {
      twid_re = costab[j];
      twid_im = sintab[j];
      i = istart;
      ihi = istart + iheight;
      while (i < ihi) {
        temp_re = twid_re * y[i + iy].re - twid_im * y[i + iy].im;
        temp_im = twid_re * y[i + iy].im + twid_im * y[i + iy].re;
        y[i + iy].re = y[i].re - temp_re;
        y[i + iy].im = y[i].im - temp_im;
        y[i].re += temp_re;
        y[i].im += temp_im;
        i += ix;
      }

      istart++;
    }

    ju /= 2;
    iy = ix;
    ix += ix;
    iheight -= iy;
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

//
// EXTRACT_FFTZ Extracts FFT from signal at sampling rate 4000Hz, towards
// detection of a peak at 1kHz
//  INPUTS:
//  X: raw signal voltage, (4*Fs, 1);
//  OUTPUTS:
//  fsel: selected frequencies
//  Y: FFT(including frequencies f=900:1100)
//  Fs = 4000;
//  [bh, ah] = butter(3, 2*2/Fs, 'high');
// Arguments    : const double X[16000]
//                double Y[800]
// Return Type  : void
//
void extract_fftz(const double X[16000], double Y[800])
{
  static double Xfilt[16000];
  static creal_T Yfft[16000];
  int i;
  static creal_T b_Yfft[16000];
  double P1[8001];

  //  f = Fs*(0:(L/2))/L; fsel = f(3601:4400);
  filtfilt(X, Xfilt);
  fft(Xfilt, Yfft);
  for (i = 0; i < 16000; i++) {
    if (Yfft[i].im == 0.0) {
      b_Yfft[i].re = Yfft[i].re / 16000.0;
      b_Yfft[i].im = 0.0;
    } else if (Yfft[i].re == 0.0) {
      b_Yfft[i].re = 0.0;
      b_Yfft[i].im = Yfft[i].im / 16000.0;
    } else {
      b_Yfft[i].re = Yfft[i].re / 16000.0;
      b_Yfft[i].im = Yfft[i].im / 16000.0;
    }
  }

  b_abs(b_Yfft, Xfilt);
  memcpy(&P1[0], &Xfilt[0], 8001U * sizeof(double));
  for (i = 0; i < 7999; i++) {
    P1[1 + i] = 2.0 * Xfilt[1 + i];
  }

  memcpy(&Y[0], &P1[3600], 800U * sizeof(double));
}

//
// Arguments    : void
// Return Type  : void
//
void extract_fftz_initialize()
{
  rt_InitInfAndNaN(8U);
}

//
// Arguments    : void
// Return Type  : void
//
void extract_fftz_terminate()
{
  // (no terminate code required)
}

//
// File trailer for extract_fftz.cpp
//
// [EOF]
//
