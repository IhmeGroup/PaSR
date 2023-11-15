#include "DropletArray.h"

int main () {
  DropletArray da;
  int N = 1000;
  double R0 = 0.0013801623456003516;
  double k_e = 0.9/4.;
  double Re_scale = 29059.82905982906;
  double c = 20.;
  double omega = 0.22;

  da.initialize(N, R0, k_e, Re_scale, c, omega);

  da.solve_to_time(0.04);

  // SUNContext ctx;
  // SUNContext_Create(NULL, &ctx);

  // N_Vector u = N_VNew_Serial(2, ctx);
  // NV_Ith_S(u, 1) = 1.;
  // N_VPrint_Serial(u);

  return 0;
}