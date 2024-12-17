#include "DropletArray.h"

int main () {
  DropletArray da;
  int N = 1000;
  double R0 = 0.0013801623456003516;
  double k_e = 0.9/4.;
  double Re_scale = 29059.82905982906;
  double c = 20.;
  double omega = 0.22;
  double rho_l = 1000.;
  double Tsat = 371.53327674824317;
  double kv = 0.01918528617968421;
  double rhof = 614.2155646136666;
  double rhov = 3.4709589577683153;
  double muv = 7.215094092485614e-06;
  double hfg = 316884.8797578494;


  da.initialize(N, R0, k_e, Re_scale, c, omega, rho_l,
                 Tsat, kv, rhof, rhov, muv, hfg, 1000.,
                 10.0);

  da.solve_to_time(0.04);

  // SUNContext ctx;
  // SUNContext_Create(NULL, &ctx);

  // N_Vector u = N_VNew_Serial(2, ctx);
  // NV_Ith_S(u, 1) = 1.;
  // N_VPrint_Serial(u);

  return 0;
}