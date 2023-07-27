/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include "gtest/gtest.h"

#include "IAPWS95.H"

#include "FiniteDiffDeriv.H"

void getDerivs(RegionPT * r1,
               RegionRT * r2,
               const double delta,
               const double x,
               const double y,
               const Vars vi,
               std::vector<double> & vcenter,
               double & dfdx,
               double & dfdy,
               double & d2fdx2,
               double & d2fdy2,
               double & d2fdxdy)
{
  std::vector<double> f(5,0.),df(5,0.),fx(5,0.);
  std::vector<double> v;

  //
  // Compute 5x5 stencil for derivatives
  //
  for (int i=0;i<5;i++) {
    double xv = x*(1.-delta*(i-2));
    for (int j=0;j<5;j++) {
      double yv = y*(1.-delta*(j-2));
      if (r1 != NULL) r1->getvars(xv,yv,v);
      else if (r2 != NULL) r2->getvars(xv,yv,v);
      else throw std::runtime_error("No valid pointer in getDerivs");
      f[j] = v[vi];
      //
      // save values for x derivatives
      //
      if (j == 2) fx[i] = f[j];

      //
      // save center point
      //
      if (j == 2 && i == 2) vcenter = v;
    }
    //
    // save df/dy|x for cross derivative
    //
    df[i] = FivePointStencilDeriv(f,delta*y,2,1);

    //
    // save y derivatives at center point
    //
    if (i == 2) {
      dfdy = df[i];
      d2fdy2 = FivePointStencilDeriv(f,delta*y,2,2);
    }
  }
  //
  // Compute x derivs
  //
  dfdx = FivePointStencilDeriv(fx,delta*x,2,1);
  d2fdx2 = FivePointStencilDeriv(fx,delta*x,2,2);
  d2fdxdy = FivePointStencilDeriv(df,delta*x,2,1);

}

void getRTDerivVars(const std::vector<double> & vc,
                    const double dfdr,
                    const double dfdt,
                    const double d2fdr2,
                    const double d2fdt2,
                    const double d2fdrdt,
                    std::vector<double> & v)
{
  //
  // Compute the variables for F(R,T) function given derivatives and F,R,T values in vc
  //
  v = vc;
  v[V_P] = vc[V_R]*vc[V_R]*dfdr;
  v[V_G] = vc[V_F]+vc[V_R]*dfdr;
  v[V_H] = vc[V_F]-vc[V_T]*dfdt+vc[V_R]*dfdr;
  v[V_E] = vc[V_F]-vc[V_T]*dfdt;
  v[V_S] = -dfdt;
  v[V_CP] = -vc[V_T]*(d2fdt2-d2fdrdt*d2fdrdt*vc[V_R]/(2.*dfdr+vc[V_R]*d2fdr2));
  v[V_CV] = -vc[V_T]*d2fdt2;
  v[V_CS] = sqrt(1000.*(vc[V_R]*(2.*dfdr+vc[V_R]*d2fdr2)-vc[V_R]*vc[V_R]*d2fdrdt*d2fdrdt/d2fdt2));
  v[V_KT] = vc[V_R]*vc[V_R]*(2.*dfdr+vc[V_R]*d2fdr2);
}

void getRTDerivVars2(const std::vector<double> & vc,
                     const double dfdr,
                     const double dfdt,
                     const double d2fdr2,
                     const double d2fdt2,
                     const double d2fdrdt,
                     std::vector<double> & v)
{
  //
  // Compute the variables for F(R,T) function given derivatives and F,R,T values in vc
  //
  v.resize(6);
  // P
  v[0] = vc[V_R]*vc[V_R]*dfdr;
  // E
  v[1] = vc[V_F]-vc[V_T]*dfdt;
  // dPdR
  v[2] = vc[V_R]*(2.*dfdr+vc[V_R]*d2fdr2);
  // dPdT
  v[3] = vc[V_R]*vc[V_R]*d2fdrdt;
  // dEdR
  v[4] = dfdr-vc[V_T]*d2fdrdt;
  // dEdT
  v[5] = -vc[V_T]*d2fdt2;
}

//
// Test getvars against standard defined points
//
TEST(IAPWS95, getvars)
{
  
  IAPWS95 w;

  std::vector<double> v;

  double relerr = 3.2e-9;

  w.getvars( 0.9965560e3, 300, v );

  EXPECT_NEAR( 0.992418352e2 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 4.13018112    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 1501.51914    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 0.393062643   ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.1005308e4, 300, v );

  EXPECT_NEAR( 0.200022515e5 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 4.06798347    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 1534.92501    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 0.387405401   ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.1188202e4, 300, v );

  EXPECT_NEAR( 0.700004704e6 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 3.46135580    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 2443.57992    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 0.132609616   ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.4350000e0, 500, v );

  EXPECT_NEAR( 0.999679423e2 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 1.50817541    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 548.314253    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 7.94488271    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.4532000e1, 500, v );

  EXPECT_NEAR( 0.999938125e3 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 1.66991025    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 535.739001    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 6.82502725    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.8380250e3, 500, v );

  EXPECT_NEAR( 0.100003858e5 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 3.22106219    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 1271.28441    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 2.56690919    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.1084564e4, 500, v );

  EXPECT_NEAR( 0.700000405e6 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 3.07437693    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 2412.00877    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 2.03237509    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.3580000e3, 647, v );

  EXPECT_NEAR( 0.220384756e5 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 6.18315728    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 252.145078    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 4.32092307    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.2410000e0, 900, v );

  EXPECT_NEAR( 0.100062559e3 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 1.75890657    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 724.027147    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 9.16653194    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.5261500e2, 900, v );

  EXPECT_NEAR( 0.200000690e5 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 1.93510526    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 698.445674    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 6.59070225    ,  v[V_S],  v[V_S]*relerr );

  w.getvars( 0.8707690e3, 900, v );

  EXPECT_NEAR( 0.700000006e6 ,  v[V_P],  v[V_P]*relerr );
  EXPECT_NEAR( 2.66422350    , v[V_CV], v[V_CV]*relerr );
  EXPECT_NEAR( 2019.33608    , v[V_CS], v[V_CS]*relerr );
  EXPECT_NEAR( 4.17223802    ,  v[V_S],  v[V_S]*relerr );

}

//
// Test getvars and getvars2 V_E code via finite differences
//
TEST(IAPWS95, getvarsDerivs)
{

  IAPWS95 w;

  double t,rho;
  std::vector<double> v,vc;
  double dfdt,d2fdt2,dfdr,d2fdr2,d2fdrdt;

  rho = 0.9965560e3; t = 300;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.0007,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*6.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*30.e-8 );

  rho = 0.1005308e4; t = 300;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.1188202e4; t = 300;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.5e-8 );

  rho = 0.4350000e0; t = 500;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.4532000e1; t = 500;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.8380250e3; t = 500;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.1084564e4; t = 500;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.3580000e3; t = 647;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.00009,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*8200.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*7700.e-8 );

  rho = 0.2410000e0; t = 900;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.4e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*40.e-8 );

  rho = 0.5261500e2; t = 900;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

  rho = 0.8707690e3; t = 900;
  getDerivs(NULL,dynamic_cast<RegionRT*>(&w),0.001,rho,t,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );
  getRTDerivVars2(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  w.getvars2( rho, t, V_E, vc );
  for (int i=0;i<6;i++) EXPECT_NEAR( vc[i], v[i], fabs(v[i])*1.e-8 );

}

//
// Test temperature look up routine for consistency
//
TEST(IAPWS95, getvarsPT)
{

  IAPWS95 w;

  std::vector<double> vi,v;
  double rho,t;

  //
  // test invalid phase
  //  
  EXPECT_THROW( w.getvarsPT(1,1,IAPWS95::PH_LGMIX,v) , std::runtime_error );

  //
  // test out of bounds
  //
  EXPECT_THROW( w.getvarsPT(1.e7,1,IAPWS95::PH_LIQUID,v) , std::runtime_error );
  EXPECT_THROW( w.getvarsPT(1.e-1,1,IAPWS95::PH_LIQUID,v) , std::runtime_error );

  rho = 0.9965560e3; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*300000.e-15 );

  rho = 0.1005308e4; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*2000.e-15 );

  rho = 0.1188202e4; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*270.e-15 );

  rho = 0.4350000e0; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  rho = 0.4532000e1; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  rho = 0.8380250e3; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*110.e-15 );

  rho = 0.1084564e4; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*100.e-15 );

  rho = 0.3580000e3; t = 647;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*40000.e-15 );

  rho = 0.2410000e0; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  rho = 0.5261500e2; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  rho = 0.8707690e3; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPT(vi[V_P],t,IAPWS95::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*4.1e-15 );

}


//
// Test main look up routine for consistency
//
TEST(IAPWS95, getvarsPX)
{

  IAPWS95 w;

  std::vector<double> vi,v1,v2,v3;
  double rho,t;

  //
  // test invalid type
  //  
  EXPECT_THROW( w.getvarsPX(1,1,IAPWS95::PH_LGMIX,V_M,vi) , std::runtime_error );

  //
  // test out of bounds pressure
  //
  EXPECT_THROW( w.getvarsPT(1.e7,1,IAPWS95::PH_LIQUID,vi) , std::runtime_error );
  EXPECT_THROW( w.getvarsPT(1.e-1,1,IAPWS95::PH_LIQUID,vi) , std::runtime_error );

  rho = 0.9965560e3; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*30000.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*10000.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*12000.e-14 );
  }

  rho = 0.1005308e4; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*500.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*500.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*550.e-14 );
  }

  rho = 0.1188202e4; t = 300;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*2.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*4.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*7.e-14 );
  }

  rho = 0.4350000e0; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

  rho = 0.4532000e1; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

  rho = 0.8380250e3; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*8.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*5.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*15.e-14 );
  }

  rho = 0.1084564e4; t = 500;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*2.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

  rho = 0.3580000e3; t = 647;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*200.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*330.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1200.e-14 );
  }

  rho = 0.2410000e0; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

  rho = 0.5261500e2; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

  rho = 0.8707690e3; t = 900;
  w.getvars(rho,t,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWS95::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWS95::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWS95::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-14 );
  }

}

//
// Test maxwell construction calculations against standard defined points
//
TEST(IAPWS95, getSatProps) {

  IAPWS95 w;

  double p,t,rho1,rho2;
  std::vector<double> v1,v2;

  double relerr = 4.e-9;

  p = 0.698451167e0;
  EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
  w.getvars( rho1, t, v1 );
  w.getvars( rho2, t, v2 );
  EXPECT_NEAR( 275 , t , t*relerr );
  EXPECT_NEAR( v1[V_P] , v2[V_P] , v2[V_P]*relerr*24. );
  EXPECT_NEAR( 0.698451167e0  , v2[V_P] , v2[V_P]*relerr*24. );
  EXPECT_NEAR( 0.698451167e0  , v1[V_P] , v1[V_P]*relerr );
  EXPECT_NEAR( 0.999887406e3  , v2[V_R] , v2[V_R]*relerr );
  EXPECT_NEAR( 0.550664919e-2 , v1[V_R] , v1[V_R]*relerr );
  EXPECT_NEAR( 0.775972202e1  , v2[V_H] , v2[V_H]*relerr );
  EXPECT_NEAR( 0.250428995e4  , v1[V_H] , v1[V_H]*relerr );
  EXPECT_NEAR( 0.283094670e-1 , v2[V_S] , v2[V_S]*relerr );
  EXPECT_NEAR( 0.910660121e1  , v1[V_S] , v1[V_S]*relerr );

  p = 0.932203564e3;
  EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
  w.getvars( rho1, t, v1 );
  w.getvars( rho2, t, v2 );

  EXPECT_NEAR( 450 , t , t*relerr );
  EXPECT_NEAR( v1[V_P] , v2[V_P] , v2[V_P]*relerr );
  EXPECT_NEAR( 0.932203564e3 , v2[V_P] , v2[V_P]*relerr );
  EXPECT_NEAR( 0.932203564e3 , v1[V_P] , v1[V_P]*relerr );
  EXPECT_NEAR( 0.890341250e3 , v2[V_R] , v2[V_R]*relerr );
  EXPECT_NEAR( 0.481200360e1 , v1[V_R] , v1[V_R]*relerr );
  EXPECT_NEAR( 0.749161585e3 , v2[V_H] , v2[V_H]*relerr );
  EXPECT_NEAR( 0.277441078e4 , v1[V_H] , v1[V_H]*relerr );
  EXPECT_NEAR( 0.210865845e1 , v2[V_S] , v2[V_S]*relerr );
  EXPECT_NEAR( 0.660921221e1 , v1[V_S] , v1[V_S]*relerr );

  p = 0.169082693e5;
  EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
  w.getvars( rho1, t, v1 );
  w.getvars( rho2, t, v2 );

  EXPECT_NEAR( 625 , t , t*relerr );
  EXPECT_NEAR( v1[V_P], v2[V_P], v2[V_P]*relerr );
  EXPECT_NEAR( 0.169082693e5 , v2[V_P] , v2[V_P]*relerr );
  EXPECT_NEAR( 0.169082693e5 , v1[V_P] , v1[V_P]*relerr );
  EXPECT_NEAR( 0.567090385e3 , v2[V_R] , v2[V_R]*relerr );
  EXPECT_NEAR( 0.118290280e3 , v1[V_R] , v1[V_R]*relerr );
  EXPECT_NEAR( 0.168626976e4 , v2[V_H] , v2[V_H]*relerr );
  EXPECT_NEAR( 0.255071625e4 , v1[V_H] , v1[V_H]*relerr );
  EXPECT_NEAR( 0.380194683e1 , v2[V_S] , v2[V_S]*relerr );
  EXPECT_NEAR( 0.518506121e1 , v1[V_S] , v1[V_S]*relerr );

}

//
// Test maxwell construction T, rho solver
//
TEST(IAPWS95, getSatPropsTrho) {

  IAPWS95 w;
  Region4 r4;

  std::vector<double> v1,v2;

  for (int i=0;i<100;i++) {
    //
    // Test for temperature close to Region 4
    //
    double p = 0.611657+(22064.-0.611657)*i/100.;
    double rw1,rw2,tw;
    w.getSatProps(p,rw1,rw2,tw);
    double tr4 = r4.getT(p);
    EXPECT_NEAR( tr4, tw, tw*2.e-5 );

    //
    // Test that the densities give the same P and G
    //
    w.getvars(rw1,tw,v1);
    w.getvars(rw2,tw,v2);
    EXPECT_NEAR( p, v1[V_P], v1[V_P]*1.e-13 );
    if (p < 1000.) {
      EXPECT_NEAR( v1[V_G], v2[V_G], fabs(v2[V_G]+1.)*1.e-9 );
      EXPECT_NEAR( p, v2[V_P], (v2[V_P]+1.)*1.e-8 );
    }
    else {
      EXPECT_NEAR( v1[V_G], v2[V_G], fabs(v2[V_G]+1.)*5.e-14 );
      EXPECT_NEAR( p, v2[V_P], v2[V_P]*5.e-12 );
    }
  }
}

//
// Test maxwell construction consistency for H,E,S lookups
//
TEST(IAPWS95, getMaxwellPXconsistency) {

  IAPWS95 w;

  double plow = 0.611657;
  double phigh = 22064.;

  for (int i=0;i<201;i++) {
    double p = plow+(phigh-plow)*i/200.;
    double r1,r2,t;
    w.getSatProps(p,r1,r2,t);
    std::vector<double> v1,v2,vm1,vm2,vm3,vm4;
    w.getvars(r1,t,v1);
    w.getvars(r2,t,v2);

    Vars xtype = V_H;
    double xvar = (v1[V_H]+v2[V_H])/2.;
    w.getMaxwellPX( p, xvar, xtype, vm1 );

    xtype = V_E;
    xvar = vm1[V_E];
    w.getMaxwellPX( p, xvar, xtype, vm2 );

    xtype = V_S;
    xvar = vm1[V_S];
    w.getMaxwellPX( p, xvar, xtype, vm3 );

    xtype = V_V;
    xvar = vm1[V_V];
    w.getMaxwellPX( p, xvar, xtype, vm4 );

    //
    // Check that all non-derivative quantities are consistent
    //
    EXPECT_NEAR( vm1[V_V], vm2[V_V], fabs(vm2[V_V])*1.e-15);
    EXPECT_NEAR( vm1[V_V], vm3[V_V], fabs(vm3[V_V])*1.e-15);
    EXPECT_NEAR( vm1[V_V], vm4[V_V], fabs(vm4[V_V])*1.e-15);
    EXPECT_NEAR( vm1[V_R], vm2[V_R], fabs(vm2[V_R])*1.e-15);
    EXPECT_NEAR( vm1[V_R], vm3[V_R], fabs(vm3[V_R])*1.e-15);
    EXPECT_NEAR( vm1[V_R], vm4[V_R], fabs(vm4[V_R])*1.e-15);
    EXPECT_NEAR( vm1[V_T], vm2[V_T], fabs(vm2[V_T])*1.e-15);
    EXPECT_NEAR( vm1[V_T], vm3[V_T], fabs(vm3[V_T])*1.e-15);
    EXPECT_NEAR( vm1[V_T], vm4[V_T], fabs(vm4[V_T])*1.e-15);
    EXPECT_NEAR( vm1[V_P], vm2[V_P], fabs(vm2[V_P])*1.e-15);
    EXPECT_NEAR( vm1[V_P], vm3[V_P], fabs(vm3[V_P])*1.e-15);
    EXPECT_NEAR( vm1[V_P], vm4[V_P], fabs(vm4[V_P])*1.e-15);
    EXPECT_NEAR( vm1[V_G], vm2[V_G], fabs(vm2[V_G])*1.e-15);
    EXPECT_NEAR( vm1[V_G], vm3[V_G], fabs(vm3[V_G])*1.e-15);
    EXPECT_NEAR( vm1[V_G], vm4[V_G], fabs(vm4[V_G])*1.e-15);
    EXPECT_NEAR( vm1[V_F], vm2[V_F], fabs(vm2[V_F])*1.e-15);
    EXPECT_NEAR( vm1[V_F], vm3[V_F], fabs(vm3[V_F])*1.e-15);
    EXPECT_NEAR( vm1[V_F], vm4[V_F], fabs(vm4[V_F])*1.e-15);
    EXPECT_NEAR( vm1[V_H], vm2[V_H], fabs(vm2[V_H])*1.e-15);
    EXPECT_NEAR( vm1[V_H], vm3[V_H], fabs(vm3[V_H])*1.e-15);
    EXPECT_NEAR( vm1[V_H], vm4[V_H], fabs(vm4[V_H])*1.e-15);
    EXPECT_NEAR( vm1[V_E], vm2[V_E], fabs(vm2[V_E])*1.e-15);
    EXPECT_NEAR( vm1[V_E], vm3[V_E], fabs(vm3[V_E])*1.e-15);
    EXPECT_NEAR( vm1[V_E], vm4[V_E], fabs(vm4[V_E])*1.e-15);
    EXPECT_NEAR( vm1[V_S], vm2[V_S], fabs(vm2[V_S])*1.e-15);
    EXPECT_NEAR( vm1[V_S], vm3[V_S], fabs(vm3[V_S])*1.e-15);
    EXPECT_NEAR( vm1[V_S], vm4[V_S], fabs(vm4[V_S])*1.e-15);
    EXPECT_NEAR( vm1[V_CV], vm2[V_CV], fabs(vm2[V_CV])*1.e-8);
    EXPECT_NEAR( vm1[V_CV], vm3[V_CV], fabs(vm3[V_CV])*1.e-8);
    EXPECT_NEAR( vm1[V_CV], vm4[V_CV], fabs(vm4[V_CV])*1.e-8);
    EXPECT_NEAR( vm1[V_CS], vm2[V_CS], fabs(vm2[V_CS])*1.e-8);
    EXPECT_NEAR( vm1[V_CS], vm3[V_CS], fabs(vm3[V_CS])*1.e-8);
    EXPECT_NEAR( vm1[V_CS], vm4[V_CS], fabs(vm4[V_CS])*1.e-8);
    // these derivatives are constant
    EXPECT_NEAR( vm1[V_CP], vm2[V_CP], fabs(vm2[V_CP])*1.e-15);
    EXPECT_NEAR( vm1[V_CP], vm3[V_CP], fabs(vm3[V_CP])*1.e-15);
    EXPECT_NEAR( vm1[V_CP], vm4[V_CP], fabs(vm4[V_CP])*1.e-15);
    EXPECT_NEAR( vm1[V_KT], vm2[V_KT], fabs(vm2[V_KT])*1.e-15);
    EXPECT_NEAR( vm1[V_KT], vm3[V_KT], fabs(vm3[V_KT])*1.e-15);
    EXPECT_NEAR( vm1[V_KT], vm4[V_KT], fabs(vm4[V_KT])*1.e-15);
  }
}

void getMaxwellDerivs(IAPWS95 * w,
                      const double delta,
                      const double p,
                      const double vol,
                      std::vector<double> & vcenter,
                      double & dfdr,
                      double & dfdt,
                      double & d2fdr2,
                      double & d2fdt2,
                      double & d2fdrdt)
{
  std::vector<double> f(5,0.),df(5,0.),fx(5,0.),tv(5,0.);
  std::vector<double> v;

  //
  // Compute 5x5 stencil for derivatives
  //
  std::vector<double> pv(5);
  for (int i=0;i<5;i++) {
    pv[i] = p*(1.-delta*(i-2));

    std::vector<double> rv(5);
    for (int j=0;j<5;j++) {
      double vv = vol*(1.-delta*(j-2));
      w->getMaxwellPX( pv[i], vv, V_V, v );

      f[j] = v[V_F];
      rv[j] = v[V_R];
      //
      // save values for t derivatives
      //
      if (j == 2) {
        fx[i] = f[j];
        tv[i] = v[V_T];
      }

      //
      // save center point
      //
      if (j == 2 && i == 2) vcenter = v;
    }
    //
    // save df/dy|x for cross derivative
    //
    //std::cout << " f " << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << f[4] << std::endl;
    //std::cout << " r " << rv[0] << " " << rv[1] << " " << rv[2] << " " << rv[3] << " " << rv[4] << std::endl;
    df[i] = FivePointStencilFirstDeriv( rv[0]-rv[2],rv[1]-rv[2],rv[3]-rv[2],rv[4]-rv[2],
                                        f[2],f[0],f[1],f[3],f[4] );
    //std::cout << " df[" << i << "] " << df[i] << std::endl;
    //
    // save y derivatives at center point
    //
    if (i == 2) {
      dfdr = df[i];
      d2fdr2 = FivePointStencilSecondDeriv( rv[0]-rv[2],rv[1]-rv[2],rv[3]-rv[2],rv[4]-rv[2],
                                            f[2],f[0],f[1],f[3],f[4] );
      //std::cout << " df2[" << i << "] " << d2fdr2 << std::endl;
    }
  }
  //
  // Compute t derivs
  //
  dfdt = FivePointStencilFirstDeriv( tv[0]-tv[2],tv[1]-tv[2],tv[3]-tv[2],tv[4]-tv[2],
                                     fx[2],fx[0],fx[1],fx[3],fx[4] );
  d2fdt2 = FivePointStencilSecondDeriv( tv[0]-tv[2],tv[1]-tv[2],tv[3]-tv[2],tv[4]-tv[2],
                                        fx[2],fx[0],fx[1],fx[3],fx[4] );
  d2fdrdt = FivePointStencilFirstDeriv( tv[0]-tv[2],tv[1]-tv[2],tv[3]-tv[2],tv[4]-tv[2],
                                        df[2],df[0],df[1],df[3],df[4] );
  //std::cout << " t " << tv[0] << " " << tv[1] << " " << tv[2] << " " << tv[3] << " " << tv[4] << std::endl;
}

//
// Test maxwell construction derivative consistency via finite difference
//
TEST(IAPWS95, getMaxwellPXconsistencyD)
{

  IAPWS95 w;

  double plow = 0.611657*1.02;
  double phigh = 22064.*0.99;

  for (int i=0;i<201;i++) {
    double p = plow+(phigh-plow)*i/200.;
    double r1,r2,t;
    w.getSatProps(p,r1,r2,t);
    std::vector<double> v1,v2,vm1,vm2,vm3;
    w.getvars(r1,t,v1);
    w.getvars(r2,t,v2);

    double vv = (v1[V_V]+v2[V_V])/2.;
    std::vector<double> vc,v;
    double dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt;

    getMaxwellDerivs(&w,0.001,p,vv,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
    //std::cout << " dfdr " << dfdr << " dfdt " << dfdt << " d2fdr2 " << d2fdr2 << " d2fdt2 " << d2fdt2 << " d2fdrdt " << d2fdrdt << std::endl;
    getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);

    EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
    EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*800.e-8 );
    EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
    EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
    EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
    EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*150.e-8 );
    EXPECT_NEAR( vc[V_CS]*sqrt(1000.), v[V_CS], fabs(v[V_CS])*100.e-8 );
  }
}

/*
TEST(IAPWS95, asdf) {
  
  IAPWS95 w;
  double p,t,rho1,rho2;
  std::vector<double> v1,v2;

  double relerr = 4.e-9;

  std::cout << std::setprecision(16);
  for (int i=0;i<100001;i++) {
    p = 0.611657+(22064.-0.611657)*i/100000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    w.getvars( rho1, t, v1 );
    w.getvars( rho2, t, v2 );
    std::cout << " p " << p << " dG " << v1[V_G]-v2[V_G] << std::endl;
  }
  for (int i=0;i<100001;i++) {
    p = 22063.+(22064.-22063.)*i/100000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    w.getvars( rho1, t, v1 );
    w.getvars( rho2, t, v2 );
    std::cout << " p " << p << " dG " << v1[V_G]-v2[V_G] << std::endl;
  }
}
*/
/*
TEST(IAPWS95, asdf) {
  
  IAPWS95 w;
  double p,t,rho1,rho2;
  std::vector<double> v1,v2;

  double relerr = 4.e-9;

  std::cout << std::setprecision(16);
  for (int i=0;i<10001;i++) {
    p = 0.611657+(22064.-0.611657)*i/10000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    EXPECT_NO_THROW( w.getMaxwellDerivs(p,v1,v2) );
    std::cout << " p " << p << " t " << t << " r1 " << rho1 << " r2 " << rho2;
    std::cout << " vd1p " << v1[V_P] << " vd2p " << v2[V_P];
    std::cout << " vd1v " << v1[V_V] << " vd2v " << v2[V_V];
    std::cout << " vd1e " << v1[V_E] << " vd2e " << v2[V_E];
    std::cout << std::endl;
  }
  for (int i=0;i<10001;i++) {
    p = 22063.+(22064.-22063.)*i/10000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    EXPECT_NO_THROW( w.getMaxwellDerivs(p,v1,v2) );
    std::cout << " p " << p << " t " << t << " r1 " << rho1 << " r2 " << rho2;
    std::cout << " vd1p " << v1[V_P] << " vd2p " << v2[V_P];
    std::cout << " vd1v " << v1[V_V] << " vd2v " << v2[V_V];
    std::cout << " vd1e " << v1[V_E] << " vd2e " << v2[V_E];
    std::cout << std::endl;
  }
}
*/
/*
TEST(IAPWS95, asdf) {
  
  IAPWS95 w;
  double p,t,rho1,rho2;
  std::vector<double> v1,v2,v3;

  double relerr = 4.e-9;

  std::cout << std::setprecision(16);
  for (int i=999999;i<10001;i++) {
    p = 0.621657+(22064.-0.621657)*i/10000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    w.getvars( rho1, t, v1 );
    w.getMaxwellPX( p,v1[V_E],V_E,v2);
    w.getvars( rho2, t, v1 );
    w.getMaxwellPX( p,v1[V_E],V_E,v3);
    std::cout << " p " << v2[V_P] << " v " << v2[V_V] << " t " << v2[V_T] << " cv " << v2[V_CV] << " cs " << v2[V_CS] << " cv2 " << v3[V_CV] << " cs2 " << v3[V_CS] << std::endl;
  }
  for (int i=0;i<101;i++) {
    p = 22063.+(22064.-22063.)*i/10000.;
    EXPECT_NO_THROW( w.getSatProps(p,rho1,rho2,t) );
    w.getvars( rho1, t, v1 );
    w.getMaxwellPX( p,v1[V_E],V_E,v2);
    w.getvars( rho2, t, v1 );
    w.getMaxwellPX( p,v1[V_E],V_E,v3);
    std::cout << " p " << v2[V_P] << " v " << v2[V_V] << " t " << v2[V_T] << " cv " << v2[V_CV] << " cs " << v2[V_CS] << " cv2 " << v3[V_CV] << " cs2 " << v3[V_CS] << std::endl;
  }
}
*/
/*
//
// Test maxwell construction saturation temperature calculation
//
TEST(IAPWS95, getSatT) {

  IAPWS95 w;
  double t,rho1,rho2,tsat;
  std::vector<double> v;

  for (int i=0;i<1600;i++) {
    t = 273.16+(647.096-273.16)*i/1600.;
    w.getMaxwellR(t,rho1,rho2);
    w.getvars(rho1,t,v);
    tsat = w.getSatT(v[V_P]);
    EXPECT_NEAR( t , tsat , tsat*1.e-13 );
  }
}

*/
