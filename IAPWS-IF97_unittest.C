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

#include "IAPWS-IF97.H"

//
// Test region 1
//
TEST(IAPWSIF97, Region1) {
  
  Region1 r1;

  std::vector<double> v;
  
  r1.getvars(3000,300,v);

  EXPECT_NEAR( 0.100215168e-2 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.115331273e3  , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.112324818e3  , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.392294792    , v[V_S],  v[V_S]*1.03e-9 );
  EXPECT_NEAR( 0.417301218e1  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.150773921e4  , v[V_CS], v[V_CS]*1.e-9 );

  r1.getvars(80000,300,v);

  EXPECT_NEAR( 0.971180894e-3 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.184142828e3  , v[V_H],  v[V_H]*1.45e-9 );
  EXPECT_NEAR( 0.106448356e3  , v[V_E],  v[V_E]*2.e-9 );
  EXPECT_NEAR( 0.368563852    , v[V_S],  v[V_S]*1.09e-9 );
  EXPECT_NEAR( 0.401008987e1  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.163469054e4  , v[V_CS], v[V_CS]*2.59e-9 );

  r1.getvars(3000,500,v);

  EXPECT_NEAR( 0.120241800e-2 , v[V_V],  v[V_V]*2.81e-9 );
  EXPECT_NEAR( 0.975542239e3  , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.971934985e3  , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.258041912e1  , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.465580682e1  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.124071337e4  , v[V_CS], v[V_CS]*2.5e-9 );

}

//
// Test region 2
//
TEST(IAPWSIF97, Region2) {
  
  Region2 r2;

  std::vector<double> v;
  
  r2.getvars(3.5,300,v);

  EXPECT_NEAR( 0.394913866e2 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.254991145e4 , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.241169160e4 , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.852238967e1 , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.191300162e1 , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.427920172e3 , v[V_CS], v[V_CS]*1.e-9 );

  r2.getvars(3.5,700,v);

  EXPECT_NEAR( 0.923015898e2 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.333568375e4 , v[V_H],  v[V_H]*1.15e-9 );
  EXPECT_NEAR( 0.301262819e4 , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.101749996e2 , v[V_S],  v[V_S]*2.11e-9 );
  EXPECT_NEAR( 0.208141274e1 , v[V_CP], v[V_CP]*1.78e-9 );
  EXPECT_NEAR( 0.644289068e3 , v[V_CS], v[V_CS]*1.e-9 );

  r2.getvars(30000,700,v);

  EXPECT_NEAR( 0.542946619e-2 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.263149474e4  , v[V_H],  v[V_H]*1.85e-9 );
  EXPECT_NEAR( 0.246861076e4  , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.517540298e1  , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.103505092e2  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.480386523e3  , v[V_CS], v[V_CS]*1.e-9 );

}

//
// Test region 3
//
TEST(IAPWSIF97, Region3) {
  
  Region3 r3;

  std::vector<double> v;
  
  r3.getvars(500,650,v);

  EXPECT_NEAR( 0.255837018e5 , v[V_P],  v[V_P]*2.91e-9 );
  EXPECT_NEAR( 0.186343019e4 , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.181226279e4 , v[V_E],  v[V_E]*1.7e-9 );
  EXPECT_NEAR( 0.405427273e1 , v[V_S],  v[V_S]*1.1e-9 );
  EXPECT_NEAR( 0.138935717e2 , v[V_CP], v[V_CP]*13.1e-9 );
  EXPECT_NEAR( 0.502005554e3 , v[V_CS], v[V_CS]*1.e-9 );

  r3.getvars(200,650,v);

  EXPECT_NEAR( 0.222930643e5 , v[V_P],  v[V_P]*1.93e-9 );
  EXPECT_NEAR( 0.237512401e4 , v[V_H],  v[V_H]*1.91e-9 );
  EXPECT_NEAR( 0.226365868e4 , v[V_E],  v[V_E]*1.85e-9 );
  EXPECT_NEAR( 0.485438792e1 , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.446579342e2 , v[V_CP], v[V_CP]*2.21e-9 );
  EXPECT_NEAR( 0.383444594e3 , v[V_CS], v[V_CS]*1.e-9 );

  r3.getvars(500,750,v);

  EXPECT_NEAR( 0.783095639e5 , v[V_P],  v[V_P]*1.e-9 );
  EXPECT_NEAR( 0.225868845e4 , v[V_H],  v[V_H]*2.e-9 );
  EXPECT_NEAR( 0.210206932e4 , v[V_E],  v[V_E]*1.12e-9 );
  EXPECT_NEAR( 0.446971906e1 , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.634165359e1 , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.760696041e3 , v[V_CS], v[V_CS]*1.e-9 );

}

//
// Test region 5
//
TEST(IAPWSIF97, Region5) {
  
  Region5 r5;

  std::vector<double> v;
  
  r5.getvars(500,1500,v);

  EXPECT_NEAR( 0.138455090e1 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.521976855e4 , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.452749310e4 , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.965408875e1 , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.261609445e1 , v[V_CP], v[V_CP]*1.51e-9 );
  EXPECT_NEAR( 0.917068690e3 , v[V_CS], v[V_CS]*1.e-9 );

  r5.getvars(30000,1500,v);

  EXPECT_NEAR( 0.230761299e-1 , v[V_V],  v[V_V]*2.1e-9 );
  EXPECT_NEAR( 0.516723514e4  , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.447495124e4  , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.772970133e1  , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.272724317e1  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.928548002e3  , v[V_CS], v[V_CS]*1.e-9 );

  r5.getvars(30000,2000,v);

  EXPECT_NEAR( 0.311385219e-1 , v[V_V],  v[V_V]*1.e-9 );
  EXPECT_NEAR( 0.657122604e4  , v[V_H],  v[V_H]*1.e-9 );
  EXPECT_NEAR( 0.563707038e4  , v[V_E],  v[V_E]*1.e-9 );
  EXPECT_NEAR( 0.853640523e1  , v[V_S],  v[V_S]*1.e-9 );
  EXPECT_NEAR( 0.288569882e1  , v[V_CP], v[V_CP]*1.e-9 );
  EXPECT_NEAR( 0.106736948e4  , v[V_CS], v[V_CS]*1.15e-9 );

}

//
// Test region 2-3 boundary
//
TEST(IAPWSIF97, Region23Boundary) {
  Region23Boundary r23b;

  double p = r23b.getP(623.15);
  EXPECT_NEAR( 0.165291643e5 , p , p*2.87e-9 );

  double t = r23b.getT(0.165291643e5);
  EXPECT_NEAR( 623.15 , t , t*1.e-9 );

}

//
// Test region 4 boundary
//
TEST(IAPWSIF97, Region4) {
  Region4 r4;

  double p,t;

  p = r4.getP(300.0);
  EXPECT_NEAR( 0.353658941e1 , p , p*1.e-9 );
  p = r4.getP(500.0);
  EXPECT_NEAR( 0.263889776e4 , p , p*1.42e-9 );
  p = r4.getP(600.0);
  EXPECT_NEAR( 0.123443146e5 , p , p*1.76e-9 );

  t = r4.getT(100.0);
  EXPECT_NEAR( 0.372755919e3 , t , t*1.05e-9 );
  t = r4.getT(1000.0);
  EXPECT_NEAR( 0.453035632e3 , t , t*1.e-9 );
  t = r4.getT(10000.0);
  EXPECT_NEAR( 0.584149488e3 , t , t*1.e-9 );

}

#include "FiniteDiffDeriv.H"

void getDerivs( RegionPT * r1,
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
                double & d2fdxdy )
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

void getPTDerivVars( const std::vector<double> & vc,
                     const double dgdp,
                     const double dgdt,
                     const double d2gdp2,
                     const double d2gdt2,
                     const double d2gdpdt,
                     std::vector<double> & v)
{
  //
  // Compute the variables for G(P,T) function given derivatives and G,P,T values in vc
  //
  v = vc;
  v[V_R] = 1./dgdp;
  v[V_V] = dgdp;
  v[V_F] = vc[V_G]-vc[V_P]*dgdp;
  v[V_H] = vc[V_G]-vc[V_T]*dgdt;
  v[V_E] = vc[V_G]-vc[V_T]*dgdt-vc[V_P]*dgdp;
  v[V_S] = -dgdt;
  v[V_CP] = -vc[V_T]*d2gdt2;
  v[V_CV] = -vc[V_T]*(d2gdt2-d2gdpdt*d2gdpdt/d2gdp2);
  v[V_CS] = sqrt(1000.*dgdp*dgdp/(d2gdpdt*d2gdpdt/d2gdt2-d2gdp2));
  v[V_KT] = -dgdp/d2gdp2;
  v[V_G2P] = d2gdp2;
}

void getRTDerivVars( const std::vector<double> & vc,
                     const double dfdr,
                     const double dfdt,
                     const double d2fdr2,
                     const double d2fdt2,
                     const double d2fdrdt,
                     std::vector<double> & v )
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

//
// Test Region1 getvars via finite differences
//
TEST(Region1, getvars) {
  Region1 r1;
  std::vector<double> v,vc;
  double dgdt,d2gdt2,dgdp,d2gdp2,d2gdpdt;

  //
  // Compare finite difference values with exact
  //
  getDerivs(dynamic_cast<RegionPT*>(&r1),NULL,0.01,3000.,300.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*7.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*7.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*7.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*5.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*15.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*15.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*20.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*20.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r1),NULL,0.01,80000.,300.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*4.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*7.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*7.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*4.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*2.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r1),NULL,0.01,3000.,500.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*5.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*5.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*4.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*12.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*243.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*145.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*35.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*35.e-8 );

}

//
// Test Region2 getvars via finite differences
//
TEST(Region2, getvars) {
  Region2 r2;
  std::vector<double> v,vc;
  double dgdt,d2gdt2,dgdp,d2gdp2,d2gdpdt;

  //
  // Compare finite difference values with exact
  //
  getDerivs(dynamic_cast<RegionPT*>(&r2),NULL,0.001,3.5,300.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r2),NULL,0.001,3.5,700.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r2),NULL,0.001,30000.,700.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*16.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*109.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*62.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

}

//
// Test Region3 getvars via finite differences
//
TEST(Region3, getvars) {
  Region3 r3;
  std::vector<double> v,vc;
  double dfdt,d2fdt2,dfdr,d2fdr2,d2fdrdt;

  //
  // Compare finite difference values with exact
  //
  getDerivs(NULL,dynamic_cast<RegionRT*>(&r3),0.001,500.,650.,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
  EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*2.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*3.e-8 );

  getDerivs(NULL,dynamic_cast<RegionRT*>(&r3),0.001,200.,650.,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
  EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*12.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*13.e-8 );

  getDerivs(NULL,dynamic_cast<RegionRT*>(&r3),0.001,500.,750.,V_F,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
  getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);
  EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
  EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );

}

//
// Test Region5 getvars via finite differences
//
TEST(Region5, getvars) {
  Region5 r5;
  std::vector<double> v,vc;
  double dgdt,d2gdt2,dgdp,d2gdp2,d2gdpdt;

  //
  // Compare finite difference values with exact
  //
  getDerivs(dynamic_cast<RegionPT*>(&r5),NULL,0.001,500.,1500.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*2.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*2.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r5),NULL,0.001,30000.,1500.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

  getDerivs(dynamic_cast<RegionPT*>(&r5),NULL,0.001,30000.,2000.,V_G,vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt);
  getPTDerivVars(vc,dgdp,dgdt,d2gdp2,d2gdt2,d2gdpdt,v);
  EXPECT_NEAR( vc[V_R], v[V_R], fabs(v[V_R])*1.e-8 );
  EXPECT_NEAR( vc[V_V], v[V_V], fabs(v[V_V])*1.e-8 );
  EXPECT_NEAR( vc[V_F], v[V_F], fabs(v[V_F])*1.e-8 );
  EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
  EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
  EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
  EXPECT_NEAR( vc[V_CP], v[V_CP], fabs(v[V_CP])*1.e-8 );
  EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*1.e-8 );
  EXPECT_NEAR( vc[V_CS], v[V_CS], fabs(v[V_CS])*1.e-8 );
  EXPECT_NEAR( vc[V_KT], v[V_KT], fabs(v[V_KT])*1.e-8 );
  EXPECT_NEAR( vc[V_G2P], v[V_G2P], fabs(v[V_G2P])*1.e-8 );

}

//
// Test maxwell construction T for regions 1-2
//
TEST(IAPWSIF97, getR12satT) {
  IAPWSIF97 w;
  Region4 r4;
  Region23Boundary r23b;

  double plow = 0.611657*0.999;
  double phigh = r23b.getP(623.15)*1.001;

  for (int i=0;i<1501;i++) {
    double p = plow+(phigh-plow)*i/1500.;
    double tw = w.getR12satT(p);
    double tr4 = r4.getT(p);
    EXPECT_NEAR( tr4, tw, tw*6.e-6 );
  }
}

//
// Test maxwell construction T,rho for region 3
//
TEST(IAPWSIF97, getR3satTmaxwellR) {
  IAPWSIF97 w;
  Region4 r4;
  Region3 r3;
  Region23Boundary r23b;

  double plow = r23b.getP(623.15);
  double phigh = 22064.;
  std::vector<double> v1,v2;

  for (int i=0;i<101;i++) {
    //
    // Test for the temperature close to Region 4
    //
    double p = plow+(phigh-plow)*i/100.;
    double tw = w.getR3satT(p);
    double tr4 = r4.getT(p);
    EXPECT_NEAR( tr4, tw, tw*4.e-6 );

    //
    // test that the densities give the same pressure and G
    //
    double rho1,rho2;
    w.getR3maxwellR(tw,rho1,rho2);
    r3.getvars(rho1,tw,v1);
    EXPECT_NEAR( p, v1[V_P], v1[V_P]*2.e-13 );
    r3.getvars(rho2,tw,v2);
    EXPECT_NEAR( p, v2[V_P], v2[V_P]*1.e-12 );
    EXPECT_NEAR( v1[V_G], v2[V_G], fabs(v2[V_G])*3.e-14 );
  }
}

//
// Test maxwell construction consistency in Region 1,2  for H,E,S lookups
//
TEST(IAPWSIF97, getR12maxwellPXconsistency) {

  IAPWSIF97 w;
  Region23Boundary r23b;
  Region1 r1;
  Region2 r2;

  double plow = 0.611657;
  double phigh = r23b.getP(623.15);

  for (int i=0;i<151;i++) {
    double p = plow+(phigh-plow)*i/150.;
    double t = w.getR12satT(p);
    std::vector<double> v1,v2,vm1,vm2,vm3,vm4;
    r1.getvars(p,t,v1);
    r2.getvars(p,t,v2);

    Vars xtype = V_H;
    double xvar = (v1[V_H]+v2[V_H])/2.;
    w.getR12maxwellPX( p, t, xvar, xtype, vm1 );

    xtype = V_E;
    xvar = vm1[V_E];
    w.getR12maxwellPX( p, t, xvar, xtype, vm2 );

    xtype = V_S;
    xvar = vm1[V_S];
    w.getR12maxwellPX( p, t, xvar, xtype, vm3 );

    xtype = V_V;
    xvar = vm1[V_V];
    w.getR12maxwellPX( p, t, xvar, xtype, vm4 );

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

//
// Test maxwell construction consistency in Region 3 for H,E,S lookups
//
TEST(IAPWSIF97, getR3maxwellPXconsistency) {

  IAPWSIF97 w;
  Region4 r4;
  Region3 r3;

  Region23Boundary r23b;

  double plow = r23b.getP(623.15);
  double phigh = 22064.;

  for (int i=0;i<11;i++) {
    double p = plow+(phigh-plow)*i/10.;
    double t = w.getR3satT(p);
    double rho1,rho2;
    w.getR3maxwellR(t,rho1,rho2);
    std::vector<double> v1,v2,vm1,vm2,vm3,vm4;
    r3.getvars(rho1,t,v1);
    r3.getvars(rho2,t,v2);

    Vars xtype = V_H;
    double xvar = (v1[V_H]+v2[V_H])/2.;
    w.getR3maxwellPX( p, t, xvar, xtype, vm1 );

    xtype = V_E;
    xvar = vm1[V_E];
    w.getR3maxwellPX( p, t, xvar, xtype, vm2 );

    xtype = V_S;
    xvar = vm1[V_S];
    w.getR3maxwellPX( p, t, xvar, xtype, vm3 );

    xtype = V_V;
    xvar = vm1[V_V];
    w.getR3maxwellPX( p, t, xvar, xtype, vm4 );

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

void getMaxwellDerivs( IAPWSIF97 * w,
                       const int region,
                       const double delta,
                       const double p,
                       const double vol,
                       std::vector<double> & vcenter,
                       double & dfdr,
                       double & dfdt,
                       double & d2fdr2,
                       double & d2fdt2,
                       double & d2fdrdt )
{
  std::vector<double> f(5,0.),df(5,0.),fx(5,0.);
  std::vector<double> v;

  //
  // Compute 5x5 stencil for derivatives
  //
  std::vector<double> tv(5);
  for (int i=0;i<5;i++) {
    double pv = p*(1.-delta*(i-2));
    //
    // get the temperature corresponding to this pressure
    //
    if (region == 0) tv[i] = w->getR12satT(pv);
    else if (region == 1) tv[i] = w->getR3satT(pv);
    else throw std::runtime_error("Invalid region specifier in getMaxwellDerivs");

    std::vector<double> rv(5);
    for (int j=0;j<5;j++) {
      double vv = vol*(1.-delta*(j-2));
      if (region == 0) w->getR12maxwellPX( pv, tv[i], vv, V_V, v );
      else if (region == 1) w->getR3maxwellPX( pv, tv[i], vv, V_V, v );
      else throw std::runtime_error("Invalid region specifier in getMaxwellDerivs");

      f[j] = v[V_F];
      rv[j] = v[V_R];
      //
      // save values for t derivatives
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
// Test maxwell construction derivative consistency in Region 1,2 via finite difference
//
TEST(IAPWSIF97, getR12maxwellPXconsistencyD) {

  IAPWSIF97 w;
  Region23Boundary r23b;
  Region1 r1;
  Region2 r2;

  double plow = 0.611657*1.02;
  double phigh = r23b.getP(623.15)*0.99;

  for (int i=0;i<151;i++) {
    double p = plow+(phigh-plow)*i/150.;
    double t = w.getR12satT(p);
    std::vector<double> v1,v2,vm1,vm2,vm3;
    r1.getvars(p,t,v1);
    r2.getvars(p,t,v2);

    double vv = (v1[V_V]+v2[V_V])/2.;
    std::vector<double> vc,v;
    double dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt;

    getMaxwellDerivs(&w,0,0.001,p,vv,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
    //std::cout << " dfdr " << dfdr << " dfdt " << dfdt << " d2fdr2 " << d2fdr2 << " d2fdt2 " << d2fdt2 << " d2fdrdt " << d2fdrdt << std::endl;
    getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);

    EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
    EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*390.e-8 );
    EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
    EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
    EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
    EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*10.e-8 );
    EXPECT_NEAR( vc[V_CS]*sqrt(1000.), v[V_CS], fabs(v[V_CS])*4.5e-8 );
  }
}

//
// Test maxwell construction derivative consistency in Region 3 via finite difference
//
TEST(IAPWSIF97, getR3maxwellPXconsistencyD) {

  IAPWSIF97 w;
  Region4 r4;
  Region3 r3;

  Region23Boundary r23b;

  double plow = r23b.getP(623.15)*1.02;
  double phigh = 22064.*0.99;

  for (int i=0;i<11;i++) {
    double p = plow+(phigh-plow)*i/10.;
    double t = w.getR3satT(p);
    double rho1,rho2;
    w.getR3maxwellR(t,rho1,rho2);
    std::vector<double> v1,v2,vm1,vm2,vm3;
    r3.getvars(rho1,t,v1);
    r3.getvars(rho2,t,v2);

    double vv = (v1[V_V]+v2[V_V])/2.;
    std::vector<double> vc,v;
    double dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt;

    getMaxwellDerivs(&w,1,0.001,p,vv,vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt);
    //std::cout << " dfdr " << dfdr << " dfdt " << dfdt << " d2fdr2 " << d2fdr2 << " d2fdt2 " << d2fdt2 << " d2fdrdt " << d2fdrdt << std::endl;
    getRTDerivVars(vc,dfdr,dfdt,d2fdr2,d2fdt2,d2fdrdt,v);

    EXPECT_NEAR( vc[V_P], v[V_P], fabs(v[V_P])*1.e-8 );
    EXPECT_NEAR( vc[V_G], v[V_G], fabs(v[V_G])*1.e-8 );
    EXPECT_NEAR( vc[V_H], v[V_H], fabs(v[V_H])*1.e-8 );
    EXPECT_NEAR( vc[V_E], v[V_E], fabs(v[V_E])*1.e-8 );
    EXPECT_NEAR( vc[V_S], v[V_S], fabs(v[V_S])*1.e-8 );
    EXPECT_NEAR( vc[V_CV], v[V_CV], fabs(v[V_CV])*300.e-8 );
    EXPECT_NEAR( vc[V_CS]*sqrt(1000.), v[V_CS], fabs(v[V_CS])*150.e-8 );
  }
}

//
// Test temperature look up routine for consistency
//
TEST(IAPWSIF97, getvarsPT) {

  IAPWSIF97 w;
  Region1 r1;
  Region2 r2;
  Region3 r3;
  Region5 r5;

  std::vector<double> vi,v;

  //
  // test invalid phase
  //  
  EXPECT_THROW( w.getvarsPT(3000,300,IAPWSIF97::PH_LGMIX,v) , std::runtime_error );

  //
  // test out of bounds
  //
  EXPECT_THROW( w.getvarsPT(3000000,300,IAPWSIF97::PH_LIQUID,v) , std::runtime_error );
  EXPECT_THROW( w.getvarsPT(0.3,300,IAPWSIF97::PH_LIQUID,v) , std::runtime_error );
  EXPECT_THROW( w.getvarsPT(60000,1100,IAPWSIF97::PH_LIQUID,v) , std::runtime_error );

  r1.getvars(3000,300,vi);
  w.getvarsPT(3000,300,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r1.getvars(80000,300,vi);
  w.getvarsPT(80000,300,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r1.getvars(3000,500,vi);
  w.getvarsPT(3000,500,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r2.getvars(3.5,300,vi);
  w.getvarsPT(3.5,300,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r2.getvars(3.5,700,vi);
  w.getvarsPT(3.5,700,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r2.getvars(30000,700,vi);
  w.getvarsPT(30000,700,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r3.getvars(500,650,vi);
  w.getvarsPT(vi[V_P],650,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_G2P;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*4.e-13 );

  r3.getvars(200,650,vi);
  w.getvarsPT(vi[V_P],650,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_G2P;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*2.e-13 );

  r3.getvars(500,750,vi);
  w.getvarsPT(vi[V_P],750,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_G2P;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*2.e-13 );

  r3.getvars(130,630,vi);
  w.getvarsPT(vi[V_P],630,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_G2P;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*2.e-13 );

  r3.getvars(550,630,vi);
  w.getvarsPT(vi[V_P],630,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_G2P;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*9.e-13 );

  r5.getvars(500,1500,vi);
  w.getvarsPT(500,1500,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r5.getvars(30000,1500,vi);
  w.getvarsPT(30000,1500,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

  r5.getvars(30000,2000,vi);
  w.getvarsPT(30000,2000,IAPWSIF97::PH_LIQUID,v);
  for (int i=0;i<V_M;i++) EXPECT_NEAR( vi[i], v[i], fabs(v[i])*1.e-15 );

}

//
// Test main look up routine for consistency
//
TEST(IAPWSIF97, getvarsPXlgmix) {

  IAPWSIF97 w;
  Region23Boundary r23b;
  Region1 r1;
  Region2 r2;
  Region3 r3;

  double pbound = r23b.getP(623.15);

  double p = (pbound+22064.)/2.;
  double t = w.getR3satT(p);
  double rho1,rho2;
  w.getR3maxwellR(t,rho1,rho2);
  std::vector<double> v1,v2;
  r3.getvars(rho1,t,v1);
  r3.getvars(rho2,t,v2);

  std::vector<double> vm1,vm2,vm3,vm4;

  Vars xtype = V_H;
  double xvar = (v1[V_H]+v2[V_H])/2.;
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm1 );

  xtype = V_E;
  xvar = vm1[V_E];
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm2 );

  xtype = V_S;
  xvar = vm1[V_S];
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm3 );

  //
  // Check that all quantities are consistent
  //
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vm1[i], vm2[i], fabs(vm2[i])*1.e-15);
    EXPECT_NEAR( vm1[i], vm3[i], fabs(vm3[i])*1.e-15);
  }

  p = (pbound+0.611657)/2.;
  t = w.getR12satT(p);
  r1.getvars(p,t,v1);
  r2.getvars(p,t,v2);
  
  xtype = V_H;
  xvar = (v1[V_H]+v2[V_H])/2.;
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm1 );

  xtype = V_E;
  xvar = vm1[V_E];
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm2 );

  xtype = V_S;
  xvar = vm1[V_S];
  w.getvarsPX( p, xvar, IAPWSIF97::PH_LGMIX, xtype, vm3 );

  //
  // Check that all quantities are consistent
  //
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vm1[i], vm2[i], fabs(vm2[i])*1.e-15);
    EXPECT_NEAR( vm1[i], vm3[i], fabs(vm3[i])*1.e-15);
  }

}

TEST(IAPWSIF97, getvarsPXlr1) {
  
  IAPWSIF97 w;
  Region1 r1;

  std::vector<double> vi,v1,v2,v3;
  
  r1.getvars(3000,300,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_LIQUID,V_S,v3);
  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-13 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-13 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.5e-13 );
  }

  r1.getvars(80000,300,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_LIQUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*5.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*5.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*8.e-14 );
  }

  r1.getvars(3000,500,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_LIQUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_FLUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_FLUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_FLUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

}

//
// Test region 2
//
TEST(IAPWSIF97, getvarsPXgr2) {
  
  IAPWSIF97 w;
  Region2 r2;

  std::vector<double> vi,v1,v2,v3;
  
  r2.getvars(3.5,300,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*6.e-13 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*8.e-13 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*5.e-13 );
  }

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_FLUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_FLUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_FLUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*6.e-13 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*8.e-13 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*5.e-13 );
  }

  r2.getvars(3.5,700,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*2.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*2.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*2.e-15 );
  }

  r2.getvars(30000,700,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.5e-13 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

}

//
// Test region 3
//
TEST(IAPWSIF97, getvarsPXr3) {
  
  IAPWSIF97 w;
  Region3 r3;

  std::vector<double> vi,v1,v2,v3;
  
  r3.getvars(500,650,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*12.e-13 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*6.e-13 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*9.e-13 );
  }

  r3.getvars(200,650,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*6.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*3.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-13 );
  }

  r3.getvars(500,750,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*4.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*3.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*3.e-14 );
  }

  r3.getvars(200,640,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*17.e-12 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*20.e-12 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*10.e-12 );
  }

  r3.getvars(450,640,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_LIQUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*6.e-12 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*6.e-12 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*5.e-12 );
  }

  r3.getvars(150,625,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*750.e-14 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*60.e-14 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*45.e-14 );
  }

  r3.getvars(570,625,vi);
  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_LIQUID,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_LIQUID,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_LIQUID,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*6.e-12 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*6.e-12 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*5.e-12 );
  }

}

TEST(IAPWSIF97, getvarsPXgr5) {

  IAPWSIF97 w;
  Region5 r5;

  std::vector<double> vi,v1,v2,v3;
  
  r5.getvars(500,1500,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

  r5.getvars(30000,1500,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

  r5.getvars(30000,2000,vi);

  w.getvarsPX(vi[V_P],vi[V_H],IAPWSIF97::PH_GAS,V_H,v1);
  w.getvarsPX(vi[V_P],vi[V_E],IAPWSIF97::PH_GAS,V_E,v2);
  w.getvarsPX(vi[V_P],vi[V_S],IAPWSIF97::PH_GAS,V_S,v3);

  for (int i=0;i<V_M;i++) {
    EXPECT_NEAR( vi[i], v1[i], fabs(v1[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v2[i], fabs(v2[i])*1.e-15 );
    EXPECT_NEAR( vi[i], v3[i], fabs(v3[i])*1.e-15 );
  }

}

TEST(IAPWSIF97, badinputs) {

  IAPWSIF97 w;

  // pressure out of range for Region 3 maxwell constructions
  EXPECT_THROW( w.getR3satT(0.), std::runtime_error );
  EXPECT_THROW( w.getR3satT(22100.), std::runtime_error );

  // pressure out of range for Region 1-2 maxwell constructions
  EXPECT_THROW( w.getR12satT(0.), std::runtime_error );
  EXPECT_THROW( w.getR12satT(22100.), std::runtime_error );

  // temperature out of range for Region 3 maxwell construction densities
  double rho1,rho2;
  EXPECT_THROW( w.getR3maxwellR(0.,rho1,rho2), std::runtime_error );
  EXPECT_THROW( w.getR3maxwellR(22100.,rho1,rho2), std::runtime_error );

  // bad variable type
  std::vector<double> v;
  EXPECT_THROW( w.getR12maxwellPX(1.,1.,1.,V_R,v), std::runtime_error );
  EXPECT_THROW( w.getR3maxwellPX(1.,1.,1.,V_R,v), std::runtime_error );
  EXPECT_THROW( w.getvarsPX(1.,1.,IAPWSIF97::PH_NONE,V_R,v), std::runtime_error );
  
  // pressure out of range
  EXPECT_THROW( w.getvarsPX(0.,1.,IAPWSIF97::PH_NONE,V_S,v), std::runtime_error );
  EXPECT_THROW( w.getvarsPX(1000000.,1.,IAPWSIF97::PH_NONE,V_S,v), std::runtime_error );

  // invalid phase
  EXPECT_THROW( w.getvarsPX(100.,1.,IAPWSIF97::PH_NONE,V_S,v), std::runtime_error );

}

TEST(IAPWSIF97, getvarsPXbounds) {

  IAPWSIF97 w;
  std::vector<double> v;

  // lgmix: pressure above critical point
  EXPECT_THROW( w.getvarsPX(22222.,1.,IAPWSIF97::PH_LGMIX,V_S,v), std::runtime_error );

  // temperature above region 5
  EXPECT_THROW( w.getvarsPX(10000.,1.e99,IAPWSIF97::PH_FLUID,V_S,v), std::runtime_error );
  EXPECT_THROW( w.getvarsPX(40000.,1.e99,IAPWSIF97::PH_FLUID,V_S,v), std::runtime_error );

  // pressure above region 5
  EXPECT_THROW( w.getvarsPX(60000.,1.e99,IAPWSIF97::PH_FLUID,V_S,v), std::runtime_error );

}
