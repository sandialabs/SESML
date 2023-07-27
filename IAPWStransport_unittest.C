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

#include "IAPWStransport.H"

//
// Test u0 and u1 contribution
//
TEST(IAPWSvisc2008, getViscosityU0U1) {
  
  IAPWSvisc2008 w;

  double v[11];
  v[0]  = w.getViscosityU0U1(  998.,  298.15 );
  v[1]  = w.getViscosityU0U1( 1200.,  298.15 );
  v[2]  = w.getViscosityU0U1( 1000.,  373.15 );
  v[3]  = w.getViscosityU0U1(    1.,  433.15 );
  v[4]  = w.getViscosityU0U1( 1000.,  433.15 );
  v[5]  = w.getViscosityU0U1(    1.,  873.15 );
  v[6]  = w.getViscosityU0U1(  100.,  873.15 );
  v[7]  = w.getViscosityU0U1(  600.,  873.15 );
  v[8]  = w.getViscosityU0U1(    1., 1173.15 );
  v[9]  = w.getViscosityU0U1(  100., 1173.15 );
  v[10] = w.getViscosityU0U1(  400., 1173.15 );
  EXPECT_NEAR(  889.735100 ,  v[0] ,  v[0]*1.e-9);
  EXPECT_NEAR( 1437.649467 ,  v[1] ,  v[1]*1.e-9);
  EXPECT_NEAR(  307.883622 ,  v[2] ,  v[2]*2.e-9);
  EXPECT_NEAR(   14.538324 ,  v[3] ,  v[3]*40.e-9);
  EXPECT_NEAR(  217.685358 ,  v[4] ,  v[4]*2.e-9);
  EXPECT_NEAR(   32.619287 ,  v[5] ,  v[5]*1.e-9);
  EXPECT_NEAR(   35.802262 ,  v[6] ,  v[6]*8.e-9);
  EXPECT_NEAR(   77.430195 ,  v[7] ,  v[7]*3.e-9);
  EXPECT_NEAR(   44.217245 ,  v[8] ,  v[8]*12.e-9);
  EXPECT_NEAR(   47.640433 ,  v[9] ,  v[9]*2.e-9);
  EXPECT_NEAR(   64.154608 , v[10] , v[10]*3.e-9);
}

//
// Test full viscosity
//
#include "IAPWS95.H"
#include <iomanip>
TEST(IAPWSvisc2008, getViscosity) {
  
  IAPWS95 w;
  IAPWSvisc2008 wv;

  std::vector<double> v1,v2;
  std::cout <<std::setprecision(15);
  double r = 122.;
  double t = 647.35;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  double mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 25.520677 , mu , mu*1.e-8 );

  r = 222.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 31.337589 , mu , mu*1.e-8 );

  r = 272.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 36.228143 , mu , mu*1.e-8 );

  r = 322.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 42.961579 , mu , mu*1.e-8 );

  r = 372.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 45.688204 , mu , mu*2.e-8 );

  r = 422.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  EXPECT_NEAR( 49.436256 , mu , mu*1.e-8 );
}

//
// Test l0 and l1 contribution
//
TEST(IAPWSthcon2011, getThermalConductivityL0L1) {
  
  IAPWSthcon2011 w;

  double v;

  v = w.getThermalConductivityL0L1(    0. , 298.15 );
  EXPECT_NEAR( 18.4341883 ,  v ,  v*3.e-9);

  v = w.getThermalConductivityL0L1(  998. , 298.15 );
  EXPECT_NEAR( 607.712868 ,  v ,  v*1.e-9);

  v = w.getThermalConductivityL0L1( 1200. , 298.15 );
  EXPECT_NEAR( 799.038144 ,  v ,  v*1.e-9);

  v = w.getThermalConductivityL0L1(    0. , 873.15 );
  EXPECT_NEAR( 79.1034659 ,  v ,  v*1.e-9);

}

//
// Test full conductivity
//
#include "IAPWS95.H"
#include <iomanip>
TEST(IAPWSthcon2011, getThermalConductivity) {
  
  IAPWS95 w;
  IAPWSvisc2008 wv;
  IAPWSthcon2011 wt;

  std::vector<double> v1,v2;
  std::cout << std::setprecision(15);
  double r = 1.;
  double t = 647.35;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  double mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  double l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 51.9298924 , l , l*1.e-9 );

  r = 122.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 130.922885 , l , l*2.e-9 );

  r = 222.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 367.787459 , l , l*1.e-9 );

  r = 272.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 757.959776 , l , l*1.e-9 );

  r = 322.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 1443.75556 , l , l*2.e-9 );

  r = 372.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 650.319402 , l , l*1.e-9 );

  r = 422.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 448.883487 , l , l*1.e-9 );

  r = 750.;
  w.getvars(r,t,v1);
  w.getvars(r,647.096*1.5,v2);
  mu = wv.getViscosity(r,t,v1[V_KT],v2[V_KT]);
  l = wt.getThermalConductivity(r,t,v1[V_KT],v2[V_KT],v1[V_CV],v1[V_CP],mu);
  EXPECT_NEAR( 600.961346 , l , l*1.e-9 );

}
