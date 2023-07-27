/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>

#include "gtest/gtest.h"

#include "FiniteDiffDeriv.H"

double testfun1(double x)
{
  return x*x*x*x*x*(x-1.);
}

//
// Verification tests for five point stencil derivatives
//
TEST(FiniteDiffDeriv, FivePointStencilDeriv) {

  std::vector<double> f(5);

  double fp = -1.;

  //
  // zero at midpoint
  //
  for (int i=0;i<5;i++) f[i] = testfun1(0.0001*(2-i));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,1) );
  EXPECT_NEAR( fp, 0., 4.e-16 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,2) );
  EXPECT_NEAR( fp, 0., 8.1e-16 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,3) );
  EXPECT_NEAR( fp, 0., 3.e-7 );

  for (int i=0;i<5;i++) f[i] = testfun1(1.*(1.+0.0001*(2-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,1) );
  EXPECT_NEAR( fp, 1., 1.2e-13 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,2) );
  EXPECT_NEAR( fp, 10., 6.e-12 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,2,3) );
  EXPECT_NEAR( fp, 60., 1.6e-6 );

  for (int i=0;i<5;i++) f[i] = testfun1(2.*(1.+0.0001*(2-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,2,1) );
  EXPECT_NEAR( fp, 112., 112.*4.3e-13 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,2,2) );
  EXPECT_NEAR( fp, 320., 320.*3.e-11 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,2,3) );
  EXPECT_NEAR( fp, 720., 720.*4.3e-6 );

  //
  // zero offset by 1
  //
  for (int i=0;i<5;i++) f[i] = testfun1(0.0001*(1-i));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,1) );
  EXPECT_NEAR( fp, 0., 6.1e-16 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,2) );
  EXPECT_NEAR( fp, 0., 1.1e-11 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,3) );
  EXPECT_NEAR( fp, 0., 3.1e-7 );

  for (int i=0;i<5;i++) f[i] = testfun1(1.*(1.+0.0001*(1-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,1) );
  EXPECT_NEAR( fp, 1., 1.2e-13 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,2) );
  EXPECT_NEAR( fp, 10., 4.5e-11 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,1,3) );
  EXPECT_NEAR( fp, 60., 1.6e-6 );

  for (int i=0;i<5;i++) f[i] = testfun1(2.*(1.+0.0001*(1-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,1,1) );
  EXPECT_NEAR( fp, 112., 112.*4.5e-13 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,1,2) );
  EXPECT_NEAR( fp, 320., 5.8e-8 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,1,3) );
  EXPECT_NEAR( fp, 720., 720.*6.e-6 );

  //  
  // zero at edge
  //
  for (int i=0;i<5;i++) f[i] = testfun1(0.0001*(0-i));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,1) );
  EXPECT_NEAR( fp, 0., 2.5e-15 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,2) );
  EXPECT_NEAR( fp, 0., 1.1e-10 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,3) );
  EXPECT_NEAR( fp, 0., 2.2e-6 );

  for (int i=0;i<5;i++) f[i] = testfun1(1.*(1.+0.0001*(0-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,1) );
  EXPECT_NEAR( fp, 1., 1.3e-13 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,2) );
  EXPECT_NEAR( fp, 10., 5.e-10 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,0.0001,0,3) );
  EXPECT_NEAR( fp, 60., 1.1e-5 );

  for (int i=0;i<5;i++) f[i] = testfun1(2.*(1.+0.0001*(0-i)));
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,0,1) );
  EXPECT_NEAR( fp, 112., 1.1e-10 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,0,2) );
  EXPECT_NEAR( fp, 320., 1.8e-6 );
  EXPECT_NO_THROW( fp = FivePointStencilDeriv(f,2.*0.0001,0,3) );
  EXPECT_NEAR( fp, 720., 720.*9.e-6 );

}

//
// Verification tests for five point stencil first derivatives
//
TEST(FiniteDiffDeriv, FivePointStencilXDeriv) {

  std::vector<double> x(5);
  std::vector<double> f(5);

  for (int i=0;i<5;i++) {
    x[i] = 0.0001*(2-i);
    f[i] = testfun1(x[i]);
  }
  
  double fp = -1.;

  EXPECT_NO_THROW( fp = FivePointStencilFirstDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 0., 4.e-16 );
  EXPECT_NO_THROW( fp = FivePointStencilSecondDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 0., 8.e-16 );

  for (int i=0;i<5;i++) {
    x[i] = 1.*(1.+0.00001*(2-i));
    f[i] = testfun1(x[i]);
  }
  EXPECT_NO_THROW( fp = FivePointStencilFirstDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 1., 1.e-15 );
  EXPECT_NO_THROW( fp = FivePointStencilSecondDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 10., 1.e-10 );

  for (int i=0;i<5;i++) {
    x[i] = 2.*(1.+0.0001*(2-i));
    f[i] = testfun1(x[i]);
  }
  EXPECT_NO_THROW( fp = FivePointStencilFirstDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 112., 112.*1.e-12 );
  EXPECT_NO_THROW( fp = FivePointStencilSecondDeriv(x[1]-x[2],x[0]-x[2],x[3]-x[2],x[4]-x[2],f[2],f[1],f[0],f[3],f[4]) );
  EXPECT_NEAR( fp, 320., 320.*1.e-9 );

}

//
// Bad argument testing
//
TEST(FiniteDiffDeriv, FivePointStencilDerivArgs) {
  std::vector<double> f(0);

  // bad zero value
  EXPECT_THROW( FivePointStencilDeriv(f,0.,-1,0), std::logic_error );
  EXPECT_THROW( FivePointStencilDeriv(f,0.,3,0), std::logic_error );

  // bad order value
  EXPECT_THROW( FivePointStencilDeriv(f,0.,1,0), std::logic_error );
  EXPECT_THROW( FivePointStencilDeriv(f,0.,1,4), std::logic_error );

  // bad vector size
  EXPECT_THROW( FivePointStencilDeriv(f,0.,1,1), std::logic_error );

}
