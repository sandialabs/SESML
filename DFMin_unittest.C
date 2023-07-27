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

#include "DFMin.H"

//
// Helper class for testing minimization routine
//
class TestDFMin : public DFMin
{
 public:

  TestDFMin() {}
  ~TestDFMin() {}

 protected:
  double f(const double x) {
    return x*x*x+1.;
  }

};

//
// Verification tests for Brent's gradient-free minmizer
//
TEST(DFMin, minimize) {

  TestDFMin d;

  double eps = std::numeric_limits<double>::epsilon();

  double m = d.minimize(0,1,eps);

  EXPECT_NEAR( 0., m, eps*1.e11 );

}
