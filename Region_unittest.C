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

//
// use the IAPWS-IF97 regions as test instantiations of the Region classes
//
#include "IAPWS-IF97.H"
#include "Region.H"

//
// Test error handling
//
TEST(RegionPT, getTatPX) {
  Region1 r1;

  // lower and upper bound checks
  EXPECT_THROW( r1.getTatPX( 10000., 0., V_H, 273.15, 300.00 ), std::runtime_error );
  EXPECT_THROW( r1.getTatPX( 10000., 1.e99, V_H, 273.15, 300.00 ), std::runtime_error );

}

TEST(RegionRT, getRTatPX) {
  Region3 r3;
  double rho,temp;

  // invalid variable check
  EXPECT_THROW( r3.getRTatPX( 50000., 0., V_R, 1.0, 1.0, rho, temp ), std::runtime_error );

}

TEST(RegionRT, getRatPT) {
  Region3 r3;

  // side check
  EXPECT_THROW( r3.getRatPT( 22000., 630., 600., 650.00, 0 ), std::runtime_error );

  // lower and upper bound checks
  EXPECT_THROW( r3.getRatPT( 2000., 650., 300., 950.00, 0 ), std::runtime_error );
  EXPECT_THROW( r3.getRatPT( 90000., 650., 50., 60.00, 0 ), std::runtime_error );

}
