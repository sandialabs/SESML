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

#include "EOSModel.H"

//
// Test that the constructor works
//
TEST(EOSModel, Constructor) {

  EXPECT_NO_THROW( EOSModel m );

}

//
// Test that the getColdCurve method fails on the stub class
//
TEST(EOSModel, getColdCurve) {
  EOSModel m;
  EOSData d(0,0);

  EXPECT_THROW( m.getColdCurve(d), std::logic_error );

}

//
// Test that the getThermalEOS_F method fails on the stub class
//
TEST(EOSModel, getThermalEOS_F) {
  EOSModel m;
  EOSData d(0,0);

  EXPECT_THROW( m.getThermalEOS_F(d), std::logic_error );

}

//
// Test that the getThermalEOS_S method fails on the stub class
//
TEST(EOSModel, getThermalEOS_S) {
  EOSModel m;
  EOSData d(0,0);

  EXPECT_THROW( m.getThermalEOS_S(d), std::logic_error );

}

//
// Test that the getThermalEOS_G method fails on the stub class
//
TEST(EOSModel, getThermalEOS_G) {
  EOSModel m;
  EOSData d(0,0);

  EXPECT_THROW( m.getThermalEOS_G(d), std::logic_error );

}

//
// Helper class for testing protected members
//
class EOSModelTest : public EOSModel
{
 public:
  EOSModelTest(){}
  ~EOSModelTest(){}

  const EOSParam * test_getParams( const Parameters & params,
                                   const std::string name ) {
    return getParams(params,name);
  }

};

//
// Test parameter searching
//
TEST(EOSModel, getParams) {
  EOSModelTest t;

  Parameters p;
  
  p["test1"] = new EOSParam(1,0,0);
  p["test2"] = new EOSParam(0,1,0);

  EXPECT_THROW( t.test_getParams(p,"test"), std::runtime_error );
  
  EXPECT_NO_THROW( t.test_getParams(p,"test1") );
  EXPECT_EQ( p["test1"], t.test_getParams(p,"test1") );

}

