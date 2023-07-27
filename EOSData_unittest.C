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

#include "EOSData.H"

//
// Test that constructor initializes correctly
//
TEST(EOSData, Constructor) {
  const EOSData d(3,5);

  ASSERT_EQ( 3U, d.inputs.size()  );
  ASSERT_EQ( 5U, d.outputs.size() );

  for (int i=0;i<3;i++) {
    EXPECT_EQ( 0.0,             d.inputs[i] );
  }

  for (int i=0;i<5;i++) {
    EXPECT_EQ( 0.0,             d.outputs[i] );
  }

}
