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

#include "EOSParam.H"

//
// Test that the constructor initializes correctly
//
TEST(EOSParam, Constructor) {
  const EOSParam p(1,2,5);

  ASSERT_EQ( 1U, p.params_double.size() );
  ASSERT_EQ( 2U, p.params_int.size()    );
  ASSERT_EQ( 5U, p.params_string.size() );

  EXPECT_EQ( 0.0, p.params_double[0] );

  EXPECT_EQ( 0, p.params_int[0] );
  EXPECT_EQ( 0, p.params_int[0] );

  for (int i=0;i<5;i++) EXPECT_EQ( "", p.params_string[i] );

}
