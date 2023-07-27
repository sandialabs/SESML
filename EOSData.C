/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include "EOSData.H"

//
// Constructor
//
EOSData::EOSData( const unsigned int ninput,
                  const unsigned int noutput )
{

  //
  // initialize storage arrays
  //
  inputs.resize(ninput,0.0);
  outputs.resize(noutput,0.0);

}

//
// Destructor
//
EOSData::~EOSData()
{

}
