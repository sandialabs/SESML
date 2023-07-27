/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef DFMIN_H
#define DFMIN_H

class DFMin
{
public:

  //
  // Constructor
  //
  DFMin() {}

  //
  // Destructor
  //
  virtual ~DFMin() {}

  //
  // Minimization routine:
  //
  //   Return x in [ax,bx] that minimizes the function f given the tolerance tol
  //
  double minimize( double ax, double bx, double tol );

protected:

  //
  // function to evaulate
  //
  virtual double f( const double x ) = 0;

};

#endif
