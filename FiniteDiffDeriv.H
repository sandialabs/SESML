/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef FINITE_DIFF_DERIV_H
#define FINITE_DIFF_DERIV_H

double FivePointStencilFirstDeriv( const double dx1,
                                   const double dx2,
                                   const double dx3,
                                   const double dx4,
                                   const double y0,
                                   const double y1,
                                   const double y2,
                                   const double y3,
                                   const double y4 );

double FivePointStencilSecondDeriv( const double dx1,
                                    const double dx2,
                                    const double dx3,
                                    const double dx4,
                                    const double y0,
                                    const double y1,
                                    const double y2,
                                    const double y3,
                                    const double y4 );

#include <vector>

double SevenPointStencilDeriv( std::vector<double> & f,
                               double a,
                               int zero,
                               int order );

double FivePointStencilDeriv( std::vector<double> & f,
                              double a,
                              int zero,
                              int order );

#endif
