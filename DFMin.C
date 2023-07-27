/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

// DFMin.C

// translation to C++ of GAMS: module DFMIN in NMS
// optimized and modified to allow calling function with data

/*
  C***BEGIN PROLOGUE  DFMIN
  C***CATEGORY NO.    G1A2
  C***KEYWORD(S)  ONE-DIMENSIONAL MINIMIZATION, UNIMODAL FUNCTION
  C***AUTHOR  R. BRENT
  C***DATE WRITTEN    730101  (YYMMDD)
  C***PURPOSE
  C     An approximation to the point where  F  attains a minimum  on
  C     the interval (AX,BX) is determined as the value of the function
  C     DFMIN.
  C
  C PARAMETERS
  C
  C INPUT
  C
  C  AX    (double precision)  left endpoint of initial interval
  C  BX    (double precision) right endpoint of initial interval
  C  F     function subprogram which evaluates F(X)  for any  X
  C        in the interval  (AX,BX)
  C  TOL   (double precision) desired length of the interval of uncertainty 
  C        of the final result ( .ge. 0.0)
  C
  C
  C OUTPUT
  C
  C DFMIN   abcissa approximating the minimizer of F
  C AX     lower bound for minimizer
  C BX     upper bound for minimizer
  C
  C***DESCRIPTION
  C
  C     The method used is a combination of golden section search and
  C     successive parabolic interpolation.  Convergence is never much 
  C     slower than that for a Fibonacci search.  If F has a continuous 
  C     second derivative which is positive at the minimum (which is not
  C     at AX or BX), then convergence is superlinear, and usually of the 
  C     order of about 1.324....
  C
  C     The function F is never evaluated at two points closer together
  C     than EPS*ABS(DFMIN) + (TOL/3), where EPS is approximately the 
  C     square root of the relative machine precision.  If F is a unimodal
  C     function and the computed values of F are always unimodal when
  C     separated by at least EPS*ABS(XSTAR) + (TOL/3), then DFMIN 
  C     approximates the abcissa of the global minimum of F on the 
  C     interval AX,BX with an error less than 3*EPS*ABS(DFMIN) + TOL.  
  C     If F is not unimodal, then DFMIN may approximate a local, but 
  C     perhaps non-global, minimum to the same accuracy.
  C
  C     This function subprogram is a slightly modified version of the
  C     ALGOL 60 procedure LOCALMIN given in Richard Brent, Algorithms for
  C     Minimization Without Derivatives, Prentice-Hall, Inc. (1973).
  C
  C***REFERENCE(S)
  C     Richard Brent, Algorithms for Minimization Without Derivatives,
  C     Prentice-Hall, Inc. (1973).
  C***ROUTINES CALLED   NONE
  C***END PROLOGUE
*/

#include <cmath>
#include <iostream>

#include "DFMin.H"

double DFMin::minimize( double ax, double bx, double tol )
{
  double a,b,c,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w;
  double fu,fv,fw,fx,x;
  double tol3;

  // c is the squared inverse of the golden ratio
  c = 0.5*(3.-sqrt(5.));
  // use a direct computation
  //c = 0.3819660112501051518;

  // eps is approximately the square root of the relative machine
  // precision.
  /*
  eps = 1.;
  for (;;) {
    eps /= 2.;
    tol1 = 1. + eps;
    if (tol1 <= 1.) break;
  }
  */
  // directly use value for IEEE 64bit floating point
  eps = 2.220446049250313108e-16;

  eps = sqrt(eps);
  // skip the square root
  //eps = 1.490116119384765625e-8;

  // initialization
  a = ax;
  b = bx;
  v = a + c*(b - a);
  w = v;
  x = v;
  d = 0.;
  e = 0.;
  fx = f(x);
  fv = fx;
  fw = fx;
  tol3 = tol/3.;

  // main loop starts here
  for (;;) {
    xm = 0.5*(a+b);
    tol1 = eps*fabs(x)+tol3;
    tol2 = 2.*tol1;

    // check stopping criterion
    if (fabs(x-xm) <= tol2-0.5*(b-a)) break;

    // is golden-section necessary
    if (fabs(e) > tol1) {
      // fit parabola
      r = (x - w)*(fx - fv);
      q = (x - v)*(fx - fw);
      p = (x - v)*q - (x - w)*r;
      q = 2.*(q - r);
      if (q > 0.) p = -p;
      else q = -q;
      r = e;
      e = d;

      // is parabola acceptable
      if (fabs(p) >= fabs(0.5*q*r) || p <= q*(a - x) || p >= q*(b - x)) {
        // a golden-section step
        if (x < xm) e = b - x;
        else e = a - x;
        d = c*e;
      }
      else {
        // a parabolic interpolation step
        d = p/q;
        u = x + d;

        // f must not be evaluated too close to ax or bx
        if (u - a < tol2 || b - u < tol2) {
          if (xm >= x) d = tol1;
          else d = -tol1;
        }
      }
    }
    else {
      // a golden-section step
      if (x < xm) e = b - x;
      else e = a - x;
      d = c*e;
    }

    // f must not be evaluated too close to x
    if (fabs(d) < tol1) {
      if (d < 0.) u = x - tol1;
      else u = x + tol1;
    }
    else u = x + d;
    fu = f(u);

    // update  a, b, v, w, and x
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else {
      if (u < x) a = u;
      else b = u;

      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  // end of main loop
  return x;
}
