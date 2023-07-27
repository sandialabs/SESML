/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <cmath>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <iomanip>

#include "IAPWStransport.H"

double IAPWSvisc2008::getViscosityU0U1( const double density,
                                        const double temperature )
{
  //
  // reduced variables
  //
  double r = density/rs;
  double ti = ts/temperature;
  //std::cout << "r " << r << " ti "  << ti << std::endl;
  double u0 = 100./(sqrt(ti)*(c0[0].n+ti*(c0[1].n+ti*(c0[2].n+ti*c0[3].n))));

  double u1 = 0.;
  double lt = log(fabs(ti-1.));
  double lr = log(fabs(r-1.));
  for (std::size_t i=0;i<c1.size();i++) {
    if (fabs(r-1.) < 1.e-15) {
      if (fabs(ti-1.) < 1.e-15) {
        if (c1[i].I == 0 && c1[i].J == 0) u1 += c1[i].n;
      }
      else {
        if (c1[i].J == 0) {
          double sgn = 1.;
          if ((c1[i].I % 2) == 1 && ti < 1.) sgn *= -1.;
          u1 += sgn*c1[i].n*exp(c1[i].I*lt);
        }
      }
    }
    else if (fabs(ti-1.) < 1.e-15) {
      if (c1[i].I == 0) {
        double sgn = 1.;
        if ((c1[i].J % 2) == 1 && r < 1.) sgn *= -1.;
        u1 += sgn*c1[i].n*exp(c1[i].J*lr);
      }
    }
    else {
      double sgn = 1.;
      if ((c1[i].I % 2) == 1 && ti < 1.) sgn *= -1.;
      if ((c1[i].J % 2) == 1 && r < 1.) sgn *= -1.;
      u1 += sgn*c1[i].n*exp(c1[i].I*lt+c1[i].J*lr);
    }
    //std::cout << "k " << i << " i " << c1[i].I << " j " << c1[i].J << " H " << c1[i].n << " u " << u1 << std::endl;
  }
  u1 = exp(r*u1);

  //std::cout << "u0 " << u0 << " u1 " << u1 << std::endl;
  return u0*u1;

}



double IAPWSvisc2008::getViscosity( const double density,
                                    const double temperature,
                                    const double kt,
                                    const double kt15 )
{
  //std::cout << "getViscosity: r " << density << " t " << temperature << " kt " << kt << " kt15 " << kt15 << std::endl;
  double dchi = (1./kt-1./kt15*tr*ts/temperature)*ps/g0*density*density/rs/rs;
  double xi = 0.;
  //std::cout << "dchi " << dchi << " ldchi " << log(dchi) << std::endl;
  if (dchi > 0.) xi = xi0*exp(xie*log(dchi));
  double xic = xi*qc;
  double xid = xi*qd;

  double Y=0.;
  if (xi <= 0.3817016416) {
    Y = xic*xid*xid*xid*xid*xid/5.*(1.-xic+xic*xic-765./504.*xid*xid);
  }
  else {
    double phid = acos(1./sqrt(1.+xid*xid));
    double w = sqrt(fabs((xic-1.)/(xic+1.)))*tan(phid/2.);
    double lw;
    if (xic <= 1.) lw = 2.*atan(fabs(w));
    else lw = log((1.+w)/(1.-w));
    Y = sin(3*phid)/12.-sin(2.*phid)/xic/4.+(1.-5./4.*xic*xic)/xic/xic*sin(phid)-1./xic/xic/xic*((1.-3./2.*xic*xic)*phid-exp(3./2.*log(fabs(xic*xic-1.)))*lw);
  }
  //std::cout << "xi " << xi << " mu2 " << exp(xmu*Y) << std::endl;
  return exp(xmu*Y)*getViscosityU0U1(density,temperature);
}



double IAPWSthcon2011::getThermalConductivityL0L1( const double density,
                                                   const double temperature )
{
  //
  // reduced variables
  //
  double r = density/rs;
  double ti = ts/temperature;
  //std::cout << "r " << r << " ti "  << ti << std::endl;
  double l0 = 1./(sqrt(ti)*(c0[0].n+ti*(c0[1].n+ti*(c0[2].n+ti*(c0[3].n+ti*c0[4].n)))));

  double l1 = 0.;
  double lt = log(fabs(ti-1.));
  double lr = log(fabs(r-1.));
  for (std::size_t i=0;i<c1.size();i++) {
    if (fabs(r-1.) < 1.e-15) {
      if (fabs(ti-1.) < 1.e-15) {
        if (c1[i].I == 0 && c1[i].J == 0) l1 += c1[i].n;
      }
      else {
        if (c1[i].J == 0) {
          double sgn = 1.;
          if ((c1[i].I % 2) == 1 && ti < 1.) sgn *= -1.;
          l1 += sgn*c1[i].n*exp(c1[i].I*lt);
        }
      }
    }
    else if (fabs(ti-1.) < 1.e-15) {
      if (c1[i].I == 0) {
        double sgn = 1.;
        if ((c1[i].J % 2) == 1 && r < 1.) sgn *= -1.;
        l1 += sgn*c1[i].n*exp(c1[i].J*lr);
      }
    }
    else {
      double sgn = 1.;
      if ((c1[i].I % 2) == 1 && ti < 1.) sgn *= -1.;
      if ((c1[i].J % 2) == 1 && r < 1.) sgn *= -1.;
      l1 += sgn*c1[i].n*exp(c1[i].I*lt+c1[i].J*lr);
    }
    //std::cout << "k " << i << " i " << c1[i].I << " j " << c1[i].J << " H " << c1[i].n << " l " << l1 << std::endl;
  }
  l1 = exp(r*l1);

  //std::cout << "l0 " << l0 << " l1 " << l1 << std::endl;
  return l0*l1;

}



double IAPWSthcon2011::getThermalConductivity( const double density,
                                               const double temperature,
                                               const double kt,
                                               const double kt15,
                                               const double cv,
                                               const double cp,
                                               const double mu )
{
  double dchi = (1./kt-1./kt15*tr*ts/temperature)*ps/g0*density*density/rs/rs;
  double xi = 0.;
  if (dchi > 0.) xi = xi0*exp(xie*log(dchi));
  double xid = xi*qd;
  double kci = cv/cp;

  double Z = 0.;
  if (xid >= 1.2e-7) {
    Z = 2./M_PI/xid*((1.-kci)*atan(xid)+xid*kci-1.+exp(-1./(1./xid+xid*xid/3./density/density*rs*rs)));
  }

  double l2 = L*density/rs*cp/R*temperature/ts/mu*Z;
  //std::cout << " l2 " << l2 << std::endl;
  return l2+getThermalConductivityL0L1(density,temperature);
}



