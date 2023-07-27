/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <iomanip>

#include "Region.H"
#include "Vars.H"

double RegionPT::getTatPX( const double pressure,
                           const double xvar,
                           const int xtype,
                           const double minT,
                           const double maxT )
{
  //
  // evaluation vector
  //
  std::vector<double> v;

  //
  // check lower bound
  //
  getvars(pressure, minT, v);
  if ( v[xtype] > xvar ) {
    throw std::runtime_error("RegionPT::getTatPX: xvar too large at low temperature bound");
  }

  double ls = v[xtype];
  double lt = minT;
  double tmin = minT;

  //
  // check upper bound
  //
  getvars(pressure, maxT, v);
  if ( v[xtype] < xvar ) {
    throw std::runtime_error("RegionPT::getTatPX: xvar too small at high temperature bound");
  }

  double lls = v[xtype];
  double llt = maxT;
  double tmax = maxT;

  //
  // secant/bisection search for T
  //
  while (fabs(tmax-tmin) > 1.e-14*fabs(tmax+tmin)/2.) {
    double t = lt+(xvar-ls)*(lt-llt)/(ls-lls);
    if (t >= tmax || t <= tmin || ls-lls == 0.) {
      t = (tmax+tmin)/2.;
    }
    getvars(pressure, t, v);
    
    if (v[xtype] < xvar) tmin = t;
    else tmax = t;
    llt = lt;
    lt = t;
    lls = ls;
    ls = v[xtype];
  }
  
  return tmax;
}

void RegionRT::getRTatPX( const double pressure,
                          const double xvar,
                          const Vars xtype,
                          const double density0,
                          const double temperature0,
                          double & density,
                          double & temperature )
{
  int debug = 0;

  if (debug > 0) {
    std::cout << std::setprecision(18);
    std::cout << "RegionRT::getRTatPX: p " << pressure << " x " << xvar << " type " << xtype << " r0 " << density0 << " t0 " << temperature0 << std::endl;
    std::cout << std::setprecision(6);
  }

  //
  // check xtype validity
  //
  if (xtype != V_H && xtype != V_E && xtype != V_S) {
    throw std::runtime_error("RegionRT::getRTatPX: invalid variable type");
  }

  // starting point
  double r = density0;
  double t = temperature0;
  double alpha;
  double ra,ta;

  std::vector<double> v;
  double dp,dx;
  double norm;
  double dr(1.0),dt(1.0);

  int i=0;
  int j=0;
  getvars2(r,t,xtype,v);
  while(i<60) {
    // function values
    dp = v[0]-pressure;
    dx = v[1]-xvar;
    if (fabs(dr) < 1.e-14*r && fabs(dt) < 1.e-14*t) break;
    // jacobian
    norm = v[2]*v[5]-v[3]*v[4];
    if (fabs(norm) < 1.e-14) throw std::runtime_error("RegionRT::getRTatPX: Jacobian is singular");
    dr = (-v[5]*dp+v[3]*dx)/norm;
    dt = (v[4]*dp-v[2]*dx)/norm;
    alpha = 1.0;
    if (alpha*fabs(dr) > 100.) alpha = 100./fabs(dr);
    if (r+alpha*dr < 0.) alpha = r*0.9/fabs(dr);
    j=0;
    while(j<10) {
      ra = r+alpha*dr;
      ta = t+alpha*dt;
      double rho1 = rbounds[0];
      double rho2 = rbounds[1];
      double rhos1 = rbounds[1];
      double rhos2 = rbounds[0];
      double tmin = tbounds[0];
      double tmax = tbounds[1];
      if (ta < ts) {
        rhos2 = getSpinodalDensityBound( ta, true );
        rhos1 = getSpinodalDensityBound( ta, false );
      }
      v[2] = -9.e99;
      if (ra > rho1 && ra < rho2 && (ra < rhos1 || ra > rhos2)
          && ta > tmin && ta < tmax) getvars2(ra,ta,xtype,v);
      if (debug > 0) std::cout << "i " << i << " j " << j << " rhos1 " << rhos1 << " rhos2 " << rhos2 << " alpha " << alpha << " ra " << ra << " ta " << ta << " kt " << r*v[2] << std::endl;
      if (v[2] < 0.0 || v[0] < 0.0) alpha *= 0.1;
      else break;
      j++;
    }
    r = ra;
    t = ta;
    if (debug > 0) std::cout << "i " << i << " j " << j << " ptgt " << pressure << " p " << v[0] << " xtgt " << xvar << " x " << v[1] << " r " << r << " dr " << dr << " t " << t << " dt " << dt << " kt " << r*v[2] << std::endl;
    if (j >= 10) break;
    i++;
  }
  if (j >= 10) {
    if (debug > 0) std::cout << "dr " << fabs(dr) << " drt " << 1.e-14*r << " dt " << fabs(dt) << " dtt " << 1.e-14*t << " dp " << dp/pressure << " dx " << dx/xvar << std::endl;
    if (fabs(dp/pressure) > 1.e-5 || fabs(dx/xvar) > 1.e-5) {
      if (debug > 0) {
        std::cout << std::setprecision(18);
        std::cout << "RegionRT::getRTatPX: p " << pressure << " x " << xvar << " type " << xtype << " r0 " << density0 << " t0 " << temperature0 << std::endl;
        std::cout << std::setprecision(6);
      }
      throw std::runtime_error("RegionRT::getRTatPX: unable to take valid step");
    }
  }
  if (i >= 60) {
    if (debug > 0) std::cout << "dr " << fabs(dr) << " drt " << 1.e-14*r << " dt " << fabs(dt) << " dtt " << 1.e-14*t << " dp " << dp/pressure << " de " << dx/xvar << std::endl;
    if (fabs(dp/pressure) > 1.e-5 || fabs(dx/xvar) > 1.e-5) {
      if (debug > 0) {
        std::cout << std::setprecision(18);
        std::cout << "RegionRT::getRTatPX: p " << pressure << " x " << xvar << " type " << xtype << " r0 " << density0 << " t0 " << temperature0 << std::endl;
        std::cout << std::setprecision(6);
      }
      throw std::runtime_error("RegionRT::getRTatPX: too many iterations without convergence");
    }
  }
  
  density = r;
  temperature = t;

}

double RegionRT::getRatPT( const double pressure,
                           const double temperature,
                           const double minR,
                           const double maxR,
                           const int side )
{
  //
  // evaluation vector
  //
  std::vector<double> v;

  //
  // check for side
  //
  if (temperature < ts && side == 0 && pressure < pcrit) throw std::runtime_error("RegionRT::getRatPT: side must be specified at temperatures below ts");

  //
  // min and max densities tailored to the IAPWS-IF97 region 3
  // densities between 40-800 kg/m^3 bounds the pressure region
  // between 16.5-100 MPa from 600-900 K
  //
  double minr = minR;
  if (minr < rbounds[0]) minr = rbounds[0];
  double maxr = maxR;
  if (maxr > rbounds[1]) maxr = rbounds[1];

  //
  // check lower bound
  //
  getvars(minr, temperature, v);
  if ( pressure < v[V_P] && (side < 1 || v[V_KT] >= 0.) && minr < rs) {
    std::cout << "minr " << minr << " T " << temperature << " p " << v[V_P] << " pressure " << pressure << " kt " << v[V_KT] << std::endl;
    throw std::runtime_error("RegionRT::getRatPT: pressure too large at low density bound");
  }

  double lr = minr;
  double lp = v[V_P];
  double rmin = minr;

  //
  // check upper bound
  //
  getvars(maxr, temperature, v);

  if ( pressure > v[V_P]*(1.+1.e-9) && (side > -1 || v[V_KT] >= 0.)) {
    std::cout << "maxr " << maxr << " T " << temperature << " p " << v[V_P] << " pressure " << pressure << " kt " << v[V_KT] << " dp " << pressure-v[V_P] << " ptest " << pressure - v[V_P]*(1.+1.e-9) << std::endl;
    throw std::runtime_error("RegionRT::getRatPT: pressure too small at high density bound");
  }

  double llr = maxr;
  double llp = v[V_P];
  double rmax = maxr;

  //
  // secant/bisection search for R
  //
  while (fabs(rmax-rmin) > 1.e-14*fabs(rmax+rmin)/2.) {
    double r = lr+(pressure-lp)*(lr-llr)/(lp-llp);
    if (r >= rmax || r <= rmin || lp-llp == 0. || v[V_KT] < 0.) {
      r = (rmax+rmin)/2.;
    }
    getvars(r, temperature, v);
    
    if (v[V_KT] < 0.) {
      if (side > 0) rmin = r;
      else if (side < 0) rmax = r;
      else {
        // side == 0 and pressure >= pcrit ensured by check at start of routine
        if (pressure > v[V_P]) rmin = r;
        else rmax = r;
      }
    }
    else if (pressure > v[V_P]) {
      rmin = r;
    }
    else {
      rmax = r;
    }
    llp = lp;
    lp = v[V_P];
    llr = lr;
    lr = r;
  }
  
  //
  // failed to meet the desired pressure
  //
  /*
  if (fabs(pressure-v[V_P]) > 1.e-8*fabs(pressure+v[V_P]))
    throw std::runtime_error("RegionRT::getRatPT: failed to converge to input pressure");
  */

  return rmax;
}
