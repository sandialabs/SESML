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

#include "IAPWS95.H"
#include "FiniteDiffDeriv.H"
#include "PhaseBoundaryInfo.H"

void IAPWS95::getvars2( const double density,
                        const double temperature,
                        const Vars xtype,
                        std::vector<double> & v )
{
  // entries are:
  // 0 1  2    3    4    5
  // P X dPdR dPdT dXdR dXdT
  // for X in E,H,S,G
  v.resize(6);

  //
  // check xtype validity
  //
  if (xtype != V_H && xtype != V_E && xtype != V_S && xtype != V_G) {
    throw std::runtime_error("IAPWS95::getvars2: invalid variable type");
  }

  double s = density/rs;
  double t = ts/temperature;
  double ls = log(s);
  double lt = log(t);

  //
  // clip t values close to critical point to avoid division by zero
  //
  if (fabs(t-1.) < 1.e-14) {
    if (t < 1.) t = 1.-1.e-14;
    else t = 1.+1.e-14;
  }

  // ideal gas part density terms and first three temperature terms
  double vP = 1.;
  double vE = ci[1].n*t+ci[2].n;
  double cp1 = ci[2].n;
  double cp2 = 1.;
  double cp3 = 1.;
  double vH = 1.+ci[1].n*t+ci[2].n;
  double vS = ci[2].n*(1.-lt)-ls-ci[0].n;
  double vG = 1.+ls+ci[0].n+ci[1].n*t+ci[2].n*lt;
  // ideal gas part remaining temperature terms
  for (int i=3;i<8;i++) {
    double e = 1.-exp(-ci[i].t*t);
    vE += ci[i].n*ci[i].t*(1./e-1.)*t;
    cp1 += ci[i].n*ci[i].t*ci[i].t*(1.-e)/e/e*t*t;
    vH += ci[i].n*ci[i].t*(1./e-1.)*t;
    vS += ci[i].n*(ci[i].t*(1./e-1.)*t-log(e));
    vG += ci[i].n*log(e);
  }

  // first set of residual
  for (int i=0;i<7;i++) {
    double e = exp(cr[i].d*ls+cr[i].t*lt);
    vP += cr[i].n*e*cr[i].d;
    vE += cr[i].n*e*cr[i].t;
    cp1 -= cr[i].n*e*cr[i].t*(cr[i].t-1.);
    cp2 += cr[i].n*e*cr[i].d*(1.-cr[i].t);
    cp3 += cr[i].n*e*cr[i].d*(cr[i].d+1.);
    vH += cr[i].n*e*(cr[i].d+cr[i].t);
    vS += cr[i].n*e*(cr[i].t-1.);
    vG += cr[i].n*e*(cr[i].d+1.);
  }

  // second set of residual
  for (int i=7;i<51;i++) {
    double e = exp(cr[i].d*ls+cr[i].t*lt-exp(cr[i].c*ls));
    double dc = exp(cr[i].c*ls);
    vP += cr[i].n*e*(cr[i].d-cr[i].c*dc);
    vE += cr[i].n*e*cr[i].t;
    cp1 -= cr[i].n*e*cr[i].t*(cr[i].t-1.);
    cp2 += cr[i].n*e*(cr[i].d-cr[i].c*dc)*(1.-cr[i].t);
    cp3 += cr[i].n*e*((cr[i].d-cr[i].c*dc)*(1.+cr[i].d-cr[i].c*dc)-cr[i].c*cr[i].c*dc);
    vH += cr[i].n*e*(cr[i].d-cr[i].c*dc+cr[i].t);
    vS += cr[i].n*e*(cr[i].t-1.);
    vG += cr[i].n*e*(cr[i].d-cr[i].c*dc+1.);
  }

  // third set of residual
  for (int i=51;i<54;i++) {
    double e = exp(3.*ls+cr[i].c*lt-20.*(s-1.)*(s-1.)-cr[i].d*(t-cr[i].t)*(t-cr[i].t));
    vP += cr[i].n*e*(3.-s*40.*(s-1.));
    vE += cr[i].n*e*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t));
    cp1 -= cr[i].n*e*((cr[i].c-t*2.*cr[i].d*(t-cr[i].t))*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t))-cr[i].c-2.*t*t*cr[i].d);
    cp2 += cr[i].n*e*(3.-s*40.*(s-1.))*(1.-cr[i].c+t*2.*cr[i].d*(t-cr[i].t));
    cp3 += cr[i].n*e*(12.-40.*s*s+80.*s*(s-1.)*(20.*s*(s-1.)-4.));
    vH += cr[i].n*e*(3.-s*40.*(s-1.)+cr[i].c-t*2.*cr[i].d*(t-cr[i].t));
    vS += cr[i].n*e*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t)-1.);
    vG += cr[i].n*e*(3.-s*40.*(s-1.)+1.);
  }

  // fourth set of residual
  for (int i=54;i<56;i++) {
    double lss = log((s-1.)*(s-1.));
    double theta = 1.-t+0.32*exp(lss/2./0.3);
    double delta = theta*theta+0.2*exp(3.5*lss);
    double sodddds = s*(s-1.)/delta*(0.32*theta*2./0.3*exp((1./2./0.3-1.)*lss)+2.*0.2*3.5*exp((3.5-1.)*lss));
    double e = exp(cr[i].t*log(delta)+ls-cr[i].c*(s-1.)*(s-1.)-cr[i].d*(t-1.)*(t-1.));
    vP += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds);
    vE += -2.*cr[i].n*e*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d);
    cp1 -= cr[i].n*e*2.*t*t*(cr[i].t/delta+2.*theta*theta*cr[i].t*(cr[i].t-1.)/delta/delta+4.*theta*cr[i].t/delta*cr[i].d*(t-1.)+(2.*cr[i].d*(t-1.)*(t-1.)-1.)*cr[i].d);
    cp2 += cr[i].n*e*(1.-2.*cr[i].c*(s-1.)*s+(cr[i].t*(s-1.)*s*(1.4*exp(2.5*lss)+32./15.*exp(2./3.*lss)*theta))/delta-t*(-((32./15.*cr[i].t*exp(2./3.*lss)*(s-1.)*s)/delta)-2.*cr[i].d*(t-1.)+4.*cr[i].c*cr[i].d*(s-1.)*s*(t-1.)-(2.*cr[i].t*theta)/delta+(4.*cr[i].t*cr[i].c*(s-1.)*s*theta)/delta-1./delta/delta*64./15.*(cr[i].t-1.)*cr[i].t*(s-1.)*s*theta*(21./32.*exp(2.5*lss)+exp(2./3.*lss)*theta)-1./delta*2.8*cr[i].t*cr[i].d*(s-1.)*s*(t-1.)*(exp(2.5*lss)+32./21.*exp(2./3.*lss)*theta)));
    cp3 += cr[i].n*e*(2.*(1.+cr[i].c*s*(4.+(-5.+2.*cr[i].c*(s-1.)*(s-1.))*s))+(cr[i].t*s*(512./225.*exp(7./3.*lss)*s+exp(2.5*lss)*(-5.6+(14.-5.6*cr[i].c*(s-1.)*(s-1.))*s)+exp(2./3.*lss)*(-128./15.+(608./45.-128./15.*cr[i].c*(s-1.)*(s-1.))*s)*theta))/delta+(cr[i].t*(s-1.)*(s-1.)*s*s*(1.96*(cr[i].t-1.)*exp(5.*lss)+theta*(448./75.*(cr[i].t-1.)*exp(19./6.*lss)+1024./225.*(cr[i].t-1.)*exp(4./3.*lss)*theta)))/delta/delta);
    vH += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds-2.*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d));
    vS += cr[i].n*e*(-2.*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d)-1.);
    vG += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds+1.);
  }

  // P
  v[0] = vP*density*temperature*R;
  // dPdR
  v[2] = cp3*R*temperature;
  // dPdT
  v[3] = R*density*cp2;
  if (xtype == V_E) {
    // E
    v[1] = vE*temperature*R;
    // dEdR = (P-TdPdT)/R^2
    v[4] = R*(vP-cp2)*temperature/density;
    // dEdT
    v[5] = R*cp1;
  }
  else if (xtype == V_H) {
    // H
    v[1] = vH*temperature*R;
    // dHdR = (RdPdR-TdPdT)/R^2
    v[4] = R*(cp3-cp2)*temperature/density;
    // dHdT = dEdT+dPdT/R
    v[5] = R*(cp1+cp2);
  }
  else if (xtype == V_S) {
    // S
    v[1] = vS*R;
    // dSdR = -dPdT/R^2
    v[4] = -R*cp2/density;
    // dSdT = dE/dT/T
    v[5] = R*cp1/temperature;
  }
  else { // xtype == V_G
    // G
    v[1] = vG*temperature*R;
    // dGdR = dPdR/R
    v[4] = cp3*R*temperature/density;
    // dGdT = -S+dP/dT/R
    v[5] = (cp2-vS)*R;
  }
}

void IAPWS95::getvars( const double density,
                       const double temperature,
                       std::vector<double> & v)
{
  v.resize(V_MAX);

  double s = density/rs;
  double t = ts/temperature;
  double ls = log(s);
  double lt = log(t);

  //
  // clip t values close to critical point to avoid division by zero
  //
  if (fabs(t-1.) < 1.e-14) {
    if (t < 1.) t = 1.-1.e-14;
    else t = 1.+1.e-14;
  }

  // ideal gas part density terms and first three temperature terms
  v[V_P] = 1.;
  v[V_E] = ci[1].n*t+ci[2].n;
  v[V_F] = ls+ci[0].n+ci[1].n*t+ci[2].n*lt;
  double cp1 = ci[2].n;
  double cp2 = 1.;
  double cp3 = 1.;
  v[V_H] = 1.+ci[1].n*t+ci[2].n;
  v[V_S] = ci[2].n*(1.-lt)-ls-ci[0].n;
  v[V_G] = 1.+ls+ci[0].n+ci[1].n*t+ci[2].n*lt;
  // ideal gas part remaining temperature terms
  for (int i=3;i<8;i++) {
    double e = 1.-exp(-ci[i].t*t);
    v[V_E] += ci[i].n*ci[i].t*(1./e-1.)*t;
    v[V_F] += ci[i].n*log(e);
    cp1 += ci[i].n*ci[i].t*ci[i].t*(1.-e)/e/e*t*t;
    v[V_H] += ci[i].n*ci[i].t*(1./e-1.)*t;
    v[V_S] += ci[i].n*(ci[i].t*(1./e-1.)*t-log(e));
    v[V_G] += ci[i].n*log(e);
  }

  //std::cout << std::setprecision(9);
  //std::cout << "r " << density << " t " << temperature << " p " << v[V_P]/s << " e " << v[V_E]/t << " f " << v[V_F] << " cv " << -cp1/t/t << " cp2 " << cp2 << " cp3 " << cp3 << std::endl;

  // first set of residual
  for (int i=0;i<7;i++) {
    double e = exp(cr[i].d*ls+cr[i].t*lt);
    //std::cout << "i " << i << " e " << e << " lt " << lt << " ls " << ls << " c " << cr[i].c << " d " << cr[i].d << " t " << cr[i].t << std::endl;
    //std::cout << "i " << i << " e " << e << " n " << cr[i].n << std::endl;
    v[V_F] += cr[i].n*e;
    v[V_P] += cr[i].n*e*cr[i].d;
    v[V_E] += cr[i].n*e*cr[i].t;
    cp1 -= cr[i].n*e*cr[i].t*(cr[i].t-1.);
    cp2 += cr[i].n*e*cr[i].d*(1.-cr[i].t);
    cp3 += cr[i].n*e*cr[i].d*(cr[i].d+1.);
    v[V_H] += cr[i].n*e*(cr[i].d+cr[i].t);
    v[V_S] += cr[i].n*e*(cr[i].t-1.);
    v[V_G] += cr[i].n*e*(cr[i].d+1.);
  }
  //std::cout << "r " << density << " t " << temperature << " p " << v[V_P]/s << " e " << v[V_E]/t << " f " << v[V_F] << " cv " << -cp1/t/t << " cp2 " << cp2 << " cp3 " << cp3 << std::endl;

  // second set of residual
  for (int i=7;i<51;i++) {
    double e = exp(cr[i].d*ls+cr[i].t*lt-exp(cr[i].c*ls));
    double dc = exp(cr[i].c*ls);
    //std::cout << "i " << i << " e " << e << " lt " << lt << " ls " << ls << " c " << cr[i].c << " d " << cr[i].d << " t " << cr[i].t << std::endl;
    //std::cout << "i " << i << " e " << e << " n " << cr[i].n << std::endl;
    v[V_F] += cr[i].n*e;
    v[V_P] += cr[i].n*e*(cr[i].d-cr[i].c*dc);
    v[V_E] += cr[i].n*e*cr[i].t;
    cp1 -= cr[i].n*e*cr[i].t*(cr[i].t-1.);
    cp2 += cr[i].n*e*(cr[i].d-cr[i].c*dc)*(1.-cr[i].t);
    cp3 += cr[i].n*e*((cr[i].d-cr[i].c*dc)*(1.+cr[i].d-cr[i].c*dc)-cr[i].c*cr[i].c*dc);
    v[V_H] += cr[i].n*e*(cr[i].d-cr[i].c*dc+cr[i].t);
    v[V_S] += cr[i].n*e*(cr[i].t-1.);
    v[V_G] += cr[i].n*e*(cr[i].d-cr[i].c*dc+1.);
  }
  //std::cout << "r " << density << " t " << temperature << " p " << v[V_P]/s << " e " << v[V_E]/t << " f " << v[V_F] << " cv " << -cp1/t/t << " cp2 " << cp2 << " cp3 " << cp3 << std::endl;

  // third set of residual
  for (int i=51;i<54;i++) {
    double e = exp(3.*ls+cr[i].c*lt-20.*(s-1.)*(s-1.)-cr[i].d*(t-cr[i].t)*(t-cr[i].t));
    //std::cout << "i " << i << " e " << e << " n " << cr[i].n << std::endl;
    v[V_F] += cr[i].n*e;
    v[V_P] += cr[i].n*e*(3.-s*40.*(s-1.));
    v[V_E] += cr[i].n*e*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t));
    cp1 -= cr[i].n*e*((cr[i].c-t*2.*cr[i].d*(t-cr[i].t))*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t))-cr[i].c-2.*t*t*cr[i].d);
    cp2 += cr[i].n*e*(3.-s*40.*(s-1.))*(1.-cr[i].c+t*2.*cr[i].d*(t-cr[i].t));
    cp3 += cr[i].n*e*(12.-40.*s*s+80.*s*(s-1.)*(20.*s*(s-1.)-4.));
    v[V_H] += cr[i].n*e*(3.-s*40.*(s-1.)+cr[i].c-t*2.*cr[i].d*(t-cr[i].t));
    v[V_S] += cr[i].n*e*(cr[i].c-t*2.*cr[i].d*(t-cr[i].t)-1.);
    v[V_G] += cr[i].n*e*(3.-s*40.*(s-1.)+1.);
  }
  //std::cout << "r " << density << " t " << temperature << " p " << v[V_P]/s << " e " << v[V_E]/t << " f " << v[V_F] << " cv " << -cp1/t/t << " cp2 " << cp2 << " cp3 " << cp3 << std::endl;

  // fourth set of residual
  for (int i=54;i<56;i++) {
    double lss = log((s-1.)*(s-1.));
    double theta = 1.-t+0.32*exp(lss/2./0.3);
    double delta = theta*theta+0.2*exp(3.5*lss);
    double sodddds = s*(s-1.)/delta*(0.32*theta*2./0.3*exp((1./2./0.3-1.)*lss)+2.*0.2*3.5*exp((3.5-1.)*lss));
    double e = exp(cr[i].t*log(delta)+ls-cr[i].c*(s-1.)*(s-1.)-cr[i].d*(t-1.)*(t-1.));
    //std::cout << "i " << i << " e " << e << " n " << cr[i].n << std::endl;
    v[V_F] += cr[i].n*e;
    v[V_P] += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds);
    v[V_E] += -2.*cr[i].n*e*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d);
    cp1 -= cr[i].n*e*2.*t*t*(cr[i].t/delta+2.*theta*theta*cr[i].t*(cr[i].t-1.)/delta/delta+4.*theta*cr[i].t/delta*cr[i].d*(t-1.)+(2.*cr[i].d*(t-1.)*(t-1.)-1.)*cr[i].d);
    cp2 += cr[i].n*e*(1.-2.*cr[i].c*(s-1.)*s+(cr[i].t*(s-1.)*s*(1.4*exp(2.5*lss)+32./15.*exp(2./3.*lss)*theta))/delta-t*(-((32./15.*cr[i].t*exp(2./3.*lss)*(s-1.)*s)/delta)-2.*cr[i].d*(t-1.)+4.*cr[i].c*cr[i].d*(s-1.)*s*(t-1.)-(2.*cr[i].t*theta)/delta+(4.*cr[i].t*cr[i].c*(s-1.)*s*theta)/delta-1./delta/delta*64./15.*(cr[i].t-1.)*cr[i].t*(s-1.)*s*theta*(21./32.*exp(2.5*lss)+exp(2./3.*lss)*theta)-1./delta*2.8*cr[i].t*cr[i].d*(s-1.)*s*(t-1.)*(exp(2.5*lss)+32./21.*exp(2./3.*lss)*theta)));
    cp3 += cr[i].n*e*(2.*(1.+cr[i].c*s*(4.+(-5.+2.*cr[i].c*(s-1.)*(s-1.))*s))+(cr[i].t*s*(512./225.*exp(7./3.*lss)*s+exp(2.5*lss)*(-5.6+(14.-5.6*cr[i].c*(s-1.)*(s-1.))*s)+exp(2./3.*lss)*(-128./15.+(608./45.-128./15.*cr[i].c*(s-1.)*(s-1.))*s)*theta))/delta+(cr[i].t*(s-1.)*(s-1.)*s*s*(1.96*(cr[i].t-1.)*exp(5.*lss)+theta*(448./75.*(cr[i].t-1.)*exp(19./6.*lss)+1024./225.*(cr[i].t-1.)*exp(4./3.*lss)*theta)))/delta/delta);
    v[V_H] += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds-2.*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d));
    v[V_S] += cr[i].n*e*(-2.*t*(theta*cr[i].t/delta+(t-1.)*cr[i].d)-1.);
    v[V_G] += cr[i].n*e*(1.-2.*cr[i].c*s*(s-1.)+cr[i].t*sodddds+1.);
  }
  //std::cout << "r " << density << " t " << temperature << " p " << v[V_P]/s << " e " << v[V_E]/t << " f " << v[V_F] << " cv " << -cp1/t/t << " cp2 " << cp2 << " cp3 " << cp3 << std::endl;

  v[V_P] *= density*temperature*R;
  v[V_H] *= temperature*R;
  v[V_E] *= temperature*R;
  v[V_S] *= R;

  v[V_CV] = R*cp1;
  v[V_CP] = R*(cp1+cp2*cp2/cp3);
  v[V_KT] = cp3*R*temperature*density;

  // convert to m^2/s^2 from kJ/kg
  v[V_CS] = R*temperature*(cp3+cp2*cp2/cp1)*1000.;
  //if (v[V_CS] < 0.) throw std::runtime_error("IAPWS95::getvars: sound speed imaginary");
  v[V_CS] = sqrt(v[V_CS]);

  v[V_F] *= temperature*R;
  v[V_G] *= temperature*R;

  v[V_R] = density;
  v[V_V] = 1./density;
  v[V_T] = temperature;
}

void IAPWS95::getvarsRT( const double density,
                         const double temperature,
                         const Phases phase,
                         std::vector<double> & v )
{
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  if (phase == PH_LGMIX) throw std::runtime_error("IAPWS95::getvarsRT: liquid-gas mixture not supported");

  getvars(density,temperature,v);
}

void IAPWS95::getvarsPT( const double pressure,
                         const double temperature,
                         const Phases phase,
                         std::vector<double> & v )
{
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  if (phase == PH_LGMIX) throw std::runtime_error("IAPWS95::getvarsPT: cannot evaluate liquid-gas mixture in P-T space");

  //
  // out of bounds pressure check
  //
  if (pressure > 1.e6*(1.+1.e-10) || pressure < 0.611657) {
    throw std::runtime_error("IAPWS95::getvarsPT: pressure outside range of validity");
  }

  double minr(0.001),maxr(10000.);
  double rho;
  if (temperature < ts) {
    //
    // deal with liquid-vapor transition
    //
    // do it approximately as this is much faster
    double rho1,rho2;
    rho2 = getSpinodalDensityBound( temperature, true );
    rho1 = getSpinodalDensityBound( temperature, false );
    double tsat(temperature);
    if (pressure < 22064.) tsat = r4.getT(pressure);
    if (temperature > tsat) {
      maxr = rho1;
      rho = getRatPT( pressure, temperature, minr, maxr, -1 );
    }
    else {
      minr = rho2;
      rho = getRatPT( pressure, temperature, minr, maxr, 1 );
    }
  }
  else {
    //
    // above liquid-vapor transition so just get the state
    //
    rho = getRatPT( pressure, temperature, minr, maxr, 0);
  }

  getvars(rho,temperature,v);
  computeTransportProperties(v);
  if (isnan(v[V_M])) {
    throw std::runtime_error("IAPWS95::getvarsPT: transport property is nan");
  }
}

double IAPWS95::getSpinodalDensityBound( const double temperature,
                                         const bool glstate )
{
  //
  // clip to critical density when close enough to critical point
  //
  if (temperature > 647.04) return rs;

  const double & x(temperature);

  if (glstate) {
    //
    // spinodal line density on liquid side -- 0.99 times to give lower bound
    //
    double ra = 47.;
    double rc = 331.;
    double tc = 647.096;
    double r0 = -1.78773275622362;
    double r1 = 0.00382177916509115;
    double r2 = -3.68875574092353e-06;
    return 0.99*(rc+sqrt(tc-x)*ra+(tc-x)*(r0+(tc-x)*(r1+(tc-x)*r2)));
  }
  else {
    //
    // spinodal line density on gas side -- 1.03 times to give upper bound
    //
    double tc = 647.096;
    double g0 = 5.75;
    double g1 = 0.00727588987900081;
    double g2 = -0.000112703444590623;
    double g3 = 1.90920000027463e-07;
    double g4 = 0.000423160550154122;
    double gs = -0.160407509437179;
    double g5 = -7.04954987320237e-06;
    return 1.03*exp((g0+sqrt(tc-x)*gs+(tc-x)*(g1+(tc-x)*(g2+(tc-x)*g3)))/(1.+(tc-x)*(g4+(tc-x)*g5)));
  }
}

void IAPWS95::getvarsPX( const double pressure,
                         const double xvar,
                         const Phases phase,
                         const Vars xtype,
                         std::vector<double> & v )
{
  //std::cerr << std::setprecision(18);
  //std::cerr << "IAPWS95::getvarsPX: p " << pressure << " xvar " << xvar << " phase " << phase << " xtype " << xtype << std::endl;
  //std::cerr << std::setprecision(6);
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  //
  // check xtype validity
  //
  if (xtype != V_H && xtype != V_E && xtype != V_S) {
    throw std::runtime_error("IAPWS95::getvarsPX: invalid variable type");
  }

  //
  // out of bounds pressure check
  //
  if (pressure > 1.e6*(1.+1.e-10) || pressure < 0.611657) {
    throw std::runtime_error("IAPWS95::getvarsPX: pressure outside range of validity");
  }

  //
  // give 50 K of extrapolation room
  //
  double minT = 273.16;
  if (Textrap) minT -= 10.;

  if (phase & (PH_GAS | PH_LIQUID | PH_FLUID)) {
    //std::cerr << " fluid phase " << phase << std::endl;
    double r,t;
    //
    // if close to the critical pressure, start the temperature closer
    // to the critical temperature to aid in convergence
    //
    if (fabs(pressure-22064.) < 0.1)
      getRTatPX(pressure,xvar,xtype,rs,ts+1.,r,t);
    else
      getRTatPX(pressure,xvar,xtype,rs,ts+100.,r,t);
    getvars(r,t,v);
    //std::cout << "cp " << v[V_CP] << " kt " << v[V_KT] << std::endl;
    computeTransportProperties(v);
    if (isnan(v[V_M])) {
      throw std::runtime_error("IAPWS95::getvarsPX: transport property is nan");
    }
  }
  else if (phase == PH_LGMIX) {
    if (pressure > 22064.*(1.+1.e-10)) {
      throw std::runtime_error("IAPWS95::getvarsPX: pressure too high for liquid-gas mixture");
    }
    else {
      getMaxwellPX( pressure, xvar, xtype, v );
    }

    //
    // transport calculation happens inside maxwell call
    //
  }
  else {
    throw std::runtime_error("IAPWS95::getvarsPX: invalid phase");
  }
}

void inverse3x3(const std::vector<double> a,std::vector<double> & b)
{
  if (a.size() != 9) throw std::runtime_error("inverse3x3: input matrix has wrong size");
  b.resize(9);
  b[0] = a[4]*a[8]-a[5]*a[7];
  b[3] = -a[3]*a[8]+a[5]*a[6];
  b[6] = a[3]*a[7]-a[4]*a[6];
  b[1] = -a[1]*a[8]+a[2]*a[7];
  b[4] = a[0]*a[8]-a[2]*a[6];
  b[7] = -a[0]*a[7]+a[1]*a[6];
  b[2] = a[1]*a[5]-a[2]*a[4];
  b[5] = -a[0]*a[5]+a[2]*a[3];
  b[8] = a[0]*a[4]-a[1]*a[3];
  double det = a[0]*b[0]+a[1]*b[3]+a[2]*b[6];
  if (fabs(det) < 1.e-14) throw std::runtime_error("inverse3x3: input matrix is singular");
  double deti = 1./det;
  for (int i=0;i<9;i++) b[i] *= deti;
}

void getcritdensities(const double temperature,double & rho1,double & rho2)
{
  rho1 = 322.-168.61*sqrt(647.096-temperature)*(1.+940./168.61*(647.096-temperature));
  rho2 = 322.+168.61*sqrt(647.096-temperature)*(1.-940./168.61*(647.096-temperature));
}

bool IAPWS95::getSatPropsNearCritical( const double pressure,
                                       double & density1,
                                       double & density2,
                                       double & temperature )
{
  std::vector<double> v;
  double rho1,rho2;

  //
  // temperature range
  //
  double tmin = 647.09585;
  double tmax = 647.096;

  //
  // check the pressure at the temperature minimum
  //
  getcritdensities(tmin,rho1,rho2);
  getvars2(rho1,tmin,V_E,v);

  //
  // return false if pressure out of range
  //
  if (pressure <= v[0]) return false;

  //
  // first guess
  //
  double lt = tmin;
  double lp = v[0];

  //
  // second guess
  //
  double llt = (tmin+tmax)*0.5;
  getcritdensities(llt,rho1,rho2);
  getvars2(rho1,llt,V_E,v);
  double llp = v[0];

  //
  // secant/bisection search for T
  //
  double t = llt;
  while (fabs(tmax-tmin) > 1.e-14*fabs(tmax+tmin)/2.) {
    t = lt+(pressure-lp)*(lt-llt)/(lp-llp);
    if (t >= tmax || t <= tmin || lp-llp == 0.0) {
      t = (tmax+tmin)/2.;
    }
    getcritdensities(t,rho1,rho2);
    getvars2(rho1,t,V_E,v);
    double p = v[0];
    
    if (p < pressure) tmin = t;
    else tmax = t;
    llt = lt;
    lt = t;
    llp = lp;
    lp = p;
  }

  temperature = t;
  density1 = rho1;
  density2 = rho2;

  return true;
}

void IAPWS95::getSatProps( const double pressure,
                           double & density1,
                           double & density2,
                           double & temperature )
{
  //
  // check for valid pr<essure
  //
  if (pressure < 0.5 || pressure > 22064.*(1.+1.e-10)) {
    throw std::runtime_error("IAPWFS95::getSatProps: pressure out of bounds");
  }

  //
  // check if in the asymptotic critical regime
  //
  if (getSatPropsNearCritical(pressure,density1,density2,temperature)) return;

  int debug = 0;

  //
  // starting temperature estimate from region 4
  //
  double t = r4.getT(pressure);

  //
  // temperature bounds
  //
  double tmin = std::max(273.16,t-1.);
  double tmax = std::min(647.096,t+1.);

  //
  // starting density at midpoint of bounds
  //
  double r1min = 1.6e-6*exp(.0287*t);
  double r1max = std::min(rs,getSpinodalDensityBound(t,false));
  double r2min = std::max(rs,getSpinodalDensityBound(t,true));
  double r2max = 1350.-t;
  double r1 = r1min;
  double r2 = r2max;
  if (debug > 0) std::cout << "start: t " << t << " r1 " << r1 << " r2 " << r2 << std::endl;

  double alpha;
  double r1a,r2a,ta;

  std::vector<double> v1,v2;
  double dp1,dp2,dg;
  double dr1(1.0),dr2(1.0),dt(1.0);

  int i=0;
  int j=0;
  getvars2(r1,t,V_G,v1);
  getvars2(r2,t,V_G,v2);
  std::vector<double> jac(9),ijac(9);
  while(i<100) {
    // function values
    dp1 = v1[0]-pressure;
    dp2 = v2[0]-pressure;
    dg = v1[1]-v2[1];
    if (debug > 0) std::cout << " dp1 " << fabs(dp1/pressure) << " dp2 " << fabs(dp2/pressure) << " dg " << fabs(dg/(v1[1]+v2[1])) << std::endl;
    if (fabs(dr1) < 1.e-9*r1 && fabs(dr2) < 1.e-9*r2 && fabs(dt) < 1.e-9*t
        && fabs(dg) < 1.e-9*fabs(v1[1]+v2[1])) break;
    // jacobian
    if (fabs(dp1/pressure) > 1.e-7 || fabs(dp2/pressure) > 1.e-7) {
      dr1 = dp1/v1[2];
      dr2 = dp2/v2[2];
      dt = 0.0;
    }
    else {
      jac[0] = v1[2];
      jac[1] = 0.0;
      jac[2] = v1[3];
      jac[3] = 0.0;
      jac[4] = v2[2];
      jac[5] = v2[3];
      jac[6] = v1[4];
      jac[7] = -v2[4];
      jac[8] = v1[5]-v2[5];
      inverse3x3(jac,ijac);
      dr1 = ijac[0]*dp1+ijac[1]*dp2+ijac[2]*dg;
      dr2 = ijac[3]*dp1+ijac[4]*dp2+ijac[5]*dg;
      dt  = ijac[6]*dp1+ijac[7]*dp2+ijac[8]*dg;
    }
    alpha = 1.0;
    j=0;
    while(j<10) {
      r1a = r1-alpha*dr1;
      r2a = r2-alpha*dr2;
      ta = t-alpha*dt;
      if (ta < tmin || ta > tmax) ta = t;
      r1min = 1.6e-6*exp(.0287*ta);
      r1max = std::min(rs,getSpinodalDensityBound(ta,false));
      r2min = std::max(rs,getSpinodalDensityBound(ta,true));
      r2max = 1350.-ta;
      v1[2] = -9.e99;
      if (r1a > r1min && r1a < r1max && r2a > r2min && r2a < r2max) {
        //&& ta > tmin && ta < tmax) {
        getvars2(r1a,ta,V_G,v1);
        getvars2(r2a,ta,V_G,v2);
      }
      if (debug > 0) std::cout << "i " << i << " j " << j << " r1min " << r1min << " r1max " << r1max << " r2min " << r2min << " r2max " << r2max << " alpha " << alpha << " r1a " << r1a << " r2a " << r2a << " ta " << ta << " kt1 " << v1[2] << " kt2 " << v2[2] << std::endl;
      if (v1[2] < 0.0 || v2[2] < 0.0 || v1[0] < 0.0 || v2[0] < 0.0) alpha *= 0.1;
      else break;
      j++;
    }
    r1 = r1a;
    r2 = r2a;
    t = ta;
    if (debug > 0) std::cout << "i " << i << " j " << j << " ptgt " << pressure << " p1 " << v1[0] << " p2 " << v2[0] << " dg " << v1[1]-v2[1] << " r1 " << r1 << " dr1 " << dr1 << " r2 " << r2 << " dr2 " << dr2 << " t " << t << " dt " << dt << " kt1 " << v1[2] << " kt2 " << v2[2] << std::endl;
    if (j >= 10) break;
    i++;
  }

  density1 = r1;
  density2 = r2;
  temperature = t;

  if (j >= 10) {
    if (debug > 0) std::cout << "dr1 " << fabs(dr1) << " dr1t " << 1.e-14*r1 << " dr2 " << fabs(dr2) << " dr2t " << 1.e-14*r2 << " dt " << fabs(dt) << " dtt " << 1.e-14*t << " dp1 " << dp1/pressure << " dp2 " << dp2/pressure << " dg " << dg/(v1[1]+v2[1]) << std::endl;
    if (fabs(dp1/pressure) > 1.e-5 || fabs(dp2/pressure) > 1.e-5 || fabs(dg/(v1[1]+v2[1])) > 1.e-14)
      throw std::runtime_error("IAPWS95::getSatProps: unable to take valid step");
  }
  if (i >= 100) {
    // reduce dg tolerance when g is near zero
    double epsmod = std::numeric_limits<double>::epsilon();
    epsmod = 2.e6*std::max(epsmod,epsmod/(fabs(v1[1]+v2[1])+epsmod));
    if (debug > 0) std::cout << "dr1 " << fabs(dr1) << " dr1t " << 1.e-14*r1 << " dr2 " << fabs(dr2) << " dr2t " << 1.e-14*r2 << " dt " << fabs(dt) << " dtt " << 1.e-14*t << " dp1 " << dp1/pressure << " dp2 " << dp2/pressure << " g2 " << v2[1] << " dg " << dg/(v1[1]+v2[1]) << " epsmod " << epsmod << std::endl;
    if (fabs(dp1/pressure) > 1.e-5 || fabs(dp2/pressure) > 1.e-5 || fabs(dg/(v1[1]+v2[1])) > epsmod)
      throw std::runtime_error("IAPWS95::getSatProps: too many iteration without convergence");
  }
}

void IAPWS95::getMaxwellPX( const double pressure,
                            //const double temperature,
                            const double xvar,
                            const Vars xtype,
                            std::vector<double> & v )
{
  //
  // initialize outputs
  //
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  if (xtype != V_S && xtype != V_H && xtype != V_E && xtype != V_V) {
    throw std::runtime_error("IAPWS95::getMaxwellPX: invalid variable type");
  }

  //
  // Compute variables at Maxwell construction boundaries
  //
  std::vector<double> v1,v2;
  double rho1,rho2,temperature;

  getSatProps(pressure,rho1,rho2,temperature);
  getvars(rho1,temperature,v1);
  getvars(rho2,temperature,v2);

  //
  // get transport properties at construction boundaries
  //
  computeTransportProperties(v1);
  computeTransportProperties(v2);
  if (isnan(v1[V_M]) || isnan(v2[V_M])) {
    throw std::runtime_error("IAPWS95::getMaxwellPX: transport property is nan");
  }

  //
  // Compute lever rule for given variable
  //
  double f = (xvar-v2[xtype])/(v1[xtype]-v2[xtype]);
  if (fabs(v1[xtype]-v2[xtype]) < 1.e-10*(v1[xtype]+v2[xtype])) f = 1.;

  //
  // Compute variables at f
  //
  v[V_V] = f/rho1+(1.-f)/rho2;
  v[V_R] = 1./v[V_V];
  v[V_T] = temperature;
  v[V_P] = pressure;
  v[V_G] = v1[V_G]*f+v2[V_G]*(1.-f);
  v[V_F] = v1[V_F]*f+v2[V_F]*(1.-f);
  v[V_H] = v1[V_H]*f+v2[V_H]*(1.-f);
  v[V_E] = v1[V_E]*f+v2[V_E]*(1.-f);
  v[V_S] = v1[V_S]*f+v2[V_S]*(1.-f);
  v[V_CP] = 1.e200;//std::numeric_limits<double>::max();
  v[V_KT] = 0.;
  // probably not correct, but for now use lever rule on transport properties
  v[V_M] = v1[V_M]*f+v2[V_M]*(1.-f);
  v[V_L] = v1[V_L]*f+v2[V_L]*(1.-f);

  std::vector<double> vd1,vd2;
  getMaxwellDerivs( pressure, vd1, vd2 );
  //getMaxwellDerivs( temperature, vd1, vd2 );

  double dfdt = (-vd2[V_V]*(v1[V_V]-v2[V_V])-(v[V_V]-v2[V_V])*(vd1[V_V]-vd2[V_V]))/(v1[V_V]-v2[V_V])/(v1[V_V]-v2[V_V]);
  
  v[V_CV] = vd1[V_E]*f+vd2[V_E]*(1-f)+(v1[V_E]-v2[V_E])*dfdt;
  v[V_CS] = temperature*v[V_V]*v[V_V]*vd1[V_P]*vd1[V_P]/v[V_CV];
  if (v[V_CS] < 0.0) v[V_CS] = 0.0;
  else v[V_CS] = sqrt(v[V_CS]);

  //
  // form matched to cv points calculated at 647.09226064 K
  // cvc is computed from the scaling relations for E and R
  // cva, cvb, and cvd from fit to 3 of the 4 cv points calculated
  // at 647.08852128 and 647.09226064 K.
  //
  double cvc = 128.5;
  double cva = 23.2552;
  double cvb = -118.1915;
  double cvd = 0.0305136;
  if (temperature > 647.09226064) {
    double cve1 = cvc+cva*sqrt(647.096-temperature)+cvb*exp(cvd*log(647.096-temperature));
    double cve2 = cvc-cva*sqrt(647.096-temperature)+cvb*exp(cvd*log(647.096-temperature));
    v[V_CV] = cve1*f+cve2*(1.-f);    
    v[V_CS] = sqrt(temperature*v[V_V]*v[V_V]*vd1[V_P]*vd1[V_P]/v[V_CV]);
  }

  //
  // overwrite xtype var
  //
  v[xtype] = xvar;

}

void IAPWS95::getMaxwellDerivs( const double pressure,
                                std::vector<double> & vd1,
                                std::vector<double> & vd2 )
{
  //
  // initialize outputs
  //
  vd1.resize(V_MAX);
  vd2.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) {
    vd1[i] = 0.;
    vd2[i] = 0.;
  }

  //
  // make temperature grid
  //
  std::vector<double> p(5,0.);

  //
  // find location in stencil
  //
  int lp = 0;
  double dp = 1.e-6;
  if (pressure*(1.+2.*dp) < 22064.) lp = 2;
  else if (pressure*(1.+dp) < 22064.) lp = 1;
  //lp = 0;

  for (int i=0;i<5;i++) p[i] = pressure*(1.-dp*(i-lp));

  //
  // compute maxwell densities
  //
  std::vector<double> r1(5),r2(5),t(5);

  for (int i=0;i<5;i++) getSatProps(p[i],r1[i],r2[i],t[i]);

  //
  // compute states on boundaries
  //
  std::vector<std::vector<double> > v1(5),v2(5);

  for (int i=0;i<5;i++) {
    getvars(r1[i],t[i],v1[i]);
    getvars(r2[i],t[i],v2[i]);
  }

  //
  // compute derivatives
  //
  double dtdp = FivePointStencilDeriv(t,dp*pressure,lp,1);
  std::vector<double> f(5);
  for (int i=0;i<V_MAX;i++) {
    for (int j=0;j<5;j++) f[j] = v1[j][i];
    vd1[i] = FivePointStencilDeriv(f,dp*pressure,lp,1)/dtdp;
    for (int j=0;j<5;j++) f[j] = v2[j][i];
    vd2[i] = FivePointStencilDeriv(f,dp*pressure,lp,1)/dtdp;
  }
}

const std::string IAPWS95::name = "IAPWS95";

const EOSParamInfo IAPWS95::ParamList[] = {
  { "TEMP_EXTRAP", "IParam", TEMP_EXTRAP }
};
const int IAPWS95::ParamListSize = sizeof(ParamList)/sizeof(EOSParamInfo);

EOSModel * IAPWS95::createModel( const Parameters & params,
                                 const std::string name,
                                 const EOSData & /*init_data*/ )
{
  return new IAPWS95(params,name);
}

IAPWS95::IAPWS95()
  : RegionRT( 322.0, 647.096, 0.46151805, 0.001, 10000., 270., 1280. ), EOSModel("IAPWS95"), ci(0), cr(0), Textrap(false)
{
  setCoefficients();
}

IAPWS95::IAPWS95( const Parameters & params,
                  const std::string name )
  : RegionRT( 322.0, 647.096, 0.46151805, 0.001, 10000., 270., 1280. ), EOSModel("IAPWS95"), ci(0), cr(0), Textrap(false)
{
  //
  // find parameters for the model
  //
  const EOSParam * myparams = getParams(params,name);

  //
  // set extrapolation flag
  //
  if (myparams->params_int[TEMP_EXTRAP] != 0) Textrap = true;

  setCoefficients();
}
/*
void IAPWS95::clip_derivatives( std::vector<double> & v)
{
  //
  // clip the heat capacity and isothermal bulk modulus as recommended
  // in the IAPWS 2011 thermal conductivity standard
  //
  if (v[V_KT] < 1.e-13*v[V_R]*22064.0/322.0) {
    //std::cout << "clipped kt from " << v[V_KT];
    v[V_KT] = 1.e-13*v[V_R]*22064.0/322.0;
    //std::cout << " to " << v[V_KT] << std::endl;
  }
  if (v[V_CP] < 0. || v[V_CP] > 0.461526*1.e13) {
    //std::cout << "clipped cp from " << v[V_CP];
    v[V_CP] = 0.461526*1.e13;
    //std::cout << " to " << v[V_CP] << std::endl;
  }


}
*/

void IAPWS95::convert_to_SI( std::vector<double> & v )
{
  //
  // var  |  input  | output
  //      |  units  | units
  //------------------------
  // V_R  | kg/m^3  | same
  // V_V  | m^3/kg  | same
  // V_T  |   K     | same
  // V_P  |  kPa    |  Pa
  // V_G  | kJ/kg   | J/kg
  // V_F  | kJ/kg   | J/kg
  // V_H  | kJ/kg   | J/kg
  // V_E  | kJ/kg   | J/kg
  // V_S  | kJ/kg/K | J/kg/K
  // V_CP | kJ/kg/K | J/kg/K
  // V_CV | kJ/kg/K | J/kg/K
  // V_CS |  m/s    | same
  // V_KT |  kPa    |  Pa
  // V_M  |  uPa*s  | Pa*s
  // V_L  | mW/m/K  | W/m/K
  //
  v[V_P] *= 1000;
  v[V_G] *= 1000;
  v[V_F] *= 1000;
  v[V_H] *= 1000;
  v[V_E] *= 1000;
  v[V_S] *= 1000;
  v[V_CP] *= 1000;
  v[V_CV] *= 1000;
  v[V_KT] *= 1000;
  v[V_M] *= 1.e-6;
  v[V_L] *= 1.e-3;
}

void IAPWS95::computeTransportProperties( std::vector<double> & v )
{
  //
  // compute the reference state
  //
  std::vector<double> vr;
  getvars(v[V_R],647.096*1.5,vr);
  
  //
  // Calculate transport properties
  //
  v[V_M] = visc.getViscosity(v[V_R],v[V_T],v[V_KT],vr[V_KT]);
  v[V_L] = thcon.getThermalConductivity(v[V_R],v[V_T],v[V_KT],vr[V_KT],v[V_CV],v[V_CP],v[V_M]);
}

void IAPWS95::getThermalEOS_PT_D( EOSData & data, const int phase )
{
  std::vector<double> v;
  getvarsPT(data.inputs[0]*0.001,data.inputs[1],static_cast<Phases>(phase),v);
  convert_to_SI(v);
  data.outputs[0] = v[V_R];
  data.outputs[1] = v[V_H];
  data.outputs[2] = v[V_E];
  data.outputs[3] = v[V_S];
  data.outputs[4] = v[V_G];
  data.outputs[5] = v[V_F];
  data.outputs[6] = v[V_CP];
  data.outputs[7] = v[V_CV];
  data.outputs[8] = v[V_CS];
  data.outputs[9] = v[V_KT];
  data.outputs[10] = phase;
  data.outputs[11] = v[V_M];
  data.outputs[12] = v[V_L];
}

void IAPWS95::getThermalEOS_PH_D( EOSData & data, const int phase )
{
  std::vector<double> v;
  getvarsPX(data.inputs[0]*0.001,data.inputs[1]*0.001,static_cast<Phases>(phase),V_H,v);
  convert_to_SI(v);
  data.outputs[0] = v[V_R];
  data.outputs[1] = v[V_T];
  data.outputs[2] = v[V_E];
  data.outputs[3] = v[V_S];
  data.outputs[4] = v[V_G];
  data.outputs[5] = v[V_F];
  data.outputs[6] = v[V_CP];
  data.outputs[7] = v[V_CV];
  data.outputs[8] = v[V_CS];
  data.outputs[9] = v[V_KT];
  data.outputs[10] = phase;
  data.outputs[11] = v[V_M];
  data.outputs[12] = v[V_L];
}

void IAPWS95::getThermalEOS_PS_D( EOSData & data, const int phase )
{
  std::vector<double> v;
  getvarsPX(data.inputs[0]*0.001,data.inputs[1]*0.001,static_cast<Phases>(phase),V_S,v);
  convert_to_SI(v);
  data.outputs[0] = v[V_R];
  data.outputs[1] = v[V_T];
  data.outputs[2] = v[V_E];
  data.outputs[3] = v[V_H];
  data.outputs[4] = v[V_G];
  data.outputs[5] = v[V_F];
  data.outputs[6] = v[V_CP];
  data.outputs[7] = v[V_CV];
  data.outputs[8] = v[V_CS];
  data.outputs[9] = v[V_KT];
  data.outputs[10] = phase;
  data.outputs[11] = v[V_M];
  data.outputs[12] = v[V_L];
}

void IAPWS95::getThermalEOS_PE_D( EOSData & data, const int phase )
{
  std::vector<double> v;
  getvarsPX(data.inputs[0]*0.001,data.inputs[1]*0.001,static_cast<Phases>(phase),V_E,v);
  convert_to_SI(v);
  data.outputs[0] = v[V_R];
  data.outputs[1] = v[V_T];
  data.outputs[2] = v[V_H];
  data.outputs[3] = v[V_S];
  data.outputs[4] = v[V_G];
  data.outputs[5] = v[V_F];
  data.outputs[6] = v[V_CP];
  data.outputs[7] = v[V_CV];
  data.outputs[8] = v[V_CS];
  data.outputs[9] = v[V_KT];
  data.outputs[10] = phase;
  data.outputs[11] = v[V_M];
  data.outputs[12] = v[V_L];
}

void IAPWS95::getThermalEOS_RT_D( EOSData & data, const int phase )
{
  std::vector<double> v;
  getvarsRT(data.inputs[0],data.inputs[1],static_cast<Phases>(phase),v);
  convert_to_SI(v);
  data.outputs[0] = v[V_P];
  data.outputs[1] = v[V_H];
  data.outputs[2] = v[V_E];
  data.outputs[3] = v[V_S];
  data.outputs[4] = v[V_G];
  data.outputs[5] = v[V_F];
  data.outputs[6] = v[V_CP];
  data.outputs[7] = v[V_CV];
  data.outputs[8] = v[V_CS];
  data.outputs[9] = v[V_KT];
  data.outputs[10] = phase;
  data.outputs[11] = v[V_M];
  data.outputs[12] = v[V_L];
}

//
// here input pressures are in Pa
//
void IAPWS95::calculatePhaseBoundaryInfo(const double lowP,const double highP)
{
  //
  // liquid-vapor dome is only phase boundary
  //

  //
  // if above the critical point, don't put on any boundaries
  //
  if (lowP >= 22064000.) return;
  
  //
  // setup bounds -- ignore density inputs and put in pressure instead of temperature
  //
  std::vector<double> v1(3),v2(3);
  v1[0] = 0.;
  v1[1] = 0.;
  v1[2] = lowP;
  v2[0] = 0.;
  v2[1] = 0.;
  v2[2] = highP;

  //
  // clip high pressure to critical point
  //
  if (highP > 22064000.) v2[2] = 22064000.;

  //
  // setup phases
  //
  std::vector<int> phases(5);
  phases[0] = PH_GAS;
  phases[1] = PH_LGMIX;
  phases[2] = PH_LIQUID;
  phases[3] = PH_LGMIX;
  phases[4] = PH_LGMIX;

  phase_boundaries.push_back(new PhaseBoundaryInfo(this,v1,v2,phases));
  
}


void IAPWS95::getPhaseCoexistence( const int phase, const double pressure, double & temperature, double * enthalpy )
{
  if (phase != PH_LGMIX) throw std::runtime_error("IAPWS95::getPhaseCoexistence: invalid phase");

  //
  // convert to internal units
  //
  double ipressure = pressure/1000.;
  
  if (ipressure < 0.611657 || ipressure > 22064.) {
    throw std::runtime_error("IAPWS95::getPhaseCoexistence: pressure out of valid range");
  }
  else {
    double rho1,rho2;
    getSatProps(ipressure,rho1,rho2,temperature);
    std::vector<double> v1,v2;
    getvars(rho1,temperature,v1);
    getvars(rho2,temperature,v2);
    convert_to_SI(v1);
    convert_to_SI(v2);
    enthalpy[0] = v1[V_H];
    enthalpy[1] = v2[V_H];
  }
}

