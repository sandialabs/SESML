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

#include "IAPWS-IF97.H"
#include "FiniteDiffDeriv.H"
#include "PhaseBoundaryInfo.H"

void Region1::getvars( const double pressure,
                       const double temperature,
                       std::vector<double> & v )
{
  v.resize(V_MAX);

  double p = pressure/ps;
  double t = ts/temperature;
  double lp = log(7.1-p);
  double lt = log(t-1.222);
  double pr = p/(7.1-p);
  double tr = t/(t-1.222);
  //std::cout << "pressure " << pressure << " temperature " << temperature << " p " << p << " t " << t << " lp " << lp << " lt " << lt << std::endl;
  v[V_V] = 0.;
  v[V_H] = 0.;
  v[V_E] = 0.;
  v[V_S] = 0.;
  v[V_CP] = 0.;
  double w1 = 0.;
  double w2 = 0.;
  v[V_G] = 0.;
  v[V_F] = 0.;
  for (int i=0;i<34;i++) {
    //std::cout << "n " << i << " i " << c[i].I << " j " << c[i].J << " n " << c[i].n << std::endl;
    double e = exp(c[i].I*lp+c[i].J*lt);
    v[V_V] -= c[i].n*c[i].I*e;
    v[V_H] += c[i].n*c[i].J*e;
    v[V_E] += c[i].n*(c[i].J*tr+c[i].I*pr)*e;
    v[V_S] += c[i].n*(c[i].J*tr-1)*e;
    v[V_CP] -= c[i].n*c[i].J*(c[i].J-1)*e;
    w1 += c[i].n*c[i].I*(c[i].J*tr-1)*e;
    w2 += c[i].n*c[i].I*(c[i].I-1)*e;
    v[V_G] += c[i].n*e;
    v[V_F] += c[i].n*(1+pr*c[i].I)*e;
  }
  
  v[V_V] *= pr;
  v[V_H] *= tr;
  v[V_CP] *= tr*tr;
  w1 *= pr;
  w2 *= pr*pr;
  
  // convert to m^2/s^2 from kJ/kg
  v[V_CS] = -R*temperature*v[V_V]*v[V_V]/(w1*w1/v[V_CP]+w2)*1000.;
  //if (v[V_CS] < 0.) throw std::runtime_error("Region1::getvars: sound speed imaginary");
  v[V_CS] = sqrt(v[V_CS]);
  
  v[V_CV] = R*(v[V_CP]+w1*w1/w2);
  v[V_KT] = -pressure*v[V_V]/w2;

  v[V_R] = pressure/temperature/R/v[V_V];
  v[V_V] *= temperature*R/pressure;
  v[V_H] *= temperature*R;
  v[V_E] *= temperature*R;
  v[V_S] *= R;
  v[V_CP] *= R;
  v[V_G] *= temperature*R;
  v[V_F] *= temperature*R;

  v[V_P] = pressure;
  v[V_T] = temperature;
  v[V_G2P] = w2*R*temperature/pressure/pressure;
}

void Region2::getvars( const double pressure,
                       const double temperature,
                       std::vector<double> & v )
{
  v.resize(V_MAX);

  if (temperature >= 2*ts) throw std::runtime_error("Region2::getvars: temperature too large");
  
  double p = pressure/ps;
  double t = ts/temperature;
  double lp = log(p);
  double lt = log(t);
  //std::cout << "pressure " << pressure << " temperature " << temperature << " p " << p << " t " << t << " lp " << lp << " lt " << lt << std::endl;
  v[V_V] = 1.;
  v[V_H] = 0.;
  v[V_E] = -1.;
  v[V_S] = -lp;
  v[V_CP] = 0.;
  double w1 = 1.;
  double w2 = -1.;
  v[V_G] = lp;
  v[V_F] = lp-1.;
  // ideal gas part
  for (int i=0;i<9;i++) {
    //std::cout << "n " << i << " i " << c[i].I << " j " << c[i].J << " n " << c[i].n << std::endl;
    double e = exp(c[i].J*lt);
    v[V_H] += c[i].n*c[i].J*e;
    v[V_E] += c[i].n*c[i].J*e;
    v[V_S] += c[i].n*(c[i].J-1)*e;
    v[V_CP] -= c[i].n*c[i].J*(c[i].J-1)*e;
    v[V_G] += c[i].n*e;
    v[V_F] += c[i].n*e;
  }
  
  lt = log(t-0.5);
  double tr = t/(t-0.5);
  //std::cout << "pressure " << pressure << " temperature " << temperature << " p " << p << " t " << t << " lp " << lp << " lt " << lt << std::endl;
  // residual part
  for (int i=9;i<52;i++) {
    //std::cout << "n " << i << " i " << c[i].I << " j " << c[i].J << " n " << c[i].n << std::endl;
    double e = exp(c[i].I*lp+c[i].J*lt);
    v[V_V] += c[i].n*c[i].I*e;
    v[V_H] += c[i].n*c[i].J*tr*e;
    v[V_E] += c[i].n*(c[i].J*tr-c[i].I)*e;
    v[V_S] += c[i].n*(c[i].J*tr-1)*e;
    v[V_CP] -= c[i].n*c[i].J*(c[i].J-1)*tr*tr*e;
    w1 += c[i].n*c[i].I*(1-c[i].J*tr)*e;
    w2 += c[i].n*c[i].I*(c[i].I-1)*e;
    v[V_G] += c[i].n*e;
    v[V_F] += c[i].n*(1-c[i].I)*e;
  }
  //std::cout << " w2 " << w2 << " cp " << v[V_CP] << " w1 " << w1 << " cs " << v[V_CS] << " v " << v[V_V] << " r " << R << std::endl;
  // convert to m^2/s^2 from kJ/kg
  v[V_CS] = -R*temperature*v[V_V]*v[V_V]/(w1*w1/v[V_CP]+w2)*1000.;
  //if (v[V_CS] < 0.) throw std::runtime_error("Region2::getvars: sound speed imaginary");
  v[V_CS] = sqrt(v[V_CS]);

  v[V_CV] = R*(v[V_CP]+w1*w1/w2);
  v[V_KT] = -pressure*v[V_V]/w2;
  
  v[V_R] = pressure/temperature/R/v[V_V];
  v[V_V] *= temperature*R/pressure;
  v[V_H] *= temperature*R;
  v[V_E] *= temperature*R;
  v[V_S] *= R;
  v[V_CP] *= R;
  v[V_G] *= temperature*R;
  v[V_F] *= temperature*R;

  v[V_P] = pressure;
  v[V_T] = temperature;
  v[V_G2P] = w2*R*temperature/pressure/pressure;

}

double Region3::getSpinodalDensityBound( const double temperature,
                                         const bool glstate )
{
  //
  // clip to critical density when close enough to critical point
  //
  if (temperature > 647.) return rs;

  if (glstate) {
    //
    // spinodal line density on liquid side -- this is a lower bound good down to 610 K
    //
    return 525.+(322.-525.)/(647.-610.)*(temperature-610.);
  }
  else {
    //
    // spinodal line density on gas side -- this is an upper bound good down to 610 K
    //
    return 150.+(322.-150.)/(647.-610.)*(temperature-610.);
  }
}

void Region3::getvars2( const double density,
                        const double temperature,
                        const Vars xtype,
                        std::vector<double> & v )
{
  // entries are:
  // 0 1  2    3    4    5
  // P X dPdR dPdT dXdR dXdT
  // for X in E,H,S
  v.resize(6);

  //
  // check xtype validity
  //
  if (xtype != V_H && xtype != V_E && xtype != V_S) {
    throw std::runtime_error("Region3::getvars2: invalid variable type");
  }

  double s = density/rs;
  double t = ts/temperature;
  double ls = log(s);
  double lt = log(t);
  
  double vP = c[0].n;
  double vH = c[0].n;
  double vE = 0.;
  double vS = -c[0].n*log(s);
  double cp1 = 0.;
  double cp2 = c[0].n;
  double cp3 = c[0].n;
  for (int i=1;i<40;i++) {
    double e = exp(c[i].I*ls+c[i].J*lt);
    vP += c[i].n*c[i].I*e;
    vH += c[i].n*(c[i].I+c[i].J)*e;
    vE += c[i].n*c[i].J*e;
    vS += c[i].n*(c[i].J-1)*e;
    cp1 -= c[i].n*c[i].J*(c[i].J-1)*e;
    cp2 += c[i].n*c[i].I*(1-c[i].J)*e;
    cp3 += c[i].n*c[i].I*(1+c[i].I)*e;
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
  else { // xtype == V_S
    // S
    v[1] = vS*R;
    // dSdR = -dPdT/R^2
    v[4] = -R*cp2/density;
    // dSdT = dE/dT/T
    v[5] = R*cp1/temperature;
  }

}

void Region3::getvars( const double density,
                       const double temperature,
                       std::vector<double> & v )
{
  v.resize(V_MAX);
  
  double s = density/rs;
  double t = ts/temperature;
  double ls = log(s);
  double lt = log(t);
  
  v[V_P] = c[0].n;
  v[V_H] = c[0].n;
  v[V_E] = 0.;
  v[V_S] = -c[0].n*log(s);
  double cp1 = 0.;
  double cp2 = c[0].n;
  double cp3 = c[0].n;
  v[V_F] = c[0].n*log(s);
  v[V_G] = c[0].n*(log(s)+1.);
  for (int i=1;i<40;i++) {
    double e = exp(c[i].I*ls+c[i].J*lt);
    v[V_P] += c[i].n*c[i].I*e;
    v[V_H] += c[i].n*(c[i].I+c[i].J)*e;
    v[V_E] += c[i].n*c[i].J*e;
    v[V_S] += c[i].n*(c[i].J-1)*e;
    cp1 -= c[i].n*c[i].J*(c[i].J-1)*e;
    cp2 += c[i].n*c[i].I*(1-c[i].J)*e;
    cp3 += c[i].n*c[i].I*(1+c[i].I)*e;
    v[V_F] += c[i].n*e;
    v[V_G] += c[i].n*(1+c[i].I)*e;
  }
  
  v[V_P] *= density*temperature*R;
  v[V_H] *= temperature*R;
  v[V_E] *= temperature*R;
  v[V_S] *= R;
  
  v[V_CV] = R*cp1;
  v[V_CP] = R*(cp1+cp2*cp2/cp3);
  v[V_KT] = cp3*R*temperature*density;
  
  // convert to m^2/s^2 from kJ/kg
  v[V_CS] = R*temperature*(cp3+cp2*cp2/cp1)*1000.;
  //if (v[V_CS] < 0.) throw std::runtime_error("Region3::getvars: sound speed imaginary");
  v[V_CS] = sqrt(v[V_CS]);

  v[V_F] *= temperature*R;
  v[V_G] *= temperature*R;

  v[V_R] = density;
  v[V_V] = 1./density;
  v[V_T] = temperature;

}

void Region5::getvars( const double pressure,
                       const double temperature,
                       std::vector<double> & v )
{
  v.resize(V_MAX);

  double p = pressure/ps;
  double t = ts/temperature;
  double lp = log(p);
  double lt = log(t);
  //std::cout << "pressure " << pressure << " temperature " << temperature << " p " << p << " t " << t << " lp " << lp << " lt " << lt << std::endl;
  v[V_V] = 1.;
  v[V_H] = 0.;
  v[V_E] = -1.;
  v[V_S] = -lp;
  v[V_CP] = 0.;
  double w1 = 1.;
  double w2 = -1.;
  v[V_G] = lp;
  v[V_F] = lp-1.;
  for (int i=0;i<12;i++) {
    //std::cout << "n " << i << " i " << c[i].I << " j " << c[i].J << " n " << c[i].n << std::endl;
    double e = exp(c[i].I*lp+c[i].J*lt);
    v[V_V] += c[i].n*c[i].I*e;
    v[V_H] += c[i].n*c[i].J*e;
    v[V_E] += c[i].n*(c[i].J-c[i].I)*e;
    v[V_S] += c[i].n*(c[i].J-1)*e;
    v[V_CP] -= c[i].n*c[i].J*(c[i].J-1)*e;
    w1 += c[i].n*c[i].I*(1-c[i].J)*e;
    w2 += c[i].n*c[i].I*(c[i].I-1)*e;
    v[V_G] += c[i].n*e;
    v[V_F] += c[i].n*(1-c[i].I)*e;
  }
  
  // convert to m^2/s^2 from kJ/kg
  v[V_CS] = -R*temperature*v[V_V]*v[V_V]/(w1*w1/v[V_CP]+w2)*1000.;
  //if (v[V_CS] < 0.) throw std::runtime_error("Region5::getvars: sound speed imaginary");
  v[V_CS] = sqrt(v[V_CS]);
  
  v[V_CV] = R*(v[V_CP]+w1*w1/w2);
  v[V_KT] = -pressure*v[V_V]/w2;

  v[V_R] = pressure/temperature/R/v[V_V];
  v[V_V] *= temperature*R/pressure;
  v[V_H] *= temperature*R;
  v[V_E] *= temperature*R;
  v[V_S] *= R;
  v[V_CP] *= R;
  v[V_G] *= temperature*R;
  v[V_F] *= temperature*R;

  v[V_P] = pressure;
  v[V_T] = temperature;
  v[V_G2P] = w2*R*temperature/pressure/pressure;
}

const std::string IAPWSIF97::name = "IAPWSIF97";

const EOSParamInfo IAPWSIF97::ParamList[] = {
  { "TEMP_EXTRAP", "IParam", TEMP_EXTRAP }
};
const int IAPWSIF97::ParamListSize = sizeof(ParamList)/sizeof(EOSParamInfo);

EOSModel * IAPWSIF97::createModel( const Parameters & params,
                                   const std::string name,
                                   const EOSData & /*init_data*/ )
{
  return new IAPWSIF97(params,name);
}

void IAPWSIF97::setRefStateCoef()
{
  krcs = 22064.0;
  krc.resize(5);
  for (int i=0;i<5;i++) krc[i].resize(6);
  krc[0][0] =  6.53786807199516;
  krc[0][1] = -5.61149954923348;
  krc[0][2] =  3.39624167361325;
  krc[0][3] = -2.27492629730878;
  krc[0][4] = 10.2631854662709 ;
  krc[0][5] =  1.97815050331519;
  krc[1][0] =  6.52717759281799;
  krc[1][1] = -6.30816983387575;
  krc[1][2] =  8.08379285492595;
  krc[1][3] = -9.82240510197603;
  krc[1][4] = 12.1358413791395 ;
  krc[1][5] = -5.54349664571295;
  krc[2][0] =  5.35500529896124;
  krc[2][1] = -3.96415689925446;
  krc[2][2] =  8.91990208918795;
  krc[2][3] = -12.0338729505790;
  krc[2][4] =  9.19494865194302;
  krc[2][5] = -2.16866274479712;
  krc[3][0] =   1.55225959906681 ;
  krc[3][1] =   0.464621290821181;
  krc[3][2] =   8.93237374861479 ;
  krc[3][3] = -11.0321960061126  ;
  krc[3][4] =   6.16780999933360 ;
  krc[3][5] =  -0.965458722086812;
  krc[4][0] =   1.11999926419994 ;
  krc[4][1] =   0.595748562571649;
  krc[4][2] =   9.88952565078920 ;
  krc[4][3] = -10.3255051147040  ;
  krc[4][4] =   4.66861294457414 ;
  krc[4][5] =  -0.503243546373828;
}

IAPWSIF97::IAPWSIF97( const Parameters & params,
                      const std::string name )
  : EOSModel("IAPWSIF97"),r1(), r2(), r3(), r4(), r5(), r23b(), Textrap(false), krcs(0.)
{
  //
  // find parameters for the model
  //
  const EOSParam * myparams = getParams(params,name);

  //
  // set extrapolation flag
  //
  if (myparams->params_int[TEMP_EXTRAP] != 0) Textrap = true;

  //
  // set transport property coefficients
  //
  setRefStateCoef();
}

IAPWSIF97::IAPWSIF97( )
  : EOSModel("IAPWSIF97"),r1(), r2(), r3(), r4(), r5(), r23b(), Textrap(false), krcs(0.)
{
  //
  // set transport property coefficients
  //
  setRefStateCoef();
}

IAPWSIF97::~IAPWSIF97()
{

}

void IAPWSIF97::clip_derivatives( std::vector<double> & v)
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

void IAPWSIF97::convert_to_SI( std::vector<double> & v )
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

double IAPWSIF97::computeRefStateK( const double r )
{
  int j;
  if (r <= 0.310559006) j = 0;
  else if (r <= 0.776397516) j = 1;
  else if (r <= 1.242236025) j = 2;
  else if (r <= 1.863354037) j = 3;
  else j = 4;

  return krcs*r*(krc[j][0]+r*(krc[j][1]+r*(krc[j][2]+r*(krc[j][3]+r*(krc[j][4]+r*krc[j][5])))));
}

void IAPWSIF97::computeTransportProperties( std::vector<double> & v )
{
  //
  // compute the reference state -- this temperature always lies in
  // region 2, but don't use it because the functional form may become
  // invalid before the desired state is reached.
  //
  // double tr = 1.5*647.096;

  double ktr = computeRefStateK( v[V_R]/322.0 );
  
  //
  // Calculate transport properties
  //
  v[V_M] = visc.getViscosity(v[V_R],v[V_T],v[V_KT],ktr);
  v[V_L] = thcon.getThermalConductivity(v[V_R],v[V_T],v[V_KT],ktr,v[V_CV],v[V_CP],v[V_M]);
}

void IAPWSIF97::getThermalEOS_PT_D( EOSData & data, const int phase )
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

void IAPWSIF97::getThermalEOS_PH_D( EOSData & data, const int phase )
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

void IAPWSIF97::getThermalEOS_PS_D( EOSData & data, const int phase )
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

void IAPWSIF97::getThermalEOS_PE_D( EOSData & data, const int phase )
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

void IAPWSIF97::getThermalEOS_RT_D( EOSData & data, const int phase )
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
void IAPWSIF97::calculatePhaseBoundaryInfo(const double lowP,const double highP)
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


void IAPWSIF97::getPhaseCoexistence( const int phase, const double pressure, double & temperature, double * enthalpy )
{
  if (phase != PH_LGMIX) throw std::runtime_error("IAPWSIF97::getPhaseCoexistence: invalid phase");

  //
  // convert to internal units
  //
  double ipressure = pressure/1000.;
  
  if (ipressure < 0.611657 || ipressure > 22064.) {
    throw std::runtime_error("IAPWSIF97::getPhaseCoexistence: pressure out of valid range");
  }
  else if (ipressure < r23b.getP(623.15)) {
    temperature = getR12satT(ipressure);
    std::vector<double> v1,v2;
    r1.getvars(ipressure,temperature,v1);
    r2.getvars(ipressure,temperature,v2);
    convert_to_SI(v1);
    convert_to_SI(v2);
    enthalpy[0] = v2[V_H];
    enthalpy[1] = v1[V_H];
  }
  else {
    temperature = getR3satT(ipressure);
    double rho1,rho2;
    getR3maxwellR(temperature,rho1,rho2);
    std::vector<double> v1,v2;
    r3.getvars(rho1,temperature,v1);
    r3.getvars(rho2,temperature,v2);
    convert_to_SI(v1);
    convert_to_SI(v2);
    enthalpy[0] = v1[V_H];
    enthalpy[1] = v2[V_H];
  }
}

double IAPWSIF97::getR3satT( const double pressure )
{
  //
  // check for valid pressure
  //
  if (pressure < r23b.getP(623.15) || pressure > 22064.*(1.+1.e-9)) {
    //std::cerr << " pressure " << pressure << " pbl " << r23b.getP(623.15) << " pbh " << 22064. << std::endl;
    throw std::runtime_error("IAPWFSIF97::getR3satT: pressure out of bounds");
  }

  //
  // temperature estimate from region 4
  //
  double Tsat = r4.getT(pressure);

  //
  // temperature bounds
  //
  double tmin = std::max(622.,Tsat-1.);
  double tmax = std::min(647.096,Tsat+1.);

  //
  // get temperature estimate
  //
  double lt = tmin;

  //
  // compute construction pressure
  //
  std::vector<double> v;
  double rho1,rho2;
  getR3maxwellR(lt,rho1,rho2);
  r3.getvars(rho1,lt,v);
  double lp = v[V_P];

  //
  // perturbed estimate
  //
  double llt = tmax;
  getR3maxwellR(llt,rho1,rho2);
  r3.getvars(rho1,llt,v);
  double llp = v[V_P];
  
  //
  // secant/bisection search for T
  //
  while (fabs(tmax-tmin) > 1.e-14*fabs(tmax+tmin)/2.) {
    double t = lt+(pressure-lp)*(lt-llt)/(lp-llp);
    if (t >= tmax || t <= tmin || lp-llp == 0.) {
      t = (tmax+tmin)/2.;
    }
    getR3maxwellR(t,rho1,rho2);
    r3.getvars(rho1,t,v);
    double p = v[V_P];
    
    if (p < pressure) tmin = t;
    else tmax = t;
    llt = lt;
    lt = t;
    llp = lp;
    lp = p;
  }

  return tmax;
}

double IAPWSIF97::getR12satT( const double pressure )
{
  //
  // check for valid pressure
  //
  if (pressure > r23b.getP(623.15)*1.001 || pressure < 0.611657*0.999) {
    throw std::runtime_error("IAPWFSIF97::getR12satT: pressure out of bounds");
  }

  //
  // temperature estimate from region 4
  //
  double Tsat = r4.getT(pressure);
  
  //
  // temperature bounds
  //
  double tmin = std::max(273.13,Tsat-1.);
  double tmax = std::min(623.25,Tsat+1.);

  //
  // get temperature estimate
  //
  double lt = tmin;

  //
  // calculate Gibb's free energy difference
  //
  std::vector<double> v1,v2;

  r1.getvars(pressure,lt,v1);
  r2.getvars(pressure,lt,v2);
  double ldG = v1[V_G]-v2[V_G];

  //
  // second estimate
  //
  double llt = tmax;
  r1.getvars(pressure,llt,v1);
  r2.getvars(pressure,llt,v2);
  double lldG = v1[V_G]-v2[V_G];

  //
  // secant/bisection search for T
  //
  while (fabs(tmax-tmin) > 1.e-14*fabs(tmax+tmin)/2.) {
    double t = lt-ldG*(lt-llt)/(ldG-lldG);
    if (t >= tmax || t <= tmin || ldG-lldG == 0.) {
      t = (tmax+tmin)/2.;
    }
    r1.getvars(pressure,t,v1);
    r2.getvars(pressure,t,v2);
    double dG = v1[V_G]-v2[V_G];
    
    if (dG < 0.) tmin = t;
    else tmax = t;
    llt = lt;
    lt = t;
    lldG = ldG;
    ldG = dG;
  }

  return tmax;
}

void IAPWSIF97::getR3maxwellR( const double temperature,
                               double & rho1,
                               double & rho2 )
{
  //
  // check temperature
  //
  if (temperature > 647.096) throw std::runtime_error("IAPWSIF97::getR3maxwellR: temperature must be below Tc");
  if (temperature < 622.) throw std::runtime_error("IAPWSIF97::getR3maxwellR: below min temperature");

  //
  // fit saturation density near critical point to rho=rhoc+/-rhoa*(tc-t)^0.5
  //
  // t in [647.09575:647.09595] gives a=93.379
  //
  // relative error for both densities about 2e-6 at the following temperature
  //
  if (temperature > 647.09585) {
    rho1 = 322.-93.379*sqrt(647.096-temperature);
    rho2 = 322.+93.379*sqrt(647.096-temperature);
    return;
  }

  //
  // deterimine spinodal densities
  //
  Region3Spinodal r3s(0.,1,r3);
  r3s.ResetT(temperature);
  r3s.ResetDir(1);
  double r2min = r3s.minimize(322.,630.,std::numeric_limits<double>::epsilon());
  r3s.ResetDir(-1);
  double r1max = r3s.minimize(30.,322.,std::numeric_limits<double>::epsilon());

  //
  // get pressure bounds from spinodal pressures
  //
  std::vector<double> v1,v2;

  r3.getvars(r1max,temperature,v1);
  r3.getvars(r2min,temperature,v2);
  double pmin = v2[V_P];
  double pmax = v1[V_P];

  //std::cout << " pmin " << pmin << " pmax " << pmax << " pest " << r4.getP(temperature) << std::endl;

  //
  // start at midpoint
  //
  double llp = (pmin+pmax)/2.;

  //
  // calculate Gibb's free energy difference
  //
  double rhol = r3.getRatPT(llp,temperature,30.,r1max,-1);
  double rhor = r3.getRatPT(llp,temperature,r2min,630.,1);
  r3.getvars(rhol,temperature,v1);
  r3.getvars(rhor,temperature,v2);
  double lldG = v1[V_G]-v2[V_G];

  //
  // second estimate
  //
  double lp = (4.*llp+pmin)/5.;
  rhol = r3.getRatPT(lp,temperature,30.,r1max,-1);
  rhor = r3.getRatPT(lp,temperature,r2min,630.,1);
  r3.getvars(rhol,temperature,v1);
  r3.getvars(rhor,temperature,v2);
  double ldG = v1[V_G]-v2[V_G];
  
  double dir = 1.;
  if (lldG < ldG) dir = -1.;

  //
  // update p bounds
  //
  if (ldG < 0.) {
    if (dir > 0.) {
      if (lldG < 0.) pmin = llp;
      else {
        pmin = lp;
        pmax = llp;
      }
    }
    else pmax = lp;
  }
  else {
    if (dir > 0.) pmax = lp;
    else {
      if (lldG > 0.) pmin = llp;
      else {
        pmin = lp;
        pmax = llp;
      }
    }
  }
  //std::cout << "pmin " << pmin << " pmax " << pmax << " ldG " << ldG << " lldG " << lldG << std::endl;

  //
  // secant/bisection search for p
  //
  while (fabs(pmax-pmin) > 1.e-14*fabs(pmax+pmin)/2.) {
    double p = lp-ldG*(lp-llp)/(ldG-lldG);
    if (p >= pmax || p <= pmin || ldG-lldG == 0.) {
      p = (pmax+pmin)/2.;
    }
    rhol = r3.getRatPT(p,temperature,30.,r1max,-1);
    rhor = r3.getRatPT(p,temperature,r2min,630.,1);
    r3.getvars(rhol,temperature,v1);
    r3.getvars(rhor,temperature,v2);
    //std::cout << "ptarget " << p << " p1 " << v1[0] << " p2 " << v2[0] << std::endl;
    double dG = v1[V_G]-v2[V_G];
    
    if (dG*dir < 0.) pmin = p;
    else pmax = p;
    llp = lp;
    lp = p;
    lldG = ldG;
    ldG = dG;
    //std::cout << "pmin " << pmin << " pmax " << pmax << " ldG " << ldG << " lldG " << lldG << std::endl;
  }

  rho1 = rhol;
  rho2 = rhor;

}

void IAPWSIF97::getvarsRT( const double density,
                           const double temperature,
                           const Phases phase,
                           std::vector<double> & v )
{
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  if (phase == PH_LGMIX) throw std::runtime_error("IAPWSIF97::getvarsRT: liquid-gas mixture not supported");

  r3.getvars(density,temperature,v);

}
  
void IAPWSIF97::getvarsPT( const double pressure,
                           const double temperature,
                           const Phases phase,
                           std::vector<double> & v )
{
  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  if (phase == PH_LGMIX) throw std::runtime_error("IAPWSIF97::getvarsPT: cannot evaluate liquid-gas mixture in P-T space");
  
  //
  // initial out of bounds pressure check
  //
  if (pressure > 100000. || pressure < 0.611657) {
    throw std::runtime_error("IAPWSIF97::getvarsPT: pressure outside range of validity for all phases");
  }

  if (pressure < r23b.getP(623.15)) {
    //std::cerr << " 125 regions " << std::endl;
    //
    // region 5
    //
    if (temperature > 1073.15) {
      r5.getvars(pressure,temperature,v);
    }
    else {
      //
      // check vapor curve
      //
      double ts = getR12satT(pressure);
      
      if (temperature < ts) r1.getvars(pressure,temperature,v);
      else r2.getvars(pressure,temperature,v);
    }
  }
  else {
    //std::cerr << "1325 regions " << std::endl;
    //
    // region 5
    //
    if (temperature > 1073.15) {
      //
      // region 5 out of bounds pressure
      //
      if (pressure > 50000.) {
        throw std::runtime_error("IAPWSIF97::getvarsPT: pressure outside range of validity for region 5");
      }

      r5.getvars(pressure,temperature,v);
    }
    else {
      //
      // check region 1-3-2 boundaries
      //
      if (temperature <= 623.15) r1.getvars(pressure,temperature,v);
      else if (temperature >= r23b.getT(pressure)) r2.getvars(pressure,temperature,v);
      else {
        //
        // in region 3
        //
        double minr(40.),maxr(800.);
        double rho;
        if (pressure < 22064.) {
          //
          // deal with liquid-vapor transition
          //
          double ts = getR3satT(pressure);
          if (temperature > ts) {
            rho = r3.getRatPT( pressure, temperature, minr, maxr, -1 );
          }
          else {
            minr = 322.;
            rho = r3.getRatPT( pressure, temperature, minr, maxr, 1 );
          }
        }
        else {
          //
          // above liquid-vapor transition so just get the state
          //
          rho = r3.getRatPT( pressure, temperature, minr, maxr, 0);
        }
        r3.getvars(rho,temperature,v);

      }
      
    }
  }

  //
  // all non-failure cases fall through to here for transport calculation
  //
  clip_derivatives(v);
  computeTransportProperties(v);

}

void IAPWSIF97::getvarsPX( const double pressure,
                           const double xvar,
                           const Phases phase,
                           const Vars xtype,
                           std::vector<double> & v )
{
  //std::cout << std::setprecision(18);
  //std::cout << "IAPWSIF97::getvarsPX: p " << pressure << " x " << xvar << " phase " << phase << " type " << xtype << std::endl;
  //std::cout << std::setprecision(6);

  v.resize(V_MAX);

  for (int i=0;i<V_MAX;i++) v[i] = 0.;

  //
  // check xtype validity
  //
  if (xtype != V_H && xtype != V_E && xtype != V_S) {
    throw std::runtime_error("IAPWSIF97::getvarsPX: invalid variable type");
  }

  //
  // working vector
  //
  std::vector<double> d;

  //
  // initial out of bounds pressure check
  //
  if (pressure > 100000.*(1.+1.e-10) || pressure < 0.611657) {
    //std::cout << " pressure " << pressure << " diff " << pressure-100000. << std::endl;
    throw std::runtime_error("IAPWSIF97::getvarsPX: pressure outside range of validity for all phases");
  }

  if (phase & (PH_GAS | PH_LIQUID | PH_FLUID)) {
    //std::cerr << " fluid phase " << phase << std::endl;
    //
    // regions 1,2,5
    //
    if (pressure < r23b.getP(623.15)) {
      //std::cerr << " 125 regions " << std::endl;
      //
      // check 2-5 boundary
      //
      //std::cerr << " 2-5 boundary check " << std::endl;
      r2.getvars(pressure,1073.15,d);
      //std::cerr << " 2-5 bound var " << d[xtype] << " target " << xvar << std::endl;
      if (d[xtype] < xvar) {
        //
        // give 50 K of extrapolation room
        //
        double maxT = 2273.15;
        if (Textrap) maxT += 50.;

        //
        // in region 5 if not out of bounds
        //
        r5.getvars(pressure,maxT,d);
        if (d[xtype] < xvar) {
          throw std::runtime_error("IAPWSIF97::getvarsPX: temperature outside range of validity for region 5");
        }

        double t = r5.getTatPX(pressure,xvar,xtype,1073.15,maxT);
        r5.getvars(pressure,t,v);
        //std::cout << " region 5 " << std::endl;
      }
      else {
        //
        // in region 1 or 2 so get 1-2 boundary
        //
        double ts = getR12satT(pressure);

        //
        // give 50 K of extrapolation room
        //
        double minT = 273.16;
        if (Textrap) minT -= 50.;

        if (phase == PH_LIQUID) {
          //
          // region 1 requested -- give a few degrees of meta-stable behavior
          //
          double t = r1.getTatPX(pressure,xvar,xtype,minT,ts+20.);
          r1.getvars(pressure,t,v);
          //std::cout << " region 1 " << std::endl;
        }
        else if (phase == PH_GAS) {
          //std::cout << "gas r2 calc " << std::endl;
          //
          // region 2 requested -- give a few degrees of meta-stable behavior
          //
          double t = r2.getTatPX(pressure,xvar,xtype,ts-50.,1073.15);
          r2.getvars(pressure,t,v);
          //std::cout << " region 2 " << std::endl;
        }
        else {
          //
          // PH_FLUID -- if xvar in meta-stable region choose closest side
          //
          r1.getvars(pressure,ts,d);
          std::vector<double> d2;
          r2.getvars(pressure,ts,d2);

          double f = (xvar-d2[xtype])/(d[xtype]-d2[xtype]);

          if (f < 0.5) {
            //
            // in or closer to region 2
            //
            double t = r2.getTatPX(pressure,xvar,xtype,ts-30.,1073.15);
            r2.getvars(pressure,t,v);
            //std::cout << " region 2 " << std::endl;
          }
          else {
            //
            // in or closer to region 1
            //
            double t = r1.getTatPX(pressure,xvar,xtype,minT,ts+20.);
            r1.getvars(pressure,t,v);
            //std::cout << " region 1 " << std::endl;
          }
        }
      }
    }
    //
    // regions 1,3,2,5
    //
    else {
      //std::cerr << " 1325 regions " << std::endl;
      //
      // check X var at region 1-3 boundary
      //
      r1.getvars(pressure,623.15,d);
      if (d[xtype] >= xvar) {
        //
        // in region 1
        //

        //
        // give 50 K of extrapolation room
        //
        double minT = 273.16;
        if (Textrap) minT -= 50.;
        double t = r1.getTatPX(pressure,xvar,xtype,minT,623.15);
        r1.getvars(pressure,t,v);
        //std::cout << " region 1 " << std::endl;
      }
      else {
        //
        // check X var at region 3-2 boundary
        //
        r2.getvars(pressure,r23b.getT(pressure),d);
        //std::cout << " xvar " << xvar << " r2r23b " << d[xtype] << " r23bt " << r23b.getT(pressure) << std::endl;
        if (d[xtype] > xvar) {
          //std::cerr << " region 3 " << std::endl;
          //
          // in region 3
          //
          double rho,temp;
          r3.getRTatPX( pressure, xvar, xtype, 322., 657.,rho,temp);
          r3.getvars(rho,temp,v);
        }
        else {
          //
          // check X var at region 2-5 boundary
          //
          r2.getvars(pressure,1073.15,d);
          //std::cout << " xvar " << xvar << " r21073 " << d[xtype] << std::endl;
          if (d[xtype] >= xvar) {
            //
            // in region 2
            //
            double t = r2.getTatPX(pressure,xvar,xtype,r23b.getT(pressure),1073.15);
            r2.getvars(pressure,t,v);
            //std::cout << " region 2 " << std::endl;
          }
          else {
            //
            // region 5 out of bounds pressure
            //
            if (pressure > 50000.) {
              //
              // try extrapolating region 2 but less than limit of 1080 K
              //
              if (Textrap) {
                r2.getvars(pressure,1079.999,d);
                //std::cout << "p " << pressure << " t " << 1123.15 << " xvar " << xvar << " calc " << d[xtype] << std::endl;
                if (d[xtype] >= xvar) {
                  //
                  // in extrapolated region 2
                  //
                  double t = r2.getTatPX(pressure,xvar,xtype,r23b.getT(pressure),1079.999);
                  r2.getvars(pressure,t,v);
                  //std::cout << " region 2 " << std::endl;
                }
                else {
                  throw std::runtime_error("IAPWSIF97::getvarsPX: pressure outside of extrapolation range for region 2");
                }
              }
              else {
                throw std::runtime_error("IAPWSIF97::getvarsPX: pressure outside range of validity for region 5");
              }
            }
            else {
              //
              // check X var at region 5 high T boundary with 50 K of extrapolation room
              //
              double maxT = 2273.15;
              if (Textrap) maxT += 50.;

              r5.getvars(pressure,maxT,d);
              if (d[xtype] < xvar) {
                throw std::runtime_error("IAPWSIF97::getvarsPX: temperature outside range of validity for region 5");
              }

              //
              // allow some extrapolation at low T boundary
              //
              double minT = 1073.15;
              if (Textrap) minT -= 10.;
              double t = r5.getTatPX(pressure,xvar,xtype,minT,maxT);
              r5.getvars(pressure,t,v);
              //std::cout << " region 5 " << std::endl;
            }
          }
        }
      }
    }

    //
    // all non-failure cases fall through to here for transport calculation
    //
    clip_derivatives(v);
    computeTransportProperties(v);
  }
  else if (phase == PH_LGMIX) {
    if (pressure > 22064.*(1.+1.e-9)) {
      throw std::runtime_error("IAPWSIF97::getvarsPX: pressure too high for liquid-gas mixture");
    }
    else if (pressure < r23b.getP(623.15)) {
      //
      // lg boundary between region 1 and 2
      //
      getR12maxwellPX( pressure, getR12satT(pressure), xvar, xtype, v );

    }
    else {
      //
      // lg boundary in region 3
      //
      getR3maxwellPX( pressure, getR3satT(pressure), xvar, xtype, v );

    }

    //
    // note that transport calculation happen inside maxwell routines instead of here
    //
  }
  else {
    throw std::runtime_error("IAPWSIF97::getvarsPX: invalid phase");
  }
}

void IAPWSIF97::getR12maxwellDerivs( const double pressure,
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
  // make pressure grid
  //
  std::vector<double> p(5,0.);

  for (int i=0;i<5;i++) p[i] = pressure*(1.-0.0001*(i-2));

  //
  // compute maxwell temperatures
  //
  std::vector<double> t(5,0.);

  for (int i=0;i<5;i++) t[i] = getR12satT(p[i]);

  //
  // compute states on boundaries
  //
  std::vector<std::vector<double> > v1(5),v2(5);

  for (int i=0;i<5;i++) {
    r1.getvars(p[i],t[i],v1[i]);
    r2.getvars(p[i],t[i],v2[i]);
  }

  //
  // compute derivatives
  //
  double dpdt = (v1[2][V_S]-v2[2][V_S])/(v1[2][V_V]-v2[2][V_V]);
  std::vector<double> f(5);
  for (int i=0;i<V_MAX;i++) {
    for (int j=0;j<5;j++) f[j] = v1[j][i];
    vd1[i] = FivePointStencilDeriv(f,0.0001*pressure,2,1)*dpdt;    
    for (int j=0;j<5;j++) f[j] = v2[j][i];
    vd2[i] = FivePointStencilDeriv(f,0.0001*pressure,2,1)*dpdt;    
  }

}

void IAPWSIF97::getR12maxwellPX( const double pressure,
                                 const double temperature,
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
    throw std::runtime_error("IAPWSIF97::getR12maxwellPX: invalid variable type");
  }

  //
  // Compute variables at Maxwell construction boundaries
  //
  std::vector<double> v1,v2;

  r1.getvars(pressure,temperature,v1);
  r2.getvars(pressure,temperature,v2);
  
  //
  // get transport properties at constuction boundaries
  //
  clip_derivatives(v1);
  clip_derivatives(v2);
  computeTransportProperties(v1);
  computeTransportProperties(v2);

  //
  // Compute lever rule for given variable
  //
  double f = (xvar-v2[xtype])/(v1[xtype]-v2[xtype]);

  //
  // Compute variables at f
  //
  v[V_V] = v1[V_V]*f+v2[V_V]*(1.-f);
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
  getR12maxwellDerivs( pressure, vd1, vd2 );

  double dfdt = (-vd2[V_V]*(v1[V_V]-v2[V_V])-(v[V_V]-v2[V_V])*(vd1[V_V]-vd2[V_V]))/(v1[V_V]-v2[V_V])/(v1[V_V]-v2[V_V]);
  v[V_CV] = vd1[V_E]*f+vd2[V_E]*(1-f)+(v1[V_E]-v2[V_E])*dfdt;
  v[V_CS] = sqrt(temperature*v[V_V]*v[V_V]*vd1[V_P]*vd1[V_P]/v[V_CV]);

  //
  // overwrite xtype var
  //
  v[xtype] = xvar;
}

void IAPWSIF97::getR3maxwellDerivs( const double temperature,
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
  std::vector<double> t(5,0.);

  //
  // find location in stencil
  //
  int lt = 0;
  double dt = 0.0001;
  if (temperature > 647.09) dt = 0.000001;
  else if (temperature > 647.) dt = exp(log(0.0001)+(log(0.000001)-log(0.0001))*(temperature-647.)/0.09);
  if (temperature*(1.+2.*dt) < 647.096) lt = 2;
  else if (temperature*(1.+dt) < 647.096) lt = 1;
  lt=0;

  for (int i=0;i<5;i++) t[i] = temperature*(1.-dt*(i-lt));

  //
  // compute maxwell densities
  //
  std::vector<double> r1(5),r2(5);

  for (int i=0;i<5;i++) getR3maxwellR(t[i],r1[i],r2[i]);

  //
  // compute states on boundaries
  //
  std::vector<std::vector<double> > v1(5),v2(5);

  for (int i=0;i<5;i++) {
    r3.getvars(r1[i],t[i],v1[i]);
    r3.getvars(r2[i],t[i],v2[i]);
  }

  //
  // compute derivatives
  //
  std::vector<double> f(5);
  for (int i=0;i<V_MAX;i++) {
    for (int j=0;j<5;j++) f[j] = v1[j][i];
    vd1[i] = FivePointStencilDeriv(f,dt*temperature,lt,1);    
    for (int j=0;j<5;j++) f[j] = v2[j][i];
    vd2[i] = FivePointStencilDeriv(f,dt*temperature,lt,1);
  }

  
}

void IAPWSIF97::getR3maxwellPX( const double pressure,
                                const double temperature,
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
    throw std::runtime_error("IAPWSIF97::getR3maxwellPX: invalid variable type");
  }

  //
  // Compute variables at Maxwell construction boundaries
  //
  std::vector<double> v1,v2;
  double rho1,rho2;

  getR3maxwellR( temperature, rho1, rho2 );
  r3.getvars(rho1,temperature,v1);
  r3.getvars(rho2,temperature,v2);

  //
  // get transport properties at constuction boundaries
  //
  clip_derivatives(v1);
  clip_derivatives(v2);
  computeTransportProperties(v1);
  computeTransportProperties(v2);
  if (isnan(v1[V_M]) || isnan(v2[V_M])) {
    throw std::runtime_error("IAPWS-IF97::getR3maxwellPX: transport property is nan");
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

  //std::cerr << " v " << v[V_V] << " v1 " << 1./(322.-93.379*sqrt(647.096-temperature)) << " v2 " << 1./(322.+93.379*sqrt(647.096-temperature)) << std::endl;

  std::vector<double> vd1,vd2;
  getR3maxwellDerivs( temperature, vd1, vd2 );

  //std::cerr << " vd1 " << vd1[V_V] << " " << -93.379/2./sqrt(647.096-temperature)/(322.-93.379*sqrt(647.096-temperature))/(322.-93.379*sqrt(647.096-temperature)) << std::endl;
  //std::cerr << " vd2 " << vd2[V_V] << " " << 93.379/2./sqrt(647.096-temperature)/(322.+93.379*sqrt(647.096-temperature))/(322.+93.379*sqrt(647.096-temperature)) << std::endl;

  double dfdt = (-vd2[V_V]*(v1[V_V]-v2[V_V])-(v[V_V]-v2[V_V])*(vd1[V_V]-vd2[V_V]))/(v1[V_V]-v2[V_V])/(v1[V_V]-v2[V_V]);
  //std::cerr << " dfdt " << dfdt;
  //dfdt = (v[V_V]*(322.*322+93.379*93.379*(647.096-temperature))-322.)/4./93.379/sqrt(647.096-temperature)/(647.096-temperature);
  //std::cerr << " dfdt2 " << dfdt << std::endl;
  v[V_CV] = vd1[V_E]*f+vd2[V_E]*(1-f)+(v1[V_E]-v2[V_E])*dfdt;
  v[V_CS] = sqrt(temperature*v[V_V]*v[V_V]*vd1[V_P]*vd1[V_P]/v[V_CV]);

  /*
  double rc = 322.;
  double ra = 93.379;
  double rb = 4.1;
  double ec = 2019.0251;
  double ea = 136.498;
  double eb = 25.6;

  double cv1 = ea*ra*v[V_V];
  double cv2 = (-ra*eb-ea*rb+ea*v[V_V]*(ra*ra+2.*rb*(rc+rb*(temperature-647.096))))/ra;
  */

  double cvc = 19.95506;
  double cva = 6.3867;
  double cvb = -23.4492;
  if (temperature > 647.0948) {
    double cve1 = cvc+cva*sqrt(647.096-temperature)+cvb*(647.096-temperature);
    double cve2 = cvc-cva*sqrt(647.096-temperature)+cvb*(647.096-temperature);
    v[V_CV] = cve1*f+cve2*(1.-f);    
    v[V_CS] = sqrt(temperature*v[V_V]*v[V_V]*vd1[V_P]*vd1[V_P]/v[V_CV]);
  }
  //std::cerr << std::setprecision(15);
  //std::cerr << " t " << temperature << " cvderiv " << v[V_CV] << " cv1 " << cv1 << " cv2 " << cv2 << " cs " << v[V_CS] << std::endl;


  //
  // overwrite xtype var
  //
  v[xtype] = xvar;
}
