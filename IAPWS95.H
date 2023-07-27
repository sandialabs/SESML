/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef IAPWS95_H
#define IAPWS95_H

#include <cmath>
#include <stdexcept>

#include "Region.H"
#include "DFMin.H"
#include "Vars.H"
#include "IAPWS-IF97.H"
#include "IAPWStransport.H"

struct CoefficientData
{
  //
  // convenience constructor
  //
  CoefficientData( const int ci,
                   const int di,
                   const double ti,
                   const double ni )
    : c(ci), d(di), t(ti), n(ni) {}

  //
  // null constructor
  //
  CoefficientData()
    : c(0), d(0), t(0.0), n(0.0) {}

  //
  // two integer coefficients
  //
  int c;
  int d;

  //
  // two double coefficients
  //
  double t;
  double n;

};

class RegionRTSpinodal: public DFMin
{
public:
  //
  // constructor
  //
  RegionRTSpinodal( RegionRT * pregion,
                    const double temperature,
                    const int direction )
    : rrt(pregion), t(temperature), dir(direction)
  {
    if (dir != 1 && dir != -1) throw std::runtime_error("RegionRTSpinodal: bad direction passed to constructor");
    if (rrt == 0) throw std::runtime_error("RegionRTSpinodal: invalid region pointer");
  }

  //
  // destructor
  //
  ~RegionRTSpinodal() {}

  void ResetT( const double temperature ) {
    t = temperature;
  }

  void ResetDir( const int direction ) {
    if (direction != 1 && direction != -1) throw std::runtime_error("RegionRTSpinodal::ResetDir: bad direction input");
    dir = direction;
  }

protected:
  
  double f( const double rho ) {
    std::vector<double> v;
    rrt->getvars(rho,t,v);
    return dir*v[V_P];
  }

private:

  RegionRT * rrt;
  double t;
  int dir;

};

class IAPWS95 : public RegionRT, public EOSModel
{
public:

  //
  // start of stuff for EOSModel support
  //
  enum ParamDouble {
    TOT_PARAMDOUBLE
  };

  enum ParamInt {
    TEMP_EXTRAP,
    TOT_PARAMINT
  };

  enum ParamString {
    TOT_PARAMSTRING
  };

  //
  // model name string
  //
  static const std::string name;

  //
  // parameter parsing information
  //
  static const EOSParamInfo ParamList[];
  static const int ParamListSize;

  //
  // eos model factory method
  //
  static EOSModel * createModel( const Parameters & params,
                                 const std::string name,
                                 const EOSData & init_data );

  //
  // thermal EOS methods with derivatives
  //
  void getThermalEOS_PT_D( EOSData & data, const int phase );
  void getThermalEOS_PH_D( EOSData & data, const int phase );
  void getThermalEOS_PS_D( EOSData & data, const int phase );
  void getThermalEOS_PE_D( EOSData & data, const int phase );
  void getThermalEOS_RT_D( EOSData & data, const int phase );

  //
  // Method to request calculation of the phase boundary information for a model
  //
  // Temperature bounds give the range of temperatures within which
  // boundary calculations are desired.
  //
  void calculatePhaseBoundaryInfo(const double /*lowT*/,const double /*highT*/);

  //
  // given the pressure, get the phase coexistence temperature and enthalpy
  //
  void getPhaseCoexistence( const int phase, const double pressure, double & temperature, double * enthalpy );

  //
  // end of stuff for EOSModel support
  //

  enum Phases {
    PH_NONE = 0x0,
    PH_GAS = 0x1,
    PH_LIQUID = 0x2,
    PH_FLUID = PH_GAS | PH_LIQUID,
    PH_LGMIX = 0x4
  };

  //
  // constructors
  //
  IAPWS95( const Parameters & params,
           const std::string name );
  IAPWS95( );

  //
  // destructor
  //
  ~IAPWS95() {}

  //
  // evaluate the Helmholtz free energy
  //
  void getvars( const double density,
                const double temperature,
                std::vector<double> & v );

  void getvars2( const double density,
                 const double temperature,
                 const Vars xtype,
                 std::vector<double> & v );

  void getvarsPT( const double pressure,
                  const double temperature,
                  const Phases phase,
                  std::vector<double> & v );

  void getvarsPX( const double pressure,
                  const double xvar,
                  const Phases phase,
                  const Vars xtype,
                  std::vector<double> & v );

  void getvarsRT( const double density,
                  const double temperature,
                  const Phases phase,
                  std::vector<double> & v );

  //
  // get the saturation temperature at a given pressure
  //
  double getSatT( const double pressure );

  void getSatProps( const double pressure,
                    double & density1,
                    double & density2,
                    double & temperature );

  //
  // get the state on the Maxwell construction for the given pressure
  // and xvar. temperature must be from a call to getSatT to ensure
  // one is actually on the construction.
  //
  void getMaxwellPX( const double pressure,
                     //const double temperature,
                     const double xvar,
                     const Vars xtype,
                     std::vector<double> & v );

  //
  // get the saturation densities on the Maxwell construction
  //
  void getMaxwellR( const double temperature,
                    double & rho1,
                    double & rho2 );

  void getMaxwellDerivs( const double temperature,
                         std::vector<double> & vd1,
                         std::vector<double> & vd2 );

  void computeTransportProperties( std::vector<double> & v );

  double getSpinodalDensityBound( const double temperature,
                                  const bool glstate );

private:

  bool getSatPropsNearCritical( const double pressure,
                                double & density1,
                                double & density2,
                                double & temperature );

  void convert_to_SI( std::vector<double> & v );

  void setCoefficients() {

    // ideal gas part                    gamma_i      n_i            
    ci.push_back(CoefficientData( 0, 0,  0.0       , -8.3204464837497 ));
    ci.push_back(CoefficientData( 0, 0,  0.0       ,  6.6832105275932 ));
    ci.push_back(CoefficientData( 0, 0,  0.0       ,  3.00632         ));
    ci.push_back(CoefficientData( 0, 0,  1.28728967,  0.012436        ));
    ci.push_back(CoefficientData( 0, 0,  3.53734222,  0.97315         ));
    ci.push_back(CoefficientData( 0, 0,  7.74073708,  1.27950         ));
    ci.push_back(CoefficientData( 0, 0,  9.24437796,  0.96956         ));
    ci.push_back(CoefficientData( 0, 0, 27.5075105 ,  0.24873         ));

    // residualpart              c_i d_i  t_i           n_i
    cr.push_back(CoefficientData( 0,  1, -0.5  ,  0.12533547935523e-1  ));
    cr.push_back(CoefficientData( 0,  1,  0.875,  0.78957634722828e+1  ));
    cr.push_back(CoefficientData( 0,  1,  1    , -0.87803203303561e+1  ));
    cr.push_back(CoefficientData( 0,  2,  0.5  ,  0.31802509345418     ));
    cr.push_back(CoefficientData( 0,  2,  0.75 , -0.26145533859358     ));
    cr.push_back(CoefficientData( 0,  3,  0.375, -0.78199751687981e-2  ));
    cr.push_back(CoefficientData( 0,  4,  1.0  ,  0.88089493102134e-2  ));
    cr.push_back(CoefficientData( 1,  1,  4.0  , -0.66856572307965     ));
    cr.push_back(CoefficientData( 1,  1,  6.0  ,  0.20433810950965     ));
    cr.push_back(CoefficientData( 1,  1, 12.0  , -0.66212605039687e-4  ));
    cr.push_back(CoefficientData( 1,  2,  1.0  , -0.19232721156002     ));
    cr.push_back(CoefficientData( 1,  2,  5.0  , -0.25709043003438     ));
    cr.push_back(CoefficientData( 1,  3,  4.0  ,  0.16074868486251     ));
    cr.push_back(CoefficientData( 1,  4,  2.0  , -0.40092828925807e-1  ));
    cr.push_back(CoefficientData( 1,  4, 13.0  ,  0.39343422603254e-6  ));
    cr.push_back(CoefficientData( 1,  5,  9.0  , -0.75941377088144e-5  ));
    cr.push_back(CoefficientData( 1,  7,  3.0  ,  0.56250979351888e-3  ));
    cr.push_back(CoefficientData( 1,  9,  4.0  , -0.15608652257135e-4  ));
    cr.push_back(CoefficientData( 1, 10, 11.0  ,  0.11537996422951e-8  ));
    cr.push_back(CoefficientData( 1, 11,  4.0  ,  0.36582165144204e-6  ));
    cr.push_back(CoefficientData( 1, 13, 13.0  , -0.13251180074668e-11 ));
    cr.push_back(CoefficientData( 1, 15,  1.0  , -0.62639586912454e-9  ));
    cr.push_back(CoefficientData( 2,  1,  7.0  , -0.10793600908932     ));
    cr.push_back(CoefficientData( 2,  2,  1.0  ,  0.17611491008752e-1  ));
    cr.push_back(CoefficientData( 2,  2,  9.0  ,  0.22132295167546     ));
    cr.push_back(CoefficientData( 2,  2, 10.0  , -0.40247669763528     ));
    cr.push_back(CoefficientData( 2,  3, 10.0  ,  0.58083399985759     ));
    cr.push_back(CoefficientData( 2,  4,  3.0  ,  0.49969146990806e-2  ));
    cr.push_back(CoefficientData( 2,  4,  7.0  , -0.31358700712549e-1  ));
    cr.push_back(CoefficientData( 2,  4, 10.0  , -0.74315929710341     ));
    cr.push_back(CoefficientData( 2,  5, 10.0  ,  0.47807329915480     ));
    cr.push_back(CoefficientData( 2,  6,  6.0  ,  0.20527940895948e-1  ));
    cr.push_back(CoefficientData( 2,  6, 10.0  , -0.13636435110343     ));
    cr.push_back(CoefficientData( 2,  7, 10.0  ,  0.14180634400617e-1  ));
    cr.push_back(CoefficientData( 2,  9,  1.0  ,  0.83326504880713e-2  ));
    cr.push_back(CoefficientData( 2,  9,  2.0  , -0.29052336009585e-1  ));
    cr.push_back(CoefficientData( 2,  9,  3.0  ,  0.38615085574206e-1  ));
    cr.push_back(CoefficientData( 2,  9,  4.0  , -0.20393486513704e-1  ));
    cr.push_back(CoefficientData( 2,  9,  8.0  , -0.16554050063734e-2  ));
    cr.push_back(CoefficientData( 2, 10,  6.0  ,  0.19955571979541e-2  ));
    cr.push_back(CoefficientData( 2, 10,  9.0  ,  0.15870308324157e-3  ));
    cr.push_back(CoefficientData( 2, 12,  8.0  , -0.16388568342530e-4  ));
    cr.push_back(CoefficientData( 3,  3, 16.0  ,  0.43613615723811e-1  ));
    cr.push_back(CoefficientData( 3,  4, 22.0  ,  0.34994005463765e-1  ));
    cr.push_back(CoefficientData( 3,  4, 23.0  , -0.76788197844621e-1  ));
    cr.push_back(CoefficientData( 3,  5, 23.0  ,  0.22446277332006e-1  ));
    cr.push_back(CoefficientData( 4, 14, 10.0  , -0.62689710414685e-4  ));
    cr.push_back(CoefficientData( 6,  3, 50.0  , -0.55711118565645e-9  ));
    cr.push_back(CoefficientData( 6,  6, 44.0  , -0.19905718354408     ));
    cr.push_back(CoefficientData( 6,  6, 46.0  ,  0.31777497330738     ));
    cr.push_back(CoefficientData( 6,  6, 50.0  , -0.11841182425981     ));

    // residualpart             t_i beta_i gamma_i       n_i
    cr.push_back(CoefficientData( 0, 150, 1.21, -0.31306260323435e+2 ));
    cr.push_back(CoefficientData( 1, 150, 1.21,  0.31546140237781e+2 ));
    cr.push_back(CoefficientData( 4, 250, 1.25, -0.25213154341695e+4 ));

    // residualpart              C_i  D_i   b_i         n_i
    cr.push_back(CoefficientData( 28, 700, 0.85, -0.14874640856724 ));
    cr.push_back(CoefficientData( 32, 800, 0.95,  0.31806110878444 ));
  }

  //
  // coefficients for ideal-gas and residual parts
  //
  std::vector<CoefficientData> ci;
  std::vector<CoefficientData> cr;

  Region4 r4;
  IAPWSvisc2008 visc;
  IAPWSthcon2011 thcon;

  //
  // allow extrapolation in temperature
  //
  bool Textrap;

};

#endif
