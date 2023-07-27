/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef IAPWS_IF97_H
#define IAPWS_IF97_H

#include <cmath>
#include <stdexcept>

#include "Region.H"
#include "DFMin.H"
#include "Vars.H"
#include "EOSModel.H"
#include "EOSParam.H"
#include "EOSFactory.H"
#include "IAPWStransport.H"

class Region1 : public RegionPT
{
public:
  //
  // constructor
  //
  Region1()
    : RegionPT( 16530.0, 1386.0, 0.461526 ) {
    c.push_back(SeriesCoefficient(  0,  -2,  0.14632971213167     ));
    c.push_back(SeriesCoefficient(  0,  -1, -0.84548187169114     ));
    c.push_back(SeriesCoefficient(  0,   0, -0.37563603672040e1   ));
    c.push_back(SeriesCoefficient(  0,   1,  0.33855169168385e1   ));
    c.push_back(SeriesCoefficient(  0,   2, -0.95791963387872     ));
    c.push_back(SeriesCoefficient(  0,   3,  0.15772038513228     ));
    c.push_back(SeriesCoefficient(  0,   4, -0.16616417199501e-1  ));
    c.push_back(SeriesCoefficient(  0,   5,  0.81214629983568e-3  ));
    c.push_back(SeriesCoefficient(  1,  -9,  0.28319080123804e-3  ));
    c.push_back(SeriesCoefficient(  1,  -7, -0.60706301565874e-3  ));
    c.push_back(SeriesCoefficient(  1,  -1, -0.18990068218419e-1  ));
    c.push_back(SeriesCoefficient(  1,   0, -0.32529748770505e-1  ));
    c.push_back(SeriesCoefficient(  1,   1, -0.21841717175414e-1  ));
    c.push_back(SeriesCoefficient(  1,   3, -0.52838357969930e-4  ));
    c.push_back(SeriesCoefficient(  2,  -3, -0.47184321073267e-3  ));
    c.push_back(SeriesCoefficient(  2,   0, -0.30001780793026e-3  ));
    c.push_back(SeriesCoefficient(  2,   1,  0.47661393906987e-4  ));
    c.push_back(SeriesCoefficient(  2,   3, -0.44141845330846e-5  ));
    c.push_back(SeriesCoefficient(  2,  17, -0.72694996297594e-15 ));
    c.push_back(SeriesCoefficient(  3,  -4, -0.31679644845054e-4  ));
    c.push_back(SeriesCoefficient(  3,   0, -0.28270797985312e-5  ));
    c.push_back(SeriesCoefficient(  3,   6, -0.85205128120103e-9  ));
    c.push_back(SeriesCoefficient(  4,  -5, -0.22425281908000e-5  ));
    c.push_back(SeriesCoefficient(  4,  -2, -0.65171222895601e-6  ));
    c.push_back(SeriesCoefficient(  4,  10, -0.14341729937924e-12 ));
    c.push_back(SeriesCoefficient(  5,  -8, -0.40516996860117e-6  ));
    c.push_back(SeriesCoefficient(  8, -11, -0.12734301741641e-8  ));
    c.push_back(SeriesCoefficient(  8,  -6, -0.17424871230634e-9  ));
    c.push_back(SeriesCoefficient( 21, -29, -0.68762131295531e-18 ));
    c.push_back(SeriesCoefficient( 23, -31,  0.14478307828521e-19 ));
    c.push_back(SeriesCoefficient( 29, -38,  0.26335781662795e-22 ));
    c.push_back(SeriesCoefficient( 30, -39, -0.11947622640071e-22 ));
    c.push_back(SeriesCoefficient( 31, -40,  0.18228094581404e-23 ));
    c.push_back(SeriesCoefficient( 32, -41, -0.93537087292458e-25 ));
  }

  //
  // destructor
  //
  ~Region1() {}

  //
  // evaluate the region
  //
  void getvars( const double pressure,
                const double temperature,
                std::vector<double> & v );

};

class Region2 : public RegionPT
{
public:
  //
  // constructor
  //
  Region2()
    : RegionPT( 1000.0, 540.0, 0.461526 ) {
    c.push_back(SeriesCoefficient( 0,  0, -0.96927686500217e1  ));
    c.push_back(SeriesCoefficient( 0,  1,  0.10086655968018e2  ));
    c.push_back(SeriesCoefficient( 0, -5, -0.56087911283020e-2 ));
    c.push_back(SeriesCoefficient( 0, -4,  0.71452738081455e-1 ));
    c.push_back(SeriesCoefficient( 0, -3, -0.40710498223928    ));
    c.push_back(SeriesCoefficient( 0, -2,  0.14240819171444e1  ));
    c.push_back(SeriesCoefficient( 0, -1, -0.43839511319450e1  ));
    c.push_back(SeriesCoefficient( 0,  2, -0.28408632460772    ));
    c.push_back(SeriesCoefficient( 0,  3,  0.21268463753307e-1 ));
    c.push_back(SeriesCoefficient(  1,  0, -0.17731742473213e-2  ));
    c.push_back(SeriesCoefficient(  1,  1, -0.17834862292358e-1  ));
    c.push_back(SeriesCoefficient(  1,  2, -0.45996013696365e-1  ));
    c.push_back(SeriesCoefficient(  1,  3, -0.57581259083432e-1  ));
    c.push_back(SeriesCoefficient(  1,  6, -0.50325278727930e-1  ));
    c.push_back(SeriesCoefficient(  2,  1, -0.33032641670203e-4  ));
    c.push_back(SeriesCoefficient(  2,  2, -0.18948987516315e-3  ));
    c.push_back(SeriesCoefficient(  2,  4, -0.39392777243355e-2  ));
    c.push_back(SeriesCoefficient(  2,  7, -0.43797295650573e-1  ));
    c.push_back(SeriesCoefficient(  2, 36, -0.26674547914087e-4  ));
    c.push_back(SeriesCoefficient(  3,  0,  0.20481737692309e-7  ));
    c.push_back(SeriesCoefficient(  3,  1,  0.43870667284435e-6  ));
    c.push_back(SeriesCoefficient(  3,  3, -0.32277677238570e-4  ));
    c.push_back(SeriesCoefficient(  3,  6, -0.15033924542148e-2  ));
    c.push_back(SeriesCoefficient(  3, 35, -0.40668253562649e-1  ));
    c.push_back(SeriesCoefficient(  4,  1, -0.78847309559367e-9  ));
    c.push_back(SeriesCoefficient(  4,  2,  0.12790717852285e-7  ));
    c.push_back(SeriesCoefficient(  4,  3,  0.48225372718507e-6  ));
    c.push_back(SeriesCoefficient(  5,  7,  0.22922076337661e-5  ));
    c.push_back(SeriesCoefficient(  6,  3, -0.16714766451061e-10 ));
    c.push_back(SeriesCoefficient(  6, 16, -0.21171472321355e-2  ));
    c.push_back(SeriesCoefficient(  6, 35, -0.23895741934104e2   ));
    c.push_back(SeriesCoefficient(  7,  0, -0.59059564324270e-17 ));
    c.push_back(SeriesCoefficient(  7, 11, -0.12621808899101e-5  ));
    c.push_back(SeriesCoefficient(  7, 25, -0.38946842435739e-1  ));
    c.push_back(SeriesCoefficient(  8,  8,  0.11256211360459e-10 ));
    c.push_back(SeriesCoefficient(  8, 36, -0.82311340897998e1   ));
    c.push_back(SeriesCoefficient(  9, 13,  0.19809712802088e-7  ));
    c.push_back(SeriesCoefficient( 10,  4,  0.10406965210174e-18 ));
    c.push_back(SeriesCoefficient( 10, 10, -0.10234747095929e-12 ));
    c.push_back(SeriesCoefficient( 10, 14, -0.10018179379511e-8  ));
    c.push_back(SeriesCoefficient( 16, 29, -0.80882908646985e-10 ));
    c.push_back(SeriesCoefficient( 16, 50,  0.10693031879409     ));
    c.push_back(SeriesCoefficient( 18, 57, -0.33662250574171     ));
    c.push_back(SeriesCoefficient( 20, 20,  0.89185845355421e-24 ));
    c.push_back(SeriesCoefficient( 20, 35,  0.30629316876232e-12 ));
    c.push_back(SeriesCoefficient( 20, 48, -0.42002467698208e-5  ));
    c.push_back(SeriesCoefficient( 21, 21, -0.59056029685639e-25 ));
    c.push_back(SeriesCoefficient( 22, 53,  0.37826947613457e-5  ));
    c.push_back(SeriesCoefficient( 23, 39, -0.12768608934681e-14 ));
    c.push_back(SeriesCoefficient( 24, 26,  0.73087610595061e-28 ));
    c.push_back(SeriesCoefficient( 24, 40,  0.55414715350778e-16 ));
    c.push_back(SeriesCoefficient( 24, 58, -0.94369707241210e-6  ));
  }

  //
  // destructor
  //
  ~Region2() {}

  //
  // evaluate the region
  //
  void getvars( const double pressure,
                const double temperature,
                std::vector<double> & v );

};

class Region3 : public RegionRT
{
public:
  //
  // constructor
  //
  // densities between 40-800 kg/m^3 bounds the pressure region
  // between 16.5-100 MPa from 600-900 K
  //
  Region3()
    : RegionRT( 322.0, 647.096, 0.461526, 40., 800., 615., 900. ) {
    c.push_back(SeriesCoefficient(  0,  0,  0.10658070028513e1  ));
    c.push_back(SeriesCoefficient(  0,  0, -0.15732845290239e2  ));
    c.push_back(SeriesCoefficient(  0,  1,  0.20944396974307e2  ));
    c.push_back(SeriesCoefficient(  0,  2, -0.76867707878716e1  ));
    c.push_back(SeriesCoefficient(  0,  7,  0.26185947787954e1  ));
    c.push_back(SeriesCoefficient(  0, 10, -0.28080781148620e1  ));
    c.push_back(SeriesCoefficient(  0, 12,  0.12053369696517e1  ));
    c.push_back(SeriesCoefficient(  0, 23, -0.84566812812502e-2 ));
    c.push_back(SeriesCoefficient(  1,  2, -0.12654315477714e1  ));
    c.push_back(SeriesCoefficient(  1,  6, -0.11524407806681e1  ));
    c.push_back(SeriesCoefficient(  1, 15,  0.88521043984318    ));
    c.push_back(SeriesCoefficient(  1, 17, -0.64207765181607    ));
    c.push_back(SeriesCoefficient(  2,  0,  0.38493460186671    ));
    c.push_back(SeriesCoefficient(  2,  2, -0.85214708824206    ));
    c.push_back(SeriesCoefficient(  2,  6,  0.48972281541877e1  ));
    c.push_back(SeriesCoefficient(  2,  7, -0.30502617256965e1  ));
    c.push_back(SeriesCoefficient(  2, 22,  0.39420536879154e-1 ));
    c.push_back(SeriesCoefficient(  2, 26,  0.12558408424308    ));
    c.push_back(SeriesCoefficient(  3,  0, -0.27999329698710    ));
    c.push_back(SeriesCoefficient(  3,  2,  0.13899799569460e1  ));
    c.push_back(SeriesCoefficient(  3,  4, -0.20189915023570e1  ));
    c.push_back(SeriesCoefficient(  3, 16, -0.82147637173963e-2 ));
    c.push_back(SeriesCoefficient(  3, 26, -0.47596035734923    ));
    c.push_back(SeriesCoefficient(  4,  0,  0.43984074473500e-1 ));
    c.push_back(SeriesCoefficient(  4,  2, -0.44476435428739    ));
    c.push_back(SeriesCoefficient(  4,  4,  0.90572070719733    ));
    c.push_back(SeriesCoefficient(  4, 26,  0.70522450089767    ));
    c.push_back(SeriesCoefficient(  5,  1,  0.10770512626332    ));
    c.push_back(SeriesCoefficient(  5,  3, -0.32913623258954    ));
    c.push_back(SeriesCoefficient(  5, 26, -0.50871062041158    ));
    c.push_back(SeriesCoefficient(  6,  0, -0.22175400873096e-1 ));
    c.push_back(SeriesCoefficient(  6,  2,  0.94260751665092e-1 ));
    c.push_back(SeriesCoefficient(  6, 26,  0.16436278447961    ));
    c.push_back(SeriesCoefficient(  7,  2, -0.13503372241348e-1 ));
    c.push_back(SeriesCoefficient(  8, 26, -0.14834345352472e-1 ));
    c.push_back(SeriesCoefficient(  9,  2,  0.57922953628084e-3 ));
    c.push_back(SeriesCoefficient(  9, 26,  0.32308904703711e-2 ));
    c.push_back(SeriesCoefficient( 10,  0,  0.80964802996215e-4 ));
    c.push_back(SeriesCoefficient( 10,  1, -0.16557679795037e-3 ));
    c.push_back(SeriesCoefficient( 11, 26, -0.44923899061815e-4 ));
    pcrit = 22064.;
  }

  //
  // destructor
  //
  ~Region3() {}

  //
  // evaluate the region
  //
  virtual void getvars( const double density,
                        const double temperature,
                        std::vector<double> & v );

  virtual void getvars2( const double density,
                         const double temperature,
                         const Vars xtype,
                         std::vector<double> & v );

  double getSpinodalDensityBound( const double temperature,
                                  const bool glstate );

};

class Region5 : public RegionPT
{
public:
  //
  // constructor
  //
  Region5()
    : RegionPT( 1000.0, 1000.0, 0.461526 ) {
    c.push_back(SeriesCoefficient( 0,  0, -0.13179983674201e2  ));
    c.push_back(SeriesCoefficient( 0,  1,  0.68540841634434e1  ));
    c.push_back(SeriesCoefficient( 0, -3, -0.24805148933466e-1 ));
    c.push_back(SeriesCoefficient( 0, -2,  0.36901534980333    ));
    c.push_back(SeriesCoefficient( 0, -1, -0.31161318213925e1  ));
    c.push_back(SeriesCoefficient( 0,  2, -0.32961626538917    ));
    c.push_back(SeriesCoefficient( 1,  1,  0.15736404855259e-2 ));
    c.push_back(SeriesCoefficient( 1,  2,  0.90153761673944e-3 ));
    c.push_back(SeriesCoefficient( 1,  3, -0.50270077677648e-2 ));
    c.push_back(SeriesCoefficient( 2,  3,  0.22440037409485e-5 ));
    c.push_back(SeriesCoefficient( 2,  9, -0.41163275453471e-5 ));
    c.push_back(SeriesCoefficient( 3,  7,  0.37919454822955e-7 ));
  }

  //
  // destructor
  //
  ~Region5() {}

  //
  // evaluate the region
  //
  virtual void getvars( const double pressure,
                        const double temperature,
                        std::vector<double> & v );

};

class Region23Boundary
{
public:

  //
  // constructor
  //
  Region23Boundary()
    : c(5), ps(1000.0), ts(1.0) {
    c[0] =  0.34805185628969e3;
    c[1] = -0.11671859879975e1;
    c[2] =  0.10192970039326e-2;
    c[3] =  0.57254459862746e3;
    c[4] =  0.13918839778870e2;
  }

  //
  // destructor
  //
  ~Region23Boundary() {}

  //
  // P or T evaluation
  //
  double getT(const double pressure) {
    double p = pressure/ps;
    return ts*(c[3]+sqrt((p-c[4])/c[2]));
  }

  double getP(const double temperature) {
    double t = temperature/ts;
    return ps*(c[0]+t*(c[1]+t*c[2]));
  }

private:

  std::vector<double> c;
  double ps,ts;

};

class Region4
{
public:

  //
  // constructor
  //
  Region4()
    : c(10), ps(1000.0), ts(1.0) {
    c[0] =  0.11670521452767e4;
    c[1] = -0.72421316703206e6;
    c[2] = -0.17073846940092e2;
    c[3] =  0.12020824702470e5;
    c[4] = -0.32325550322333e7;
    c[5] =  0.14915108613530e2;
    c[6] = -0.48232657361591e4;
    c[7] =  0.40511340542057e6;
    c[8] = -0.23855557567849  ;
    c[9] =  0.65017534844798e3;
  }

  //
  // destructor
  //
  ~Region4() {}

  //
  // P or T evaluation
  //
  double getP(const double temperature) {
    double t = temperature/ts+c[8]/(temperature/ts-c[9]);
    double A = c[1]+t*(c[0]+t);
    double B = c[4]+t*(c[3]+c[2]*t);
    double C = c[7]+t*(c[6]+c[5]*t);

    double p = 2*C/(-B+sqrt(B*B-4*A*C));
    return ps*p*p*p*p;
  }

  double getT(const double pressure) {
    double b = sqrt(sqrt(pressure/ps));
    double E = c[5]+b*(c[2]+b);
    double F = c[6]+b*(c[3]+b*c[0]);
    double G = c[7]+b*(c[4]+b*c[1]);
    double D = -2*G/(F+sqrt(F*F-4.*E*G));

    return ts*(c[9]+D-sqrt((c[9]+D)*(c[9]+D)-4.*(c[8]+c[9]*D)))/2.;
  }

private:

  std::vector<double> c;
  double ps,ts;

};

class Region3Spinodal: public DFMin
{
public:
  //
  // constructor
  //
  Region3Spinodal( const double temperature,
                   const int direction,
                   Region3 & region3 )
    : r3(region3), t(temperature), dir(direction)
  {
    if (dir != 1 && dir != -1) throw std::runtime_error("Region3Spinodal: bad direction passed to constructor");
  }

  //
  // destructor
  //
  ~Region3Spinodal() {}

  void ResetT( const double temperature ) {
    t = temperature;
  }

  void ResetDir( const int direction ) {
    if (direction != 1 && direction != -1) throw std::runtime_error("Region3Spinodal::ResetDir: bad direction input");
    dir = direction;
  }

protected:
  
  double f( const double rho ) {
    std::vector<double> v;
    r3.getvars(rho,t,v);
    return dir*v[V_P];
  }

private:

  Region3 & r3;
  double t;
  int dir;

};

class IAPWSIF97 : public EOSModel
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
  IAPWSIF97( const Parameters & params,
             const std::string name );
  IAPWSIF97( );

  //
  // destructor
  //
  ~IAPWSIF97();

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
  // get the saturation temperature for regions 1 and 2 at a given pressure
  //
  double getR12satT( const double pressure );

  //
  // get the saturation temperature for region 3 at a given pressure
  //
  double getR3satT( const double pressure );

  //
  // get the state on the Maxwell construction between regions 1 and 2
  // for the given pressure and xvar. temperature must be from a call
  // to getR12satT to ensure one is actually on the construction.
  //
  void getR12maxwellPX( const double pressure,
                        const double temperature,
                        const double xvar,
                        const Vars xtype,
                        std::vector<double> & v );

  //
  // get the state on the Maxwell construction for region 3
  // for the given pressure and xvar. temperature must be from a call
  // to getR3satT to ensure one is actually on the construction.
  //
  void getR3maxwellPX( const double pressure,
                       const double temperature,
                       const double xvar,
                       const Vars xtype,
                       std::vector<double> & v );

  //
  // get the saturation densities on the Maxwell construction for Region 3
  //
  void getR3maxwellR( const double temperature,
                      double & rho1,
                      double & rho2 );

  //
  // get derivative information for vars along the l-v phase line
  //
  void getR12maxwellDerivs( const double pressure,
                            std::vector<double> & vd1,
                            std::vector<double> & vd2 );

  void getR3maxwellDerivs( const double temperature,
                           std::vector<double> & vd1,
                           std::vector<double> & vd2 );

  //
  // compute the viscosity and thermal conductivity for the given
  // state. This assumes a valid thermodynamic state has already been
  // calculated in v from a prior call to either getvarsPT or
  // getvarsPX
  //
  void computeTransportProperties( std::vector<double> & v );

  double computeRefStateK( const double r );

private:

  //
  // convert units of the getvars calls to SI
  //
  void convert_to_SI( std::vector<double> & v );

  //
  // prevent extreme unphysical states in derivatives
  //
  void clip_derivatives( std::vector<double> & v);

  void setRefStateCoef();

  Region1 r1;
  Region2 r2;
  Region3 r3;
  Region4 r4;
  Region5 r5;
  Region23Boundary r23b;
  IAPWSvisc2008 visc;
  IAPWSthcon2011 thcon;

  //
  // allow extrapolation in temperature
  //
  bool Textrap;

  //
  // coefficients for transport property reference bulk modulus computation
  //
  std::vector<std::vector<double> > krc;
  double krcs;
};

#endif
