/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef REGION_H
#define REGION_H

#include <vector>

#include "Vars.H"

struct SeriesCoefficient
{
  //
  // convenience constructor
  //
  SeriesCoefficient( const int Ii,
                     const int Ji,
                     const double ni )
    : I(Ii), J(Ji), n(ni) {}

  //
  // null constructor
  //
  SeriesCoefficient()
    : I(0), J(0), n(0.0) {}

  //
  // two series exponents
  //
  int I;
  int J;

  //
  // series coefficient value
  //
  double n;

};

class RegionPT
{
public:

  //
  // constructor
  //
  RegionPT( const double pressure_reference,
            const double temperature_reference,
            const double gas_constant )
    : c(0), ps(pressure_reference), ts(temperature_reference), R(gas_constant) {}

  //
  // destructor
  //
  virtual ~RegionPT() {}

  //
  // access method
  //
  virtual void getvars( const double pressure,
                        const double temperature,
                        std::vector<double> & v ) = 0;


  //
  // method to find temperature which gives a desired xtype within the
  // temperature bounds
  //
  double getTatPX( const double pressure,
                   const double xvar,
                   const int xtype,
                   const double minT,
                   const double maxT );

protected:

  //
  // coefficients
  //
  std::vector<SeriesCoefficient> c;

  //
  // reference pressure and temperature
  //
  double ps,ts;

  //
  // gas constant
  //
  double R;

};

class RegionRT
{
public:

  //
  // constructor
  //
  RegionRT( const double density_reference,
            const double temperature_reference,
            const double gas_constant,
            const double min_density,
            const double max_density,
            const double min_temperature,
            const double max_temperature )
    : c(0), rs(density_reference), ts(temperature_reference), pcrit(0.),
      R(gas_constant), rbounds(2,0.), tbounds(2,0.)
  {
    rbounds[0] = min_density;
    rbounds[1] = max_density;
    if (rbounds[0] >= rbounds[1] || rbounds[0] < 0.) throw std::runtime_error("RegionRT: invalid density bounds");

    tbounds[0] = min_temperature;
    tbounds[1] = max_temperature;
    if (tbounds[0] >= tbounds[1] || tbounds[0] < 0.) throw std::runtime_error("RegionRT: invalid temperature bounds");
  }

  //
  // destructor
  //
  virtual ~RegionRT() {}

  //
  // access method for main variables
  //
  virtual void getvars( const double density,
                        const double temperature,
                        std::vector<double> & v ) = 0;

  //
  // access method for variables and derivatives used in getRTatPX solver
  //
  // on return, v must contain:
  // 0 1  2    3    4    5
  // P X dPdR dPdT dXdR dXdT
  // for X in E,H,S,G
  //
  virtual void getvars2( const double density,
                         const double temperature,
                         const Vars xtype,
                         std::vector<double> & v ) = 0;

  //
  // method to find density and temperature which gives a desired
  // xtype and pressure within the temperature bounds
  //
  void getRTatPX( const double pressure,
                  const double xvar,
                  const Vars xtype,
                  const double density0,
                  const double temperature0,
                  double & density,
                  double & temperature );

  //
  // method to find density at given pressure and temperature
  //
  double getRatPT( const double pressure,
                   const double temperature,
                   const double minR,
                   const double maxR,
                   const int side = 0 );

  //
  // at the given temperature and the following values of glstate:
  // [true] - return an upper bound on the gas spinodal density
  // [false] - return a lower bound on the liquid spinodal density
  //
  virtual double getSpinodalDensityBound( const double /*temperature*/,
                                          const bool /*glstate*/ )
  {
    return rs;
  }

protected:

  //
  // coefficients
  //
  std::vector<SeriesCoefficient> c;

  //
  // reference density and temperature
  //
  double rs,ts;

  //
  // critical pressure
  //
  double pcrit;

  //
  // gas constant
  //
  double R;

  //
  // min/max density and temperature bounds for the getR* routines
  //
  std::vector<double> rbounds;
  std::vector<double> tbounds;

};

#endif
