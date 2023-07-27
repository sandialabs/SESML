/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef IAPWS_TRANSPORT_H
#define IAPWS_TRANSPORT_H

#include <cmath>
#include <stdexcept>

#include "Vars.H"
#include "Region.H"

class IAPWSvisc2008
{
public:

  //
  // constructor
  //
  IAPWSvisc2008()
    : ts(647.096), rs(322.), ps(22064.), us(1.e-6), c0(0), c1(0),
      xmu(0.068), qc(1./1.9), qd(1./1.1), xie(0.63/1.239), xi0(0.13),g0(0.06),tr(1.5) {
    // dilute-gas part
    c0.push_back(SeriesCoefficient( 0, 0,  1.67752   ));
    c0.push_back(SeriesCoefficient( 1, 0,  2.20462   ));
    c0.push_back(SeriesCoefficient( 2, 0,  0.6366564 ));
    c0.push_back(SeriesCoefficient( 3, 0, -0.241605  ));

    // finite density part
    c1.push_back(SeriesCoefficient( 0, 0,  5.20094e-1 ));
    c1.push_back(SeriesCoefficient( 1, 0,  8.50895e-2 ));
    c1.push_back(SeriesCoefficient( 2, 0, -1.08374    ));
    c1.push_back(SeriesCoefficient( 3, 0, -2.89555e-1 ));
    c1.push_back(SeriesCoefficient( 0, 1,  2.22531e-1 ));
    c1.push_back(SeriesCoefficient( 1, 1,  9.99115e-1 ));
    c1.push_back(SeriesCoefficient( 2, 1,  1.88797    ));
    c1.push_back(SeriesCoefficient( 3, 1,  1.26613    ));
    c1.push_back(SeriesCoefficient( 5, 1,  1.20573e-1 ));
    c1.push_back(SeriesCoefficient( 0, 2, -2.81378e-1 ));
    c1.push_back(SeriesCoefficient( 1, 2, -9.06851e-1 ));
    c1.push_back(SeriesCoefficient( 2, 2, -7.72479e-1 ));
    c1.push_back(SeriesCoefficient( 3, 2, -4.89837e-1 ));
    c1.push_back(SeriesCoefficient( 4, 2, -2.57040e-1 ));
    c1.push_back(SeriesCoefficient( 0, 3,  1.61913e-1 ));
    c1.push_back(SeriesCoefficient( 1, 3,  2.57399e-1 ));
    c1.push_back(SeriesCoefficient( 0, 4, -3.25372e-2 ));
    c1.push_back(SeriesCoefficient( 3, 4,  6.98452e-2 ));
    c1.push_back(SeriesCoefficient( 4, 5,  8.72102e-3 ));
    c1.push_back(SeriesCoefficient( 3, 6, -4.35673e-3 ));
    c1.push_back(SeriesCoefficient( 5, 6, -5.93264e-4 ));

  }

  //
  // destructor
  //
  ~IAPWSvisc2008() {}

  //
  // evaluate the Helmholtz free energy
  //
  double getViscosity( const double density,
                       const double temperature,
                       const double kt,
                       const double kt15 );

  //
  // evaluate the u0 and u1 terms
  //
  double getViscosityU0U1( const double density,
                           const double temperature );

private:

  //
  // reference quantities
  //
  double ts,rs,ps,us;

  //
  // coefficients for u0 and u1 parts
  //
  std::vector<SeriesCoefficient> c0;
  std::vector<SeriesCoefficient> c1;

  //
  // coefficients for u2 part
  //
  double xmu,qc,qd,xie,xi0,g0,tr;

};


class IAPWSthcon2011
{
public:

  //
  // constructor
  //
  IAPWSthcon2011()
    : ts(647.096), rs(322.), ps(22064.), ls(1.e-3), us(1.e-6), R(0.46151805), c0(0), c1(0),
      L(177.8514), qd(1./0.4), xie(0.63/1.239), xi0(0.13), g0(0.06), tr(1.5) {
    // dilute-gas part
    c0.push_back(SeriesCoefficient( 0, 0,  2.443221e-3 ));
    c0.push_back(SeriesCoefficient( 1, 0,  1.323095e-2 ));
    c0.push_back(SeriesCoefficient( 2, 0,  6.770357e-3 ));
    c0.push_back(SeriesCoefficient( 3, 0, -3.454586e-3 ));
    c0.push_back(SeriesCoefficient( 4, 0,  4.096266e-4 ));

    // finite density part
    c1.push_back(SeriesCoefficient( 0, 0,  1.60397357 ));
    c1.push_back(SeriesCoefficient( 1, 0,  2.33771842 ));
    c1.push_back(SeriesCoefficient( 2, 0,  2.19650529 ));
    c1.push_back(SeriesCoefficient( 3, 0, -1.21051378 ));
    c1.push_back(SeriesCoefficient( 4, 0, -2.7203370  ));
    c1.push_back(SeriesCoefficient( 0, 1, -0.646013523 ));
    c1.push_back(SeriesCoefficient( 1, 1, -2.78843778  ));
    c1.push_back(SeriesCoefficient( 2, 1, -4.54580785  ));
    c1.push_back(SeriesCoefficient( 3, 1,  1.60812989  ));
    c1.push_back(SeriesCoefficient( 4, 1,  4.57586331  ));
    c1.push_back(SeriesCoefficient( 0, 2,  0.111443906 ));
    c1.push_back(SeriesCoefficient( 1, 2,  1.53616167  ));
    c1.push_back(SeriesCoefficient( 2, 2,  3.55777244  ));
    c1.push_back(SeriesCoefficient( 3, 2, -0.621178141 ));
    c1.push_back(SeriesCoefficient( 4, 2, -3.18369245  ));
    c1.push_back(SeriesCoefficient( 0, 3,  0.102997357  ));
    c1.push_back(SeriesCoefficient( 1, 3, -0.463045512  ));
    c1.push_back(SeriesCoefficient( 2, 3, -1.40944978   ));
    c1.push_back(SeriesCoefficient( 3, 3,  0.0716373224 ));
    c1.push_back(SeriesCoefficient( 4, 3,  1.1168348    ));
    c1.push_back(SeriesCoefficient( 0, 4, -0.0504123634 ));
    c1.push_back(SeriesCoefficient( 1, 4,  0.0832827019 ));
    c1.push_back(SeriesCoefficient( 2, 4,  0.275418278  ));
    c1.push_back(SeriesCoefficient( 4, 4, -0.19268305   ));
    c1.push_back(SeriesCoefficient( 0, 5,  0.00609859258 ));
    c1.push_back(SeriesCoefficient( 1, 5, -0.00719201245 ));
    c1.push_back(SeriesCoefficient( 2, 5, -0.0205938816  ));
    c1.push_back(SeriesCoefficient( 4, 5,  0.012913842   ));

  }

  //
  // destructor
  //
  ~IAPWSthcon2011() {}

  //
  // evaluate the Helmholtz free energy
  //
  double getThermalConductivity( const double density,
                                 const double temperature,
                                 const double kt,
                                 const double kt15,
                                 const double cv,
                                 const double cp,
                                 const double mu );

  //
  // evaluate the u0 and u1 terms
  //
  double getThermalConductivityL0L1( const double density,
                                     const double temperature );

private:

  //
  // reference quantities
  //
  double ts,rs,ps,ls,us,R;

  //
  // coefficients for l0 and l1 parts
  //
  std::vector<SeriesCoefficient> c0;
  std::vector<SeriesCoefficient> c1;

  //
  // coefficients for l2 part
  //
  double L,qd,xie,xi0,g0,tr;

};

#endif
