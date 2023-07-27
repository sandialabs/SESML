/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef PHASE_BOUNDARY_INFO_H
#define PHASE_BOUNDARY_INFO_H

#include <vector>
#include <cmath>

#include "EOSData.H"
#include "EOSModel.H"

class PhaseBoundaryInfo
{
public:

  enum ThermoVars {
    P,
    R1,
    R2,
    E1,
    E2,
    S1,
    S2,
    F1,
    F2,
    EP1,
    EP2,
    EPP1,
    EPP2,
    TOT_THERMOVARS
  };

  //
  // Constructor for no spline case
  //
  PhaseBoundaryInfo( EOSModel * m,
                     std::vector<double> & lb,
                     std::vector<double> & ub,
                     const std::vector<int> & phases );

  //
  // Destructor
  //
  ~PhaseBoundaryInfo();

  void getExactVars( const double T,
                     double & r1,
                     double & r2 );

  void getBounds( double & T1,
                  double & T2,
                  double * R1,
                  double * R2 )
  {
    T1 = minT;
    T2 = maxT;
    R1[0] = minR[0];
    R1[1] = minR[1];
    R2[0] = maxR[0];
    R2[1] = maxR[1];
  }

  void getPhases( int & phase1,
                  int & phasemix,
                  int & phase2,
                  int & phaselow,
                  int & phasehigh )
  {
    phase1 = ph[0];
    phasemix = ph[1];
    phase2 = ph[2];
    phaselow = ph[3];
    phasehigh = ph[4];
  }

  void getCriticalVars( double & t,
                        double & b,
                        double & d,
                        std::vector<double> & ca,
                        std::vector<double> & cc )
  {
    t = Tc;
    // beta < 0 indicates the exponent was calculated but not used
    b = fabs(beta);
    d = delta;
    ca = a;
    cc = c;
  }

private:

  //
  // exponents for critical point scaling
  //
  double beta;
  double delta;

  //
  // critical temperature
  //
  double Tc;

  //
  // thermodynamic variable scaling amplitudes and critical point values
  //
  std::vector<double> a;
  std::vector<double> c;

  //
  // min and max temperature of phase boundary
  //
  double minT, maxT;
  double minR[2],maxR[2];

  //
  // list of phases that boundary represents
  //
  std::vector<int> ph;

  //
  // EOS model for use when not approximating with splines (NULL otherwise)
  //
  EOSModel * m;

};

#endif
