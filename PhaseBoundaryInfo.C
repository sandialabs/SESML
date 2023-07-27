/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <list>
#include <algorithm>

#include "PhaseBoundaryInfo.H"
#include "FiniteDiffDeriv.H"

//
// Constructor for case where splines are not used
// lb and ub contain R1,R2,T
//
PhaseBoundaryInfo::PhaseBoundaryInfo( EOSModel * model,
                                      std::vector<double> & lb,
                                      std::vector<double> & ub,
                                      const std::vector<int> & phases )
  : ph(phases), m(model)
{
  if (m == NULL) throw std::runtime_error("NULL argument passed to PhaseBoundaryInfo");
  if (lb.size() < 3 || ub.size() < 3) {
    throw std::runtime_error("PhaseBoundaryInfo: lb or ub input vectors too small");
  }
  if (lb[2] >= ub[2]) {
    throw std::runtime_error("PhaseBoundaryInfo: upper bound temperature must be greater than lower bound temperature");
  }
  if (phases.size() != 5) throw std::runtime_error("PhaseBoundaryInfo: phases size must be 5");

  //
  // set bounds
  //
  minT = lb[2];
  minR[0] = lb[0];
  minR[1] = lb[1];
  maxT = ub[2];
  maxR[0] = ub[0];
  maxR[1] = ub[1];

}

//
// Method to get the phase boundary location when not stored as spline
//
void PhaseBoundaryInfo::getExactVars( const double T,
                                      double & r1,
                                      double & r2 )
{
  if (m == NULL) throw std::runtime_error("PhaseBoundaryInfo::getExactVars: model is NULL");

  double pressure;
  double r[2];
  m->getPhaseCoexistence(ph[1],T,pressure,r);
  r1 = r[0];
  r2 = r[1];

}

PhaseBoundaryInfo::~PhaseBoundaryInfo()
{

}

