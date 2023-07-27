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
#include <sstream>
#include <limits>
#include <cmath>

#include "EOSModel.H"

//
// Constructor
//
EOSModel::EOSModel( std::string name ) :
  model_name(name), maximum_density(std::numeric_limits<double>::max()),
  minimum_density(0.0), maximum_temperature(std::numeric_limits<double>::max()),
  minimum_temperature(0.0)
{

}

//
// Destructor
//
EOSModel::~EOSModel()
{

}

//
// Stub cold curve routine
//
void EOSModel::getColdCurve( EOSData & )
{
  throw std::logic_error("getColdCurve not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal EOS routine
//
void EOSModel::getThermalEOS_F( EOSData & )
{
  throw std::logic_error("getThermalEOS_F not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal EOS routine
//
void EOSModel::getThermalEOS_S( EOSData & )
{
  throw std::logic_error("getThermalEOS_S not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal EOS routine
//
void EOSModel::getThermalEOS_G( EOSData & )
{
  throw std::logic_error("getThermalEOS_G not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_F_D( EOSData & )
{
  throw std::logic_error("getThermalEOS_F_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_S_D( EOSData & )
{
  throw std::logic_error("getThermalEOS_S_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_G_D( EOSData & )
{
  throw std::logic_error("getThermalEOS_G_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_PT_D( EOSData & , const int )
{
  throw std::logic_error("getThermalEOS_PT_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_PH_D( EOSData & , const int )
{
  throw std::logic_error("getThermalEOS_PH_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_PS_D( EOSData & , const int )
{
  throw std::logic_error("getThermalEOS_PS_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_PE_D( EOSData & , const int )
{
  throw std::logic_error("getThermalEOS_PE_D not implemented for model \""
                         +model_name+'"');
}

//
// Stub thermal derivative EOS routine
//
void EOSModel::getThermalEOS_RT_D( EOSData & , const int )
{
  throw std::logic_error("getThermalEOS_RT_D not implemented for model \""
                         +model_name+'"');
}

//
// Return a pointer to the parameters entry for the given name
//
const EOSParam * EOSModel::getParams( const Parameters & params,
                                      const std::string name )
{
  Parameters::const_iterator p = params.find(name);
  if (p == params.end()) {
    throw std::runtime_error("EOSParam for \""+name+"\" not found");
  }
  return p->second;
}

//
// Stub critical variables routine
//
void EOSModel::getCriticalVars( double &, double &, double & )
{
  throw std::logic_error("getCriticalVars not implemented for model "+model_name);
}

//
// Stub triple point variables routine
//
void EOSModel::getTriplePtVars( const int, double &, double &, double * )
{
  throw std::logic_error("getTriplePtVars not implemented for model "+model_name);
}

//
// Stub phase coexistence routine
//
void EOSModel::getPhaseCoexistence( const int, const double, double &, double * )
{
  throw std::logic_error("getPhaseCoexistence not implemented for model "+model_name);
}

