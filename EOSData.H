/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef EOSDATA_H
#define EOSDATA_H

#include <vector>

class EOSData
{
public:

  //
  // allowed input variables
  //
  enum InputVars {
    IV_R,  // density
    IV_T,  // temperature
    IV_P,  // pressure
    IV_E,  // internal energy
    IV_S,  // entropy
    IV_TOT
  };

  //
  // allowed output variables
  //
  enum OutputVars {
    OV_R,  // density
    OV_T,  // temperature
    OV_P,  // pressure
    OV_E,  // internal energy
    OV_S,  // entropy
    OV_F,  // Helmholtz free energy
    OV_G,  // Gibb's free energy
    OV_H,  // enthalpy
    OV_CS, // adiabatic sound speed
    OV_B,  // isothermal bulk modulus
    OV_B2, // second isothermal bulk modulus (R DBDR_T)
    // following derivatives use above abbrieviations
    // variables after last underscore denote quantities held constant
    OV_DCSDR_S,
    OV_DPDT_R,
    OV_DPDR_T,
    OV_DEDT_T,
    OV_DEDR_T,
    OV_D2PDT2_RR,
    OV_D2EDT2_RR,
    OV_D2EDTDR_RT,
    OV_TOT
  };

  //
  // Constructor
  //
  EOSData( const unsigned int ninput,
           const unsigned int noutput );

  //
  // Destructor
  //
  ~EOSData();

  //
  // data holders
  //
  std::vector<double> inputs;
  std::vector<double> outputs;
  
};

#endif
