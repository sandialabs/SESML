/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef EOSMODEL_H
#define EOSMODEL_H

#include <vector>
#include <string>
#include <limits>

#include "EOSData.H"
#include "EOSParam.H"

// Forward declare the PhaseBoundaryInfo class
class PhaseBoundaryInfo;

class EOSModel
{
public:
  //
  // Constructor
  //
  EOSModel( std::string name = "EOSModel stub" );

  //
  // Destructor
  //
  virtual ~EOSModel();

  //
  // Method to obtain information from cold curve model
  //
  virtual void getColdCurve( EOSData & data );

  //
  // Methods to obtain information from a thermal model. Character
  // following the underscore denotes the thermodynamic potential to
  // draw from. This implies the requirements for the first two inputs
  // and first four outputs as shown in the following table, assuming
  // getThermalEOS_X with inputs I[1-2] and outputs O[1-5] in data.
  //
  // X | I1 I2 | O1 O2 O3 O4  O5
  // --+-------+------------------
  // F |  R  T |  P  E  S  F phase
  // S |  R  E |  P  T  S  F phase
  // G |  P  T |  R  E  S  F phase
  //
  virtual void getThermalEOS_F( EOSData & data );
  virtual void getThermalEOS_S( EOSData & data );
  virtual void getThermalEOS_G( EOSData & data );

  //
  // Overloaded methods as above, except outputs are computed for the specified phase
  //
  virtual void getThermalEOS_F( EOSData & data, const int /*phase*/ ) { getThermalEOS_F(data); }
  virtual void getThermalEOS_S( EOSData & data, const int /*phase*/ ) { getThermalEOS_S(data); }
  virtual void getThermalEOS_G( EOSData & data, const int /*phase*/ ) { getThermalEOS_G(data); }

  //
  // Methods to obtain information from a thermal model. These methods
  // return higher order derivative information in addtion to the
  // first derivatives fo the potential. The first inputs and outputs
  // are as for the above methods with the additional output
  // derivatives as shown. For O6-O9, either R or T is held constant,
  // as appropriate. O10 is the adiabatic sound speed.
  //
  // O6 O7 O8 O9 O10 O11
  // -------------------
  // dP dP dE dE  CS dCS
  // dR dT dR dT     dR
  virtual void getThermalEOS_F_D( EOSData & data );
  virtual void getThermalEOS_S_D( EOSData & data );
  virtual void getThermalEOS_G_D( EOSData & data );

  //
  // Overloaded methods as above, except outputs are computed for the specified phase
  //
  virtual void getThermalEOS_F_D( EOSData & data, const int /*phase*/ ) { getThermalEOS_F_D(data); }
  virtual void getThermalEOS_S_D( EOSData & data, const int /*phase*/ ) { getThermalEOS_S_D(data); }
  virtual void getThermalEOS_G_D( EOSData & data, const int /*phase*/ ) { getThermalEOS_G_D(data); }

  //
  // Methods to obtain information from a thermal model. These models
  // are different from the above in that they compute a different set
  // of functions to support pressure based solvers and the
  // independent variables may not be defined by a usual thermodynamic
  // potential. The two characters after the first underscore denote
  // the independent variables, also listed in I[1-2] and the
  // underscore D denotes the computation of derivatives. Outputs are
  // as listed in the O[1-11] columns.
  //
  // X  | I1 I2 | O1 O2 O3 O4 O5 O6 O7 O8 O9 O10  O11
  // ---+-------+-------------------------------------
  // PT |  P  T |  R  H  E  S  G  F CP CV CS  KT phase
  // PH |  P  H |  R  T  E  S  G  F CP CV CS  KT phase
  // PS |  P  S |  R  T  E  H  G  F CP CV CS  KT phase 
  // PE |  P  E |  R  T  H  S  G  F CP CV CS  KT phase
  // RT |  R  T |  P  H  E  S  G  F CP CV CS  KT phase
  //
  virtual void getThermalEOS_PT_D( EOSData & data, const int phase );
  virtual void getThermalEOS_PH_D( EOSData & data, const int phase );
  virtual void getThermalEOS_PS_D( EOSData & data, const int phase );
  virtual void getThermalEOS_PE_D( EOSData & data, const int phase );
  virtual void getThermalEOS_RT_D( EOSData & data, const int phase );

  //
  // Method to request calculation of the phase boundary information for a model
  //
  // Temperature bounds give the range of temperatures within which
  // boundary calculations are desired.
  //
  virtual void calculatePhaseBoundaryInfo( const double /*lowT*/,
                                           const double /*highT*/ ) {
    // default behavior is to not calculate any phase boundaries
  }

  //
  // Method to obtain the phase boundary information from a model
  //
  std::vector<PhaseBoundaryInfo *> * getPhaseBoundaryInfo() {
    return &phase_boundaries;
  }

  //
  // Method to obtain the main variables at the critical point of a model, if supported
  //
  virtual void getCriticalVars( double & temperature,
                                double & pressure,
                                double & density );

  //
  // Method to obtain the main variables at the triple point of a model, if supported
  //
  virtual void getTriplePtVars( const int phase,
                                double & temperature,
                                double & pressure,
                                double * density );

  //
  // Method to set the temperature bounds
  //
  void setTemperatureBounds( const double tmin = 0.0,
                             const double tmax = std::numeric_limits<double>::max() )
  {
    minimum_temperature = tmin;
    maximum_temperature = tmax;
  }

  //
  // Method to obtain particular phase boundary information
  //
  virtual void getPhaseCoexistence( const int phase,
                                    const double temperature,
                                    double & pressure,
                                    double * density );

protected:

  //
  // Method to retrieve parameters for a model
  //
  const EOSParam * getParams( const Parameters & params,
                              const std::string name );

  const std::string model_name;

  //
  // maximium and minimum allowed densities for the model
  //
  double maximum_density;
  double minimum_density;

  //
  // maximum and minimum allowed temperature for the model
  // used for faster inversion routines, not all thermal methods may respect these
  //
  double maximum_temperature;
  double minimum_temperature;

  //
  // phase boundaries for model
  //
  std::vector<PhaseBoundaryInfo *> phase_boundaries;

};

#endif
