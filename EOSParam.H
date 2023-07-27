/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef EOSPARAM_H
#define EOSPARAM_H

#include <vector>
#include <string>
#include <map>
#include <ostream>

class EOSParam
{
public:
  //
  // Constructor
  //
  EOSParam( const unsigned int ndouble,
            const unsigned int nint,
            const unsigned int nstring );

  //
  // Destructor
  //
  ~EOSParam();

  //
  // model the parameters describe
  //
  std::string model_name;

  //
  // param holders
  //
  std::vector<double>      params_double;
  std::vector<int>         params_int;
  std::vector<std::string> params_string;

};

class Parameters : public std::map<std::string,EOSParam *>
{
public:
  Parameters() {}
  Parameters(const Parameters & p);
  ~Parameters();
};

std::ostream & operator << ( std::ostream & out, Parameters & p );

#endif
