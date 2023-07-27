/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#include <iomanip>

#include "EOSParam.H"

//
// Constructor
//
EOSParam::EOSParam( const unsigned int ndouble,
                    const unsigned int nint,
                    const unsigned int nstring )
{
  params_double.resize(ndouble,0.0);
  params_int.resize(nint,0);
  params_string.resize(nstring,"");
}

//
// Destructor
//
EOSParam::~EOSParam()
{

}

Parameters::Parameters( const Parameters & p )
  : std::map<std::string,EOSParam *>(p)
{
  for (Parameters::iterator i=begin();i!=end();i++) {
    i->second = new EOSParam( 0, 0, 0 );
    *(i->second) = *(p.find(i->first)->second);
  }
}

Parameters::~Parameters()
{
  //
  // free the model parameters
  //
  for (Parameters::iterator i=begin();i!=end();i++) delete i->second;

}

std::ostream & operator << (std::ostream & out, Parameters & p)
{
  std::streamsize prec = out.precision();
  out << "Parameters content:" << std::endl;
  for (Parameters::iterator i=p.begin();i!=p.end();++i) {
    EOSParam * pi = i->second;
    out << "  Model " << pi->model_name << std::endl;
    out << "    doubles:" << std::endl;
    for (std::size_t j=0;j<pi->params_double.size();++j) {
      out << "      " << j << " " << std::setprecision(20) << pi->params_double[j] << std::endl;
    }
    out << "    ints:" << std::endl;
    for (std::size_t j=0;j<pi->params_int.size();++j) {
      out << "      " << j << " " << pi->params_int[j] << std::endl;
    }
    out << "    strings:" << std::endl;
    for (std::size_t j=0;j<pi->params_string.size();++j) {
      out << "      " << j << " " << pi->params_string[j] << std::endl;
    }
  }
  out.precision(prec);
  
  return out;
  
}
