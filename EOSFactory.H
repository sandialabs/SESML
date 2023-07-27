/*--------------------------------------------------------------------*/
/*  Copyright (2013) Sandia Corporation. Under the terms of Contract  */
/*  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government    */
/*  retains certain rights in this software.                          */
/*                                                                    */
/*  For full terms of distribution see the accompanying LICENSE file. */
/*  To distribute this file outside of SESML, substitue the full text */
/*  of the LICENSE file for this notice.                              */
/*--------------------------------------------------------------------*/

#ifndef EOS_FACTORY_H
#define EOS_FACTORY_H

#include <vector>
#include <string>

#include "EOSModel.H"
#include "EOSData.H"
#include "EOSParam.H"

typedef EOSModel * (*EOSCreator)( const Parameters & params,
                                  const std::string name,
                                  const EOSData & init_data );

typedef struct {
  std::string name;
  std::string type;
  int index;
} EOSParamInfo;

struct EOSParsingData
{
  EOSParsingData( EOSCreator creatorfunction,
                  const EOSParamInfo * paraminfo,
                  const int nparaminfo );

  EOSCreator c;
  std::map<std::string,int> dm;
  std::map<std::string,int> im;
  std::map<std::string,int> sm;
};

typedef std::map<std::string, EOSParsingData> EOSParseMap;

class EOSFactory
{
public:
  //
  // Constructor
  //
  EOSFactory();

  //
  // Destructor
  //
  virtual ~EOSFactory();

  //
  // Factory method for models
  //
  EOSModel * getModel( const Parameters & params,
                       const std::string name )
  {
    EOSData d(0,0);
    return getModel(params,name,d);
  }

  //
  // Factory method to allocate a new model with input data
  //
  EOSModel * getModel( const Parameters & params,
                       const std::string name,
                       const EOSData & init_data );

  //
  // Method to return parameter counts for a model
  //
  void getParamCounts( const std::string modname,
                       std::vector<int> & pnum );

  //
  // Method to get location of a named parameter in its type array
  //
  int getParamLocation( const std::string modname,
                        const std::string type,
                        const std::string name );

protected:

  //
  // Method to retrieve parameters for a model
  //
  const EOSParam * getParams( const Parameters & params,
                              const std::string name );

private:

  //
  // return a static map to be used for 
  //
  static EOSParseMap & getParseMap();

  //
  // flag for map initialization
  //
  static bool ParseMapInitialized;

  //
  // initialize the parsing map
  //
  void initializeParseMap();

};

#endif
