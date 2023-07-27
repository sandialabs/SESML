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
#include <utility>

#include "EOSFactory.H"

//
// include files for all the models
//
#include "IAPWS-IF97.H"
#include "IAPWS95.H"

//
// Constructor
//
EOSParsingData::EOSParsingData( EOSCreator creatorfunction, const EOSParamInfo * paraminfo, const int nparaminfo )
  : c(creatorfunction)
{
  for (int i=0;i<nparaminfo;i++) {
    if (paraminfo[i].type == "DParam") dm[paraminfo[i].name] = paraminfo[i].index;
    else if (paraminfo[i].type == "IParam") im[paraminfo[i].name] = paraminfo[i].index;
    else if (paraminfo[i].type == "SParam") sm[paraminfo[i].name] = paraminfo[i].index;
    else throw std::runtime_error("EOSParsingData found invalid parameter type: "+paraminfo[i].type);
  }
}

//
// Constructor
//
EOSFactory::EOSFactory()
{
  if (!ParseMapInitialized) initializeParseMap();
}

//
// Destructor
//
EOSFactory::~EOSFactory()
{

}

bool EOSFactory::ParseMapInitialized = false;

EOSParseMap & EOSFactory::getParseMap()
{
  static EOSParseMap pm;
  return pm;
}

void EOSFactory::initializeParseMap()
{
  EOSParseMap & m = getParseMap();

  if (!m.insert(std::pair<std::string, EOSParsingData>(IAPWS95::name,
                                                       EOSParsingData( IAPWS95::createModel,
                                                                       IAPWS95::ParamList,
                                                                       IAPWS95::ParamListSize ))).second) {
    throw std::logic_error("EOSFactory::initializeParseMap -- duplicate model found in map when trying to insert "+IAPWS95::name);
  }

  if (!m.insert(std::pair<std::string, EOSParsingData>(IAPWSIF97::name,
                                                       EOSParsingData( IAPWSIF97::createModel,
                                                                       IAPWSIF97::ParamList,
                                                                       IAPWSIF97::ParamListSize ))).second) {
    throw std::logic_error("EOSFactory::initializeParseMap -- duplicate model found in map when trying to insert "+IAPWSIF97::name);
  }

  ParseMapInitialized = true;
}

//
// Return a newly allocated model as given by the name
//
EOSModel * EOSFactory::getModel( const Parameters & params,
                                 const std::string name,
                                 const EOSData & init_data )
{
  Parameters::const_iterator p = params.find(name);
  if (p == params.end()) {
    throw std::runtime_error("EOSFactory::getModel -- EOSParam for \""+name+"\" not found in params");
  }
  const std::string & modname = p->second->model_name;
  
  EOSParseMap::const_iterator m = getParseMap().find(modname);
  if (m == getParseMap().end()) {
    throw std::runtime_error("EOSFactory::getModel -- model type \""+modname+"\" not found in parse map");
  }
  
  return m->second.c(params,name,init_data);

}

//
// Return parameter counts for a model type
//
void EOSFactory::getParamCounts( const std::string modname,
                                 std::vector<int> & pnum )
{
  EOSParseMap::const_iterator m = getParseMap().find(modname);
  if (m == getParseMap().end()) {
    throw std::runtime_error("EOSFactory::getParamCounts -- model type \""+modname+"\" not found in parse map");
  }

  pnum[0] = m->second.dm.size();
  pnum[1] = m->second.im.size();
  pnum[2] = m->second.sm.size();

}

//
// Return parameter location for a model type
//
int EOSFactory::getParamLocation( const std::string modname,
                                  const std::string type,
                                  const std::string name )
{
  EOSParseMap::const_iterator m = getParseMap().find(modname);
  if (m == getParseMap().end()) {
    throw std::runtime_error("EOSFactory::getParamLocation -- model type \""+modname+"\" not found in parse map");
  }

  if (type == "DParam") {
    std::map<std::string,int>::const_iterator it = m->second.dm.find(name);
    if (it == m->second.dm.end()) {
      throw std::runtime_error("EOSFactory::getParamLocation -- double parameter "+name+" not found in model "+modname);
    }
    return it->second;
  }
  else if (type == "IParam") {
    std::map<std::string,int>::const_iterator it = m->second.im.find(name);
    if (it == m->second.im.end()) {
      throw std::runtime_error("EOSFactory::getParamLocation -- integer parameter "+name+" not found in model "+modname);
    }
    return it->second;
  }
  else if (type == "SParam") {
    std::map<std::string,int>::const_iterator it = m->second.sm.find(name);
    if (it == m->second.sm.end()) {
      throw std::runtime_error("EOSFactory::getParamLocation -- string parameter "+name+" not found in model "+modname);
    }
    return it->second;
  }
  else {
    throw std::runtime_error("EOSFactory::getparamLocation -- invalid parameter type "+type);
  }

}

//
// Return a pointer to the parameters entry for the given name
//
const EOSParam * EOSFactory::getParams( const Parameters & params,
                                        const std::string name )
{
  Parameters::const_iterator p = params.find(name);
  if (p == params.end()) {
    throw std::runtime_error("EOSParam for \""+name+"\" not found");
  }
  return p->second;
}
