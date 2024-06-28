#==========================================================================
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

IF(NOT TARGET DD4hep::DDCore)
  find_package ( DD4hep REQUIRED )
  include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
  include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )
  dd4hep_configure_output()
ENDIF()


dd4hep_configure_output()

dd4hep_add_plugin(ARCCELL SOURCES src/ARC_geo_cell0.cpp
  USES DD4hep::DDCore DD4hep::DDRec DD4hep::DDCond ROOT::Core ROOT::Geom ROOT::GenVector ROOT::MathCore)
install(TARGETS ARCCELL LIBRARY DESTINATION lib)

dd4hep_add_plugin(ARCBARREL SOURCES src/ARC_barrel_geo.cpp
  USES DD4hep::DDCore DD4hep::DDRec DD4hep::DDCond ROOT::Core ROOT::Geom ROOT::GenVector ROOT::MathCore)
install(TARGETS ARCBARREL LIBRARY DESTINATION lib)

dd4hep_add_plugin(ARCENDCAP SOURCES src/ARC_endcap_geo.cpp
  USES DD4hep::DDCore DD4hep::DDRec DD4hep::DDCond ROOT::Core ROOT::Geom ROOT::GenVector ROOT::MathCore)
install(TARGETS ARCENDCAP LIBRARY DESTINATION lib)

# Gaudi Algorithm for ARC
#
# dependencies
# For some reason, the find_package with DD4hep must be the first one
find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)
# dd4hep_set_compiler_flags() expects DD4hep_SET_RPATH to be set to ON
# otherwise it will not set the rpath when installing
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()
find_package(ROOT COMPONENTS RIO Tree MathCore MathMore)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(Gaudi)
find_package(GSL)
#-------------------
if(Gaudi_FOUND)
  add_subdirectory(ARCalg)
endif()
