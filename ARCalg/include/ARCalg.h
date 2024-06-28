/** ======= ARCalg ==========
 * Template Gaudi Algorithm for ARC reconstruction
 * At the moment take SimTrackerHit in ARC. Momentum can be obtained from MCParticle, to be replaced later by reconstructed track
 *
 *
 * @author Serena Pezzulo, Roberta Cardinale, Alvaro Tolosa-Delgado
 * @date   2024-06
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * <h4>Output</h4>
 * Processor produces collection of ParticleID<br>
 * @param ARC_simhits The name of input collection, type edm4hep::SimTrackerHitCollection <br>
 * (default name empty) <br>
 * @param ARC_PID The name of output collection, type edm4hep::ParticleID <br>
 * (default name "ARC_PID") <br>
 * @param GeoSvcName Geometry service name <br>
 * (default value GeoSvc)
 * @param ARC_name ARC subdetector name <br>
 * (default value ARC_DETECTORNAME) <br>
 * <br>
 */

#ifndef ARCALG_H
#define ARCALG_H

#include "Gaudi/Property.h"
#include "GaudiAlg/Transformer.h"
#include "GaudiKernel/RndmGenerators.h"

// EDM4HEP
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/ParticleIDCollection.h"
#include "edm4hep/EventHeaderCollection.h"

#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// DD4hep
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

#include <string>

using colltype_in  = edm4hep::SimTrackerHitCollection;
using colltype_out = edm4hep::ParticleIDCollection;

struct ARCalg final
    : k4FWCore::MultiTransformer<
          std::tuple<colltype_out>(
              const colltype_in&, const edm4hep::EventHeaderCollection&)> {
  ARCalg(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<colltype_out> operator()(
      const colltype_in& ,
      const edm4hep::EventHeaderCollection&   ) const override;

private:

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_ARC_name{this, "ARC_name", "ARC_DETECTORNAME", "Name of the ARC cell"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc>                        m_geoSvc;

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and throw exception
  void ThrowException(std::string s) const;

  /// Declare here variables to be initialized at when creating the algorithm and then used within the event loop
  dd4hep::Position mirrorCenter;
  double mirrorRadius;


};

DECLARE_COMPONENT(ARCalg)

#endif
