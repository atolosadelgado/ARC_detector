#include "ARCalg.h"

// ROOT
#include "TGeoBoolNode.h"

// STL
#include <iostream>
#include <sstream>

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ARCalg constructor       ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file", {"default name for the collection"}),
ARCalg::ARCalg(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                        {
                            KeyValues("ARC_simhits", {""}),
                            KeyValues("HeaderName", {"EventHeader"}),
                        },
                        {
                            KeyValues("ARC_PID", {"ARC_PID"})
                        }
                       )
{
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode ARCalg::initialize() {


	//-----------------
	// Retrieve the subdetector
	std::string ARC_name(m_ARC_name.value());
	if ( 0 == m_geoSvc->getDetector()->detectors().count(ARC_name) )
	{
		ThrowException( "Detector <<" + ARC_name + ">> does not exist." );
	}

	dd4hep::DetElement ARC_det = m_geoSvc->getDetector()->detectors().at(ARC_name);

    ///////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  retrieve mirror properties  //////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////
    auto cellDE = ARC_det; // for case of cubic single cell of ARC
    dd4hep::DetElement mirrorDE = ARC_det.child("ARC_mirror");
    auto mirrorSolid = mirrorDE.solid();
    auto m = (TGeoCompositeShape*)mirrorSolid.access();

    /// Retrieve the center of the sphere
    //-- this is just the matrix for building the intersection, translation
    auto matrix_boolean = m->GetBoolNode()->GetRightMatrix();
    //-- matrix for placing the cell
    auto matrix_cellvol = cellDE.nominal().worldTransformation();
    double local_coord [3] = {0.,0.,0.};
    double global_coord[3] = {0.,0.,0.};
    (matrix_cellvol**matrix_boolean).LocalToMaster(local_coord, global_coord);
    /// global_coord contains now the XYZ global coordinates of the center of the mirror
    mirrorCenter = dd4hep::Position( global_coord[0], global_coord[1], global_coord[2]);

    /// Retrieve the inner radius of the spherical mirror
    TGeoSphere * spherical_mirror_shape = (TGeoSphere*)m->GetBoolNode()->GetRightShape();
    mirrorRadius = spherical_mirror_shape->GetRmin();

    ///////////////////////////////////////////////////////////////////////////////////

	//-----------------
	// Retrieve the readout associated with the detector element (subdetector)
	dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(ARC_name);
	if(not dch_sd.isValid() )
		ThrowException("No valid Sensitive Detector was found for detector <<" + ARC_name + ">>.");

	dd4hep::Readout dch_readout = dch_sd.readout();
	// set the cellID decoder
	m_decoder = dch_readout.idSpec().decoder();
	//-----------------


	std::stringstream ss;
	PrintConfiguration(ss);
	info() << ss.str().c_str() <<endmsg;

	return StatusCode::SUCCESS;
}



///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple<colltype_out>
ARCalg::operator()(const colltype_in& input_sim_hits,
                         const edm4hep::EventHeaderCollection&  /*headers*/) const {


	debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

	// Create the collections we are going to return
	colltype_out output_digi_hits;

    //loop over hit collection
	for (const auto& input_sim_hit : input_sim_hits) {
      dd4hep::DDSegmentation::CellID id = input_sim_hit.getCellID();
      std::cout << "New ARC hit with cell ID: " << id << std::endl;



      // emplace output objects
      // output_digi_hits.create(....);

	}// end loop over hit collection


	/////////////////////////////////////////////////////////////////
	return std::make_tuple<colltype_out>(std::move(output_digi_hits));
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode ARCalg::finalize() { return StatusCode::SUCCESS; }

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ThrowException       ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void ARCalg::ThrowException(std::string s) const {
	error() << s.c_str()  << endmsg;
	throw std::runtime_error(s);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void ARCalg::PrintConfiguration(std::ostream& io)
{
   io << "ARCalg will use the following components:\n";
   io << "\tGeometry Service: "                  << m_geoSvcName.value().c_str()           << "\n";
   io << "\tDetector name: "                     << m_ARC_name.value().c_str()             << "\n";
   io << "\t\t|--Volume bitfield: "              << m_decoder->fieldDescription().c_str()  << "\n";
   io << "\t\t|--Mirror center (cm): "           << Form("(%g, %g, %g)",
                                                          mirrorCenter.x(),
                                                               mirrorCenter.y(),
                                                                   mirrorCenter.z())       << "\n";
   io << "\t\t|--Mirror radius (cm): "           << mirrorRadius                           << "\n";
   return;
}
