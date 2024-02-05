#define arc_raytracing_cxx
#include "arc_raytracing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1D.h>
#include "DDSegmentation/BitFieldCoder.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/BitFieldCoder.h"
#include "DDRec/CellIDPositionConverter.h"
R__LOAD_LIBRARY(libDDCore)

TH1D * hThetaTruth = new TH1D("hThetaTruth","",100,0,0.1);

TH2D * hXYtruth = new TH2D("hXYtruth","",200,-10,10,200,-10,10);
TH2D * hXYreco  = new TH2D("hXYreco","", 200,-10,10,200,-10,10);

static const double clhep2root_mm = 0.1; // cm/mm


void arc_raytracing::Loop()
{
//   In a ROOT session, you can do:
//      root> .L arc_raytracing.C
//      root> arc_raytracing t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   
   

   
auto lcdd = &(dd4hep::Detector::getInstance());
//lcdd->manager().SetDefaultUnits(TGeoManager::EDefaultUnits::kG4Units); //  xx = kG4Units, kRootUnits
lcdd->fromCompact("./compact/arc_v0.xml");
auto de = lcdd->detector("ARC_DETECTORNAME");
auto sensorde = de.child("ARC_sensor");
auto sensor2world_matrix = sensorde.nominal().worldTransformation();

auto r = lcdd->readout("ARC_HITS");
auto s = r.segmentation();

auto pixel2xyz = [&](const dd4hep::DDSegmentation::CellID& cellID) -> TVector3 {
auto p = s.position( cellID );
double local_coord [3] = {p.x()/clhep2root_mm,p.y()/clhep2root_mm,0.};
double global_coord[3] = {0.,0.,0.};
sensor2world_matrix.LocalToMaster(local_coord, global_coord);
TVector3 pixel_position(global_coord);
return pixel_position;
};









   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      TVector3 pion_p( MCParticles_momentum_x[0],  MCParticles_momentum_y[0],  MCParticles_momentum_z[0]);
      for( int i = 1; i < kMaxMCParticles; ++i )
      {
		std::cout << ARC_HITS_time[i] << '\t' << ARC_HITS_EDep[i]<< '\t' << ARC_HITS_cellID[i] << std::endl;
		TVector3 photon_p( MCParticles_momentum_x[i],  MCParticles_momentum_y[i],  MCParticles_momentum_z[i]);
		if( 0 == photon_p.Mag() ) break;
		auto cerenkov_angle = pion_p.Angle(photon_p);
		//4std::cout << MCParticles_time[i] << '\t' << cerenkov_angle << std::endl;
		hThetaTruth->Fill(cerenkov_angle);
		TVector3 pixel_position = pixel2xyz(ARC_HITS_cellID[i]);
		
		std::cout << "\t" << MCParticles_endpoint_x[i+1] << "\t" << pixel_position.x() << "\t" << pixel_position.z() << std::endl;
		//hXY->Fill(x,y);
		hXYtruth->Fill(MCParticles_endpoint_x[i+1],MCParticles_endpoint_y[i+1]);
		hXYreco->Fill(pixel_position.x(),pixel_position.y());
      }
      
		//std::cout << std::endl;
      break;
      
      
      
   }
}
