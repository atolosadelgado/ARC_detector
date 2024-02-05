//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  5 10:41:01 2024 by ROOT version 6.28/06
// from TTree events/events data tree
// found on file: arcsim_e+_50GeV_key4hep.root
//////////////////////////////////////////////////////////

#ifndef arc_raytracing_h
#define arc_raytracing_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "edm4hep/SimTrackerHitData.h"
#include "vector"
#include "podio/ObjectID.h"
#include "edm4hep/EventHeaderData.h"
#include "edm4hep/MCParticleData.h"
#include "podio/GenericParameters.h"

class arc_raytracing {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxARC_HITS = 203;
   static constexpr Int_t kMax_ARC_HITS_MCParticle = 203;
   static constexpr Int_t kMaxEventHeader = 1;
   static constexpr Int_t kMaxMCParticles = 208;
   static constexpr Int_t kMax_MCParticles_parents = 207;
   static constexpr Int_t kMax_MCParticles_daughters = 207;
   static constexpr Int_t kMax_intMap = 1;
   static constexpr Int_t kMax_floatMap = 1;
   static constexpr Int_t kMax_stringMap = 1;
   static constexpr Int_t kMax_doubleMap = 1;

   // Declaration of leaf types
   Int_t           ARC_HITS_;
   ULong_t         ARC_HITS_cellID[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_EDep[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_time[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_pathLength[kMaxARC_HITS];   //[ARC_HITS_]
   Int_t           ARC_HITS_quality[kMaxARC_HITS];   //[ARC_HITS_]
   Double_t        ARC_HITS_position_x[kMaxARC_HITS];   //[ARC_HITS_]
   Double_t        ARC_HITS_position_y[kMaxARC_HITS];   //[ARC_HITS_]
   Double_t        ARC_HITS_position_z[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_momentum_x[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_momentum_y[kMaxARC_HITS];   //[ARC_HITS_]
   Float_t         ARC_HITS_momentum_z[kMaxARC_HITS];   //[ARC_HITS_]
   Int_t           _ARC_HITS_MCParticle_;
   Int_t           _ARC_HITS_MCParticle_index[kMax_ARC_HITS_MCParticle];   //[_ARC_HITS_MCParticle_]
   UInt_t          _ARC_HITS_MCParticle_collectionID[kMax_ARC_HITS_MCParticle];   //[_ARC_HITS_MCParticle_]
   Int_t           EventHeader_;
   Int_t           EventHeader_eventNumber[kMaxEventHeader];   //[EventHeader_]
   Int_t           EventHeader_runNumber[kMaxEventHeader];   //[EventHeader_]
   ULong_t         EventHeader_timeStamp[kMaxEventHeader];   //[EventHeader_]
   Float_t         EventHeader_weight[kMaxEventHeader];   //[EventHeader_]
   Int_t           MCParticles_;
   Int_t           MCParticles_PDG[kMaxMCParticles];   //[MCParticles_]
   Int_t           MCParticles_generatorStatus[kMaxMCParticles];   //[MCParticles_]
   Int_t           MCParticles_simulatorStatus[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_charge[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_time[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_mass[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_vertex_x[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_vertex_y[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_vertex_z[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_endpoint_x[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_endpoint_y[kMaxMCParticles];   //[MCParticles_]
   Double_t        MCParticles_endpoint_z[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentum_x[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentum_y[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentum_z[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentumAtEndpoint_x[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentumAtEndpoint_y[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_momentumAtEndpoint_z[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_spin_x[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_spin_y[kMaxMCParticles];   //[MCParticles_]
   Float_t         MCParticles_spin_z[kMaxMCParticles];   //[MCParticles_]
   Int_t           MCParticles_colorFlow_a[kMaxMCParticles];   //[MCParticles_]
   Int_t           MCParticles_colorFlow_b[kMaxMCParticles];   //[MCParticles_]
   UInt_t          MCParticles_parents_begin[kMaxMCParticles];   //[MCParticles_]
   UInt_t          MCParticles_parents_end[kMaxMCParticles];   //[MCParticles_]
   UInt_t          MCParticles_daughters_begin[kMaxMCParticles];   //[MCParticles_]
   UInt_t          MCParticles_daughters_end[kMaxMCParticles];   //[MCParticles_]
   Int_t           _MCParticles_parents_;
   Int_t           _MCParticles_parents_index[kMax_MCParticles_parents];   //[_MCParticles_parents_]
   UInt_t          _MCParticles_parents_collectionID[kMax_MCParticles_parents];   //[_MCParticles_parents_]
   Int_t           _MCParticles_daughters_;
   Int_t           _MCParticles_daughters_index[kMax_MCParticles_daughters];   //[_MCParticles_daughters_]
   UInt_t          _MCParticles_daughters_collectionID[kMax_MCParticles_daughters];   //[_MCParticles_daughters_]
 //podio::GenericParameters *PARAMETERS;
   Int_t           _intMap_;
   string          _intMap_first[kMax_intMap];
   vector<int>     _intMap_second[kMax_intMap];
   Int_t           _floatMap_;
   string          _floatMap_first[kMax_floatMap];
   vector<float>   _floatMap_second[kMax_floatMap];
   Int_t           _stringMap_;
   string          _stringMap_first[kMax_stringMap];
   vector<string>  _stringMap_second[kMax_stringMap];
   Int_t           _doubleMap_;
   string          _doubleMap_first[kMax_doubleMap];
   vector<double>  _doubleMap_second[kMax_doubleMap];

   // List of branches
   TBranch        *b_ARC_HITS_;   //!
   TBranch        *b_ARC_HITS_cellID;   //!
   TBranch        *b_ARC_HITS_EDep;   //!
   TBranch        *b_ARC_HITS_time;   //!
   TBranch        *b_ARC_HITS_pathLength;   //!
   TBranch        *b_ARC_HITS_quality;   //!
   TBranch        *b_ARC_HITS_position_x;   //!
   TBranch        *b_ARC_HITS_position_y;   //!
   TBranch        *b_ARC_HITS_position_z;   //!
   TBranch        *b_ARC_HITS_momentum_x;   //!
   TBranch        *b_ARC_HITS_momentum_y;   //!
   TBranch        *b_ARC_HITS_momentum_z;   //!
   TBranch        *b__ARC_HITS_MCParticle_;   //!
   TBranch        *b__ARC_HITS_MCParticle_index;   //!
   TBranch        *b__ARC_HITS_MCParticle_collectionID;   //!
   TBranch        *b_EventHeader_;   //!
   TBranch        *b_EventHeader_eventNumber;   //!
   TBranch        *b_EventHeader_runNumber;   //!
   TBranch        *b_EventHeader_timeStamp;   //!
   TBranch        *b_EventHeader_weight;   //!
   TBranch        *b_MCParticles_;   //!
   TBranch        *b_MCParticles_PDG;   //!
   TBranch        *b_MCParticles_generatorStatus;   //!
   TBranch        *b_MCParticles_simulatorStatus;   //!
   TBranch        *b_MCParticles_charge;   //!
   TBranch        *b_MCParticles_time;   //!
   TBranch        *b_MCParticles_mass;   //!
   TBranch        *b_MCParticles_vertex_x;   //!
   TBranch        *b_MCParticles_vertex_y;   //!
   TBranch        *b_MCParticles_vertex_z;   //!
   TBranch        *b_MCParticles_endpoint_x;   //!
   TBranch        *b_MCParticles_endpoint_y;   //!
   TBranch        *b_MCParticles_endpoint_z;   //!
   TBranch        *b_MCParticles_momentum_x;   //!
   TBranch        *b_MCParticles_momentum_y;   //!
   TBranch        *b_MCParticles_momentum_z;   //!
   TBranch        *b_MCParticles_momentumAtEndpoint_x;   //!
   TBranch        *b_MCParticles_momentumAtEndpoint_y;   //!
   TBranch        *b_MCParticles_momentumAtEndpoint_z;   //!
   TBranch        *b_MCParticles_spin_x;   //!
   TBranch        *b_MCParticles_spin_y;   //!
   TBranch        *b_MCParticles_spin_z;   //!
   TBranch        *b_MCParticles_colorFlow_a;   //!
   TBranch        *b_MCParticles_colorFlow_b;   //!
   TBranch        *b_MCParticles_parents_begin;   //!
   TBranch        *b_MCParticles_parents_end;   //!
   TBranch        *b_MCParticles_daughters_begin;   //!
   TBranch        *b_MCParticles_daughters_end;   //!
   TBranch        *b__MCParticles_parents_;   //!
   TBranch        *b__MCParticles_parents_index;   //!
   TBranch        *b__MCParticles_parents_collectionID;   //!
   TBranch        *b__MCParticles_daughters_;   //!
   TBranch        *b__MCParticles_daughters_index;   //!
   TBranch        *b__MCParticles_daughters_collectionID;   //!
   TBranch        *b_PARAMETERS__intMap_;   //!
   TBranch        *b__intMap_first;   //!
   TBranch        *b__intMap_second;   //!
   TBranch        *b_PARAMETERS__floatMap_;   //!
   TBranch        *b__floatMap_first;   //!
   TBranch        *b__floatMap_second;   //!
   TBranch        *b_PARAMETERS__stringMap_;   //!
   TBranch        *b__stringMap_first;   //!
   TBranch        *b__stringMap_second;   //!
   TBranch        *b_PARAMETERS__doubleMap_;   //!
   TBranch        *b__doubleMap_first;   //!
   TBranch        *b__doubleMap_second;   //!

   arc_raytracing(TTree *tree=0);
   virtual ~arc_raytracing();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef arc_raytracing_cxx
arc_raytracing::arc_raytracing(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("arcsim_e+_50GeV_key4hep.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("arcsim_e+_50GeV_key4hep.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
   Loop();
}

arc_raytracing::~arc_raytracing()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t arc_raytracing::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t arc_raytracing::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void arc_raytracing::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ARC_HITS", &ARC_HITS_, &b_ARC_HITS_);
   fChain->SetBranchAddress("ARC_HITS.cellID", ARC_HITS_cellID, &b_ARC_HITS_cellID);
   fChain->SetBranchAddress("ARC_HITS.EDep", ARC_HITS_EDep, &b_ARC_HITS_EDep);
   fChain->SetBranchAddress("ARC_HITS.time", ARC_HITS_time, &b_ARC_HITS_time);
   fChain->SetBranchAddress("ARC_HITS.pathLength", ARC_HITS_pathLength, &b_ARC_HITS_pathLength);
   fChain->SetBranchAddress("ARC_HITS.quality", ARC_HITS_quality, &b_ARC_HITS_quality);
   fChain->SetBranchAddress("ARC_HITS.position.x", ARC_HITS_position_x, &b_ARC_HITS_position_x);
   fChain->SetBranchAddress("ARC_HITS.position.y", ARC_HITS_position_y, &b_ARC_HITS_position_y);
   fChain->SetBranchAddress("ARC_HITS.position.z", ARC_HITS_position_z, &b_ARC_HITS_position_z);
   fChain->SetBranchAddress("ARC_HITS.momentum.x", ARC_HITS_momentum_x, &b_ARC_HITS_momentum_x);
   fChain->SetBranchAddress("ARC_HITS.momentum.y", ARC_HITS_momentum_y, &b_ARC_HITS_momentum_y);
   fChain->SetBranchAddress("ARC_HITS.momentum.z", ARC_HITS_momentum_z, &b_ARC_HITS_momentum_z);
   fChain->SetBranchAddress("_ARC_HITS_MCParticle", &_ARC_HITS_MCParticle_, &b__ARC_HITS_MCParticle_);
   fChain->SetBranchAddress("_ARC_HITS_MCParticle.index", _ARC_HITS_MCParticle_index, &b__ARC_HITS_MCParticle_index);
   fChain->SetBranchAddress("_ARC_HITS_MCParticle.collectionID", _ARC_HITS_MCParticle_collectionID, &b__ARC_HITS_MCParticle_collectionID);
   fChain->SetBranchAddress("EventHeader", &EventHeader_, &b_EventHeader_);
   fChain->SetBranchAddress("EventHeader.eventNumber", EventHeader_eventNumber, &b_EventHeader_eventNumber);
   fChain->SetBranchAddress("EventHeader.runNumber", EventHeader_runNumber, &b_EventHeader_runNumber);
   fChain->SetBranchAddress("EventHeader.timeStamp", EventHeader_timeStamp, &b_EventHeader_timeStamp);
   fChain->SetBranchAddress("EventHeader.weight", EventHeader_weight, &b_EventHeader_weight);
   fChain->SetBranchAddress("MCParticles", &MCParticles_, &b_MCParticles_);
   fChain->SetBranchAddress("MCParticles.PDG", MCParticles_PDG, &b_MCParticles_PDG);
   fChain->SetBranchAddress("MCParticles.generatorStatus", MCParticles_generatorStatus, &b_MCParticles_generatorStatus);
   fChain->SetBranchAddress("MCParticles.simulatorStatus", MCParticles_simulatorStatus, &b_MCParticles_simulatorStatus);
   fChain->SetBranchAddress("MCParticles.charge", MCParticles_charge, &b_MCParticles_charge);
   fChain->SetBranchAddress("MCParticles.time", MCParticles_time, &b_MCParticles_time);
   fChain->SetBranchAddress("MCParticles.mass", MCParticles_mass, &b_MCParticles_mass);
   fChain->SetBranchAddress("MCParticles.vertex.x", MCParticles_vertex_x, &b_MCParticles_vertex_x);
   fChain->SetBranchAddress("MCParticles.vertex.y", MCParticles_vertex_y, &b_MCParticles_vertex_y);
   fChain->SetBranchAddress("MCParticles.vertex.z", MCParticles_vertex_z, &b_MCParticles_vertex_z);
   fChain->SetBranchAddress("MCParticles.endpoint.x", MCParticles_endpoint_x, &b_MCParticles_endpoint_x);
   fChain->SetBranchAddress("MCParticles.endpoint.y", MCParticles_endpoint_y, &b_MCParticles_endpoint_y);
   fChain->SetBranchAddress("MCParticles.endpoint.z", MCParticles_endpoint_z, &b_MCParticles_endpoint_z);
   fChain->SetBranchAddress("MCParticles.momentum.x", MCParticles_momentum_x, &b_MCParticles_momentum_x);
   fChain->SetBranchAddress("MCParticles.momentum.y", MCParticles_momentum_y, &b_MCParticles_momentum_y);
   fChain->SetBranchAddress("MCParticles.momentum.z", MCParticles_momentum_z, &b_MCParticles_momentum_z);
   fChain->SetBranchAddress("MCParticles.momentumAtEndpoint.x", MCParticles_momentumAtEndpoint_x, &b_MCParticles_momentumAtEndpoint_x);
   fChain->SetBranchAddress("MCParticles.momentumAtEndpoint.y", MCParticles_momentumAtEndpoint_y, &b_MCParticles_momentumAtEndpoint_y);
   fChain->SetBranchAddress("MCParticles.momentumAtEndpoint.z", MCParticles_momentumAtEndpoint_z, &b_MCParticles_momentumAtEndpoint_z);
   fChain->SetBranchAddress("MCParticles.spin.x", MCParticles_spin_x, &b_MCParticles_spin_x);
   fChain->SetBranchAddress("MCParticles.spin.y", MCParticles_spin_y, &b_MCParticles_spin_y);
   fChain->SetBranchAddress("MCParticles.spin.z", MCParticles_spin_z, &b_MCParticles_spin_z);
   fChain->SetBranchAddress("MCParticles.colorFlow.a", MCParticles_colorFlow_a, &b_MCParticles_colorFlow_a);
   fChain->SetBranchAddress("MCParticles.colorFlow.b", MCParticles_colorFlow_b, &b_MCParticles_colorFlow_b);
   fChain->SetBranchAddress("MCParticles.parents_begin", MCParticles_parents_begin, &b_MCParticles_parents_begin);
   fChain->SetBranchAddress("MCParticles.parents_end", MCParticles_parents_end, &b_MCParticles_parents_end);
   fChain->SetBranchAddress("MCParticles.daughters_begin", MCParticles_daughters_begin, &b_MCParticles_daughters_begin);
   fChain->SetBranchAddress("MCParticles.daughters_end", MCParticles_daughters_end, &b_MCParticles_daughters_end);
   fChain->SetBranchAddress("_MCParticles_parents", &_MCParticles_parents_, &b__MCParticles_parents_);
   fChain->SetBranchAddress("_MCParticles_parents.index", _MCParticles_parents_index, &b__MCParticles_parents_index);
   fChain->SetBranchAddress("_MCParticles_parents.collectionID", _MCParticles_parents_collectionID, &b__MCParticles_parents_collectionID);
   fChain->SetBranchAddress("_MCParticles_daughters", &_MCParticles_daughters_, &b__MCParticles_daughters_);
   fChain->SetBranchAddress("_MCParticles_daughters.index", _MCParticles_daughters_index, &b__MCParticles_daughters_index);
   fChain->SetBranchAddress("_MCParticles_daughters.collectionID", _MCParticles_daughters_collectionID, &b__MCParticles_daughters_collectionID);
   fChain->SetBranchAddress("_intMap", &_intMap_, &b_PARAMETERS__intMap_);
   fChain->SetBranchAddress("_intMap.first", &_intMap_first, &b__intMap_first);
   fChain->SetBranchAddress("_intMap.second", &_intMap_second, &b__intMap_second);
   fChain->SetBranchAddress("_floatMap", &_floatMap_, &b_PARAMETERS__floatMap_);
   fChain->SetBranchAddress("_floatMap.first", &_floatMap_first, &b__floatMap_first);
   fChain->SetBranchAddress("_floatMap.second", &_floatMap_second, &b__floatMap_second);
   fChain->SetBranchAddress("_stringMap", &_stringMap_, &b_PARAMETERS__stringMap_);
   fChain->SetBranchAddress("_stringMap.first", &_stringMap_first, &b__stringMap_first);
   fChain->SetBranchAddress("_stringMap.second", &_stringMap_second, &b__stringMap_second);
   fChain->SetBranchAddress("_doubleMap", &_doubleMap_, &b_PARAMETERS__doubleMap_);
   fChain->SetBranchAddress("_doubleMap.first", &_doubleMap_first, &b__doubleMap_first);
   fChain->SetBranchAddress("_doubleMap.second", &_doubleMap_second, &b__doubleMap_second);
   Notify();
}

Bool_t arc_raytracing::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void arc_raytracing::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t arc_raytracing::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef arc_raytracing_cxx
