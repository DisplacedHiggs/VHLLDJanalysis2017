#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "FWCore/Common/interface/TriggerNames.h"
#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"
#include <iomanip>
#include <bitset>
using namespace std;

// (local) variables associated with tree branchesMiniAOD
Int_t       run_;
Long64_t    event_;
Int_t       lumis_;
Bool_t      isData_;
Int_t       nVtx_;
Int_t       nGoodVtx_;
Int_t       nTrksPV_;
Bool_t      isPVGood_;
float       vtx_;
float       vty_;
float       vtz_;
float       miniAODrho_;
float       miniAODrhoCentral_;

Int_t       nTruePU_;

void lldjNtuple::branchesMiniAODGlobalEvent(TTree* tree) {

  tree->Branch("miniAOD_run",        &run_);
  tree->Branch("miniAOD_event",      &event_);
  tree->Branch("miniAOD_lumis",      &lumis_);
  tree->Branch("miniAOD_isData",     &isData_);
  tree->Branch("miniAOD_nVtx",       &nVtx_);
  tree->Branch("miniAOD_nGoodVtx",   &nGoodVtx_);
  tree->Branch("miniAOD_nTrksPV",    &nTrksPV_);
  tree->Branch("miniAOD_isPVGood",   &isPVGood_);
  tree->Branch("miniAOD_vtx",        &vtx_);
  tree->Branch("miniAOD_vty",        &vty_);
  tree->Branch("miniAOD_vtz",        &vtz_);
  tree->Branch("miniAOD_rho",        &miniAODrho_);
  tree->Branch("miniAOD_rhoCentral", &miniAODrhoCentral_);

  tree->Branch("miniAOD_nTruePU",  &nTruePU_);

}

void lldjNtuple::fillMiniAODGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  //phoPrescale_.clear();

  edm::Handle<double> miniAODrhoHandle;
  e.getByToken(miniAODrhoLabel_, miniAODrhoHandle);

  edm::Handle<double> miniAODrhoCentralHandle;
  e.getByToken(miniAODrhoCentralLabel_, miniAODrhoCentralHandle);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  miniAODrho_    = *(miniAODrhoHandle.product());
  if (miniAODrhoCentralHandle.isValid()) miniAODrhoCentral_ = *(miniAODrhoCentralHandle.product());
  else miniAODrhoCentral_ = -99.;

  nTruePU_ = -1 ;
  if (!e.isRealData()) {
   edm::Handle<vector<PileupSummaryInfo> > miniAODpuInfoHandle;
   e.getByToken(miniAODpuCollection_, miniAODpuInfoHandle);
   if ( miniAODpuInfoHandle->size() > 0 ){
    nTruePU_ = miniAODpuInfoHandle->at(1).getTrueNumInteractions() ;
   }
  }

  edm::Handle<reco::VertexCollection> miniAODvtxHandle;
  e.getByToken(miniAODvtxLabel_, miniAODvtxHandle);

  nVtx_     = -1;
  nGoodVtx_ = -1;
  if (miniAODvtxHandle.isValid())
  {
   nVtx_     = 0;
   nGoodVtx_ = 0;

   for (vector<reco::Vertex>::const_iterator v = miniAODvtxHandle->begin(); v != miniAODvtxHandle->end(); ++v)
   {
    if (nVtx_ == 0)
    {
     nTrksPV_ = v->nTracks();
     vtx_     = v->x();
     vty_     = v->y();
     vtz_     = v->z();

     isPVGood_ = false;
     if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) isPVGood_ = true;
    }

    if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
    nVtx_++;
   }
  }
  else {edm::LogWarning("lldjNtuple") << "Primary vertices info not unavailable";}

}
