#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"

using namespace std;

Int_t miniAOD_metFilters_;
float miniAOD_genMET_;
float miniAOD_genMETPhi_;
float miniAOD_pfMET_;
float miniAOD_pfMETPhi_;
float miniAOD_pfMET_T1JERUp_;
float miniAOD_pfMET_T1JERDo_;
float miniAOD_pfMET_T1JESUp_;
float miniAOD_pfMET_T1JESDo_;
float miniAOD_pfMET_T1MESUp_;
float miniAOD_pfMET_T1MESDo_;
float miniAOD_pfMET_T1EESUp_;
float miniAOD_pfMET_T1EESDo_;
float miniAOD_pfMET_T1PESUp_;
float miniAOD_pfMET_T1PESDo_;
float miniAOD_pfMET_T1TESUp_;
float miniAOD_pfMET_T1TESDo_;
float miniAOD_pfMET_T1UESUp_;
float miniAOD_pfMET_T1UESDo_;
float miniAOD_pfMET_T1TxyPhi_;
float miniAOD_pfMET_T1TxyPt_;
float miniAOD_pfMETPhi_T1JESUp_;
float miniAOD_pfMETPhi_T1JESDo_;
float miniAOD_pfMETPhi_T1UESUp_;
float miniAOD_pfMETPhi_T1UESDo_;

void lldjNtuple::branchesMiniAODMET(TTree* tree) {

  tree->Branch("miniAOD_genMET",           &miniAOD_genMET_);
  tree->Branch("miniAOD_genMETPhi",        &miniAOD_genMETPhi_);
  tree->Branch("miniAOD_metFilters",       &miniAOD_metFilters_);
  tree->Branch("miniAOD_pfMET",            &miniAOD_pfMET_);
  tree->Branch("miniAOD_pfMETPhi",         &miniAOD_pfMETPhi_);
  tree->Branch("miniAOD_pfMET_T1JERUp",    &miniAOD_pfMET_T1JERUp_);
  tree->Branch("miniAOD_pfMET_T1JERDo",    &miniAOD_pfMET_T1JERDo_);
  tree->Branch("miniAOD_pfMET_T1JESUp",    &miniAOD_pfMET_T1JESUp_);
  tree->Branch("miniAOD_pfMET_T1JESDo",    &miniAOD_pfMET_T1JESDo_);
  tree->Branch("miniAOD_pfMET_T1MESUp",    &miniAOD_pfMET_T1MESUp_);
  tree->Branch("miniAOD_pfMET_T1MESDo",    &miniAOD_pfMET_T1MESDo_);
  tree->Branch("miniAOD_pfMET_T1EESUp",    &miniAOD_pfMET_T1EESUp_);
  tree->Branch("miniAOD_pfMET_T1EESDo",    &miniAOD_pfMET_T1EESDo_);
  tree->Branch("miniAOD_pfMET_T1PESUp",    &miniAOD_pfMET_T1PESUp_);
  tree->Branch("miniAOD_pfMET_T1PESDo",    &miniAOD_pfMET_T1PESDo_);
  tree->Branch("miniAOD_pfMET_T1TESUp",    &miniAOD_pfMET_T1TESUp_);
  tree->Branch("miniAOD_pfMET_T1TESDo",    &miniAOD_pfMET_T1TESDo_);
  tree->Branch("miniAOD_pfMET_T1UESUp",    &miniAOD_pfMET_T1UESUp_);
  tree->Branch("miniAOD_pfMET_T1UESDo",    &miniAOD_pfMET_T1UESDo_);
  tree->Branch("miniAOD_pfMETPhi_T1JESUp", &miniAOD_pfMETPhi_T1JESUp_);
  tree->Branch("miniAOD_pfMETPhi_T1JESDo", &miniAOD_pfMETPhi_T1JESDo_);
  tree->Branch("miniAOD_pfMETPhi_T1UESUp", &miniAOD_pfMETPhi_T1UESUp_);
  tree->Branch("miniAOD_pfMETPhi_T1UESDo", &miniAOD_pfMETPhi_T1UESDo_);

}

void lldjNtuple::fillMiniAODMET(const edm::Event& e, const edm::EventSetup& es) {

  miniAOD_metFilters_ = 0;

  //if (addFilterInfoMINIAOD_) {
    string filterNamesToCheck[9] = {
      "Flag_HBHENoiseFilter",
      "Flag_HBHENoiseIsoFilter", 
      "Flag_globalSuperTightHalo2016Filter",
      "Flag_goodVertices",
      "Flag_eeBadScFilter",
      "Flag_EcalDeadCellTriggerPrimitiveFilter",
      "Flag_BadPFMuonFilter",
      "Flag_ecalBadCalibFilter",
      "Flag_BadChargedCandidateFilter"
    };

    edm::Handle<edm::TriggerResults> miniAODpatFilterResultsHandle;
    e.getByToken(miniAODpatTrgResultsLabel_, miniAODpatFilterResultsHandle);
    edm::TriggerResults const& miniAODpatFilterResults = *miniAODpatFilterResultsHandle;
    
    auto&& filterNames = e.triggerNames(miniAODpatFilterResults);

    // === the following lines allow us to find the filters stored in the event ! ===
    /*
    edm::TriggerNames const& triggerNames = e.triggerNames(miniAODpatFilterResults);
    for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
	  triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
      int triggerId = triggerNames.triggerIndex(*triggerName);
      if ( triggerId >= 0 && triggerId < (int)triggerNames.size() ) {
	std::string triggerDecision = ( miniAODpatFilterResultsHandle->accept(triggerId) ) ? "passed" : "failed";
	  
	std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
      }
    }
    */

    for (unsigned iF = 0; iF < 9; ++iF) {
      unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
      if ( index == filterNames.size() ) 
	LogDebug("METFilters") << filterNamesToCheck[iF] << " is missing, exiting";
      else {
	if ( !miniAODpatFilterResults.accept(index) ) {
	  miniAOD_metFilters_ += pow(2, iF+1);
	}
      }
    }
//  } 
  
  edm::Handle<edm::View<pat::MET> > miniAODpfMETHandle;
  e.getByToken(miniAODpfMETlabel_, miniAODpfMETHandle);

  miniAOD_genMET_    = -99;
  miniAOD_genMETPhi_ = -99;
  miniAOD_pfMET_     = -99;
  miniAOD_pfMETPhi_  = -99;

  if (miniAODpfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET     = &(miniAODpfMETHandle->front());
    miniAOD_pfMET_    = pfMET->et();
    miniAOD_pfMETPhi_ = pfMET->phi();
    
    // Type1MET uncertainties =======================================
    miniAOD_pfMET_T1JERUp_ = pfMET->shiftedPt(pat::MET::JetResUp);
    miniAOD_pfMET_T1JERDo_ = pfMET->shiftedPt(pat::MET::JetResDown);
    miniAOD_pfMET_T1JESUp_ = pfMET->shiftedPt(pat::MET::JetEnUp);
    miniAOD_pfMET_T1JESDo_ = pfMET->shiftedPt(pat::MET::JetEnDown);
    miniAOD_pfMET_T1MESUp_ = pfMET->shiftedPt(pat::MET::MuonEnUp);
    miniAOD_pfMET_T1MESDo_ = pfMET->shiftedPt(pat::MET::MuonEnDown);
    miniAOD_pfMET_T1EESUp_ = pfMET->shiftedPt(pat::MET::ElectronEnUp);
    miniAOD_pfMET_T1EESDo_ = pfMET->shiftedPt(pat::MET::ElectronEnDown);
    miniAOD_pfMET_T1PESUp_ = pfMET->shiftedPt(pat::MET::PhotonEnUp);
    miniAOD_pfMET_T1PESDo_ = pfMET->shiftedPt(pat::MET::PhotonEnDown);
    miniAOD_pfMET_T1TESUp_ = pfMET->shiftedPt(pat::MET::TauEnUp);
    miniAOD_pfMET_T1TESDo_ = pfMET->shiftedPt(pat::MET::TauEnDown);
    miniAOD_pfMET_T1UESUp_ = pfMET->shiftedPt(pat::MET::UnclusteredEnUp);
    miniAOD_pfMET_T1UESDo_ = pfMET->shiftedPt(pat::MET::UnclusteredEnDown);

    miniAOD_pfMETPhi_T1JESUp_ = pfMET->shiftedPhi(pat::MET::JetEnUp);
    miniAOD_pfMETPhi_T1JESDo_ = pfMET->shiftedPhi(pat::MET::JetEnDown);
    miniAOD_pfMETPhi_T1UESUp_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnUp);
    miniAOD_pfMETPhi_T1UESDo_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnDown);
    
    if (!e.isRealData()) {
      miniAOD_genMET_    = pfMET->genMET()->et();
      miniAOD_genMETPhi_ = pfMET->genMET()->phi();
    }

  } 

}
