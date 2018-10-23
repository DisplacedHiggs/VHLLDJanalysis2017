// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"

// (local) variables associated with tree branchesMiniAOD

// triggers used in VHbb
// https://indico.cern.ch/event/655639/contributions/2670493/attachments/1497203/2337123/VHbb_Approval_28_07_17_v10.pdf#page=8

// trigger fired variables
ULong64_t   AOD_HLT_Ele23Loose_;
ULong64_t   AOD_HLT_Ele27Tight_;
ULong64_t   AOD_HLT_Ele17Ele12_;
ULong64_t   AOD_HLT_Ele23Ele12_;

ULong64_t   AOD_HLT_IsoMu22_;
ULong64_t   AOD_HLT_IsoTkMu22_;
ULong64_t   AOD_HLT_Mu17Mu8_;
ULong64_t   AOD_HLT_Mu17TkMu8_;

ULong64_t   AOD_HLT_Mu17Mu8_noDZ_;
ULong64_t   AOD_HLT_Mu17TkMu8_noDZ_;
ULong64_t   AOD_HLT_IsoMu24_;
ULong64_t   AOD_HLT_IsoTkMu24_;

ULong64_t   AOD_HLT_Mu8Ele23_;
ULong64_t   AOD_HLT_Mu23Ele12_;
ULong64_t   AOD_HLT_Mu12Ele23_DZ_;
ULong64_t   AOD_HLT_Mu23Ele12_DZ_;

ULong64_t   AOD_HLT_Photon90_;
ULong64_t   AOD_HLT_Photon120_;
ULong64_t   AOD_HLT_Photon175_;
ULong64_t   AOD_HLT_Photon165_HE10_;

// trigger is prescaled variables
ULong64_t   AOD_HLT_Ele23Loose_isPS_;
ULong64_t   AOD_HLT_Ele27Tight_isPS_;
ULong64_t   AOD_HLT_Ele17Ele12_isPS_;
ULong64_t   AOD_HLT_Ele23Ele12_isPS_;

ULong64_t   AOD_HLT_IsoMu22_isPS_;
ULong64_t   AOD_HLT_IsoTkMu22_isPS_;
ULong64_t   AOD_HLT_Mu17Mu8_isPS_;
ULong64_t   AOD_HLT_Mu17TkMu8_isPS_;

ULong64_t   AOD_HLT_Mu17Mu8_noDZ_isPS_;
ULong64_t   AOD_HLT_Mu17TkMu8_noDZ_isPS_;
ULong64_t   AOD_HLT_IsoMu24_isPS_;
ULong64_t   AOD_HLT_IsoTkMu24_isPS_;

ULong64_t   AOD_HLT_Mu8Ele23_isPS_;
ULong64_t   AOD_HLT_Mu23Ele12_isPS_;
ULong64_t   AOD_HLT_Mu12Ele23_DZ_isPS_;
ULong64_t   AOD_HLT_Mu23Ele12_DZ_isPS_;

ULong64_t   AOD_HLT_Photon90_isPS_;
ULong64_t   AOD_HLT_Photon120_isPS_;
ULong64_t   AOD_HLT_Photon175_isPS_;
ULong64_t   AOD_HLT_Photon165_HE10_isPS_;


void lldjNtuple::branchesMiniAODTrigger(TTree* tree){
  tree->Branch("AOD_HLT_Ele23Loose",       &AOD_HLT_Ele23Loose_) ;
  tree->Branch("AOD_HLT_Ele27Tight",       &AOD_HLT_Ele27Tight_) ;
  tree->Branch("AOD_HLT_Ele17Ele12",       &AOD_HLT_Ele17Ele12_) ;
  tree->Branch("AOD_HLT_Ele23Ele12",       &AOD_HLT_Ele23Ele12_) ;
  
  tree->Branch("AOD_HLT_IsoMu22"  ,        &AOD_HLT_IsoMu22_)   ;
  tree->Branch("AOD_HLT_IsoTkMu22",        &AOD_HLT_IsoTkMu22_) ;
  tree->Branch("AOD_HLT_Mu17Mu8"  ,        &AOD_HLT_Mu17Mu8_)   ;
  tree->Branch("AOD_HLT_Mu17TkMu8",        &AOD_HLT_Mu17TkMu8_) ;

  tree->Branch("AOD_HLT_Mu17Mu8_noDZ"  ,   &AOD_HLT_Mu17Mu8_noDZ_)   ;
  tree->Branch("AOD_HLT_Mu17TkMu8_noDZ",   &AOD_HLT_Mu17TkMu8_noDZ_) ;
  tree->Branch("AOD_HLT_IsoMu24"  ,        &AOD_HLT_IsoMu24_)   ;
  tree->Branch("AOD_HLT_IsoTkMu24",        &AOD_HLT_IsoTkMu24_) ;

  tree->Branch("AOD_HLT_Mu8Ele23"  ,       &AOD_HLT_Mu8Ele23_)     ;
  tree->Branch("AOD_HLT_Mu23Ele12"  ,      &AOD_HLT_Mu23Ele12_)    ;
  tree->Branch("AOD_HLT_Mu12Ele23_DZ"  ,   &AOD_HLT_Mu12Ele23_DZ_) ;
  tree->Branch("AOD_HLT_Mu23Ele12_DZ"  ,   &AOD_HLT_Mu23Ele12_DZ_) ;

  tree->Branch("AOD_HLT_Photon90",         &AOD_HLT_Photon90_);
  tree->Branch("AOD_HLT_Photon120",        &AOD_HLT_Photon120_);
  tree->Branch("AOD_HLT_Photon175",        &AOD_HLT_Photon175_);
  tree->Branch("AOD_HLT_Photon165_HE10",   &AOD_HLT_Photon165_HE10_);

  //prescale
  tree->Branch("AOD_HLT_Ele23Loose_isPS",       &AOD_HLT_Ele23Loose_isPS_) ;
  tree->Branch("AOD_HLT_Ele27Tight_isPS",       &AOD_HLT_Ele27Tight_isPS_) ;
  tree->Branch("AOD_HLT_Ele17Ele12_isPS",       &AOD_HLT_Ele17Ele12_isPS_) ;
  tree->Branch("AOD_HLT_Ele23Ele12_isPS",       &AOD_HLT_Ele23Ele12_isPS_) ;
  
  tree->Branch("AOD_HLT_IsoMu22_isPS",          &AOD_HLT_IsoMu22_isPS_)   ;
  tree->Branch("AOD_HLT_IsoTkMu22_isPS",        &AOD_HLT_IsoTkMu22_isPS_) ;
  tree->Branch("AOD_HLT_Mu17Mu8_isPS",          &AOD_HLT_Mu17Mu8_isPS_)   ;
  tree->Branch("AOD_HLT_Mu17TkMu8_isPS",        &AOD_HLT_Mu17TkMu8_isPS_) ;

  tree->Branch("AOD_HLT_Mu17Mu8_noDZ_isPS",     &AOD_HLT_Mu17Mu8_noDZ_isPS_)   ;
  tree->Branch("AOD_HLT_Mu17TkMu8_noDZ_isPS",   &AOD_HLT_Mu17TkMu8_noDZ_isPS_) ;
  tree->Branch("AOD_HLT_IsoMu24_isPS",          &AOD_HLT_IsoMu24_isPS_)   ;
  tree->Branch("AOD_HLT_IsoTkMu24_isPS",        &AOD_HLT_IsoTkMu24_isPS_) ;

  tree->Branch("AOD_HLT_Mu8Ele23_isPS",         &AOD_HLT_Mu8Ele23_isPS_)     ;
  tree->Branch("AOD_HLT_Mu23Ele12_isPS",        &AOD_HLT_Mu23Ele12_isPS_)    ;
  tree->Branch("AOD_HLT_Mu12Ele23_DZ_isPS",     &AOD_HLT_Mu12Ele23_DZ_isPS_) ;
  tree->Branch("AOD_HLT_Mu23Ele12_DZ_isPS",     &AOD_HLT_Mu23Ele12_DZ_isPS_) ;

  tree->Branch("AOD_HLT_Photon90_isPS",         &AOD_HLT_Photon90_isPS_);
  tree->Branch("AOD_HLT_Photon120_isPS",        &AOD_HLT_Photon120_isPS_);
  tree->Branch("AOD_HLT_Photon175_isPS",        &AOD_HLT_Photon175_isPS_);
  tree->Branch("AOD_HLT_Photon165_HE10_isPS",   &AOD_HLT_Photon165_HE10_isPS_);
}


void lldjNtuple::fillMiniAODTrigger(const edm::Event &e, const edm::EventSetup &es){
  
  AOD_HLT_Ele23Loose_       = 0;
  AOD_HLT_Ele27Tight_       = 0;
  AOD_HLT_Ele17Ele12_       = 0;
  AOD_HLT_Ele23Ele12_       = 0;
  
  AOD_HLT_IsoMu22_          = 0;
  AOD_HLT_IsoTkMu22_        = 0;
  AOD_HLT_Mu17Mu8_          = 0;
  AOD_HLT_Mu17TkMu8_        = 0;

  AOD_HLT_Mu17Mu8_noDZ_     = 0;
  AOD_HLT_Mu17TkMu8_noDZ_   = 0;
  AOD_HLT_IsoMu24_          = 0;
  AOD_HLT_IsoTkMu24_        = 0;

  AOD_HLT_Mu8Ele23_          = 0;
  AOD_HLT_Mu23Ele12_         = 0;
  AOD_HLT_Mu12Ele23_DZ_      = 0;
  AOD_HLT_Mu23Ele12_DZ_      = 0;

  AOD_HLT_Photon90_          = 0;
  AOD_HLT_Photon120_         = 0;
  AOD_HLT_Photon175_         = 0;
  AOD_HLT_Photon165_HE10_    = 0;  

  //prescale
  AOD_HLT_Ele23Loose_isPS_       = 0;
  AOD_HLT_Ele27Tight_isPS_       = 0;
  AOD_HLT_Ele17Ele12_isPS_       = 0;
  AOD_HLT_Ele23Ele12_isPS_       = 0;
  
  AOD_HLT_IsoMu22_isPS_          = 0;
  AOD_HLT_IsoTkMu22_isPS_        = 0;
  AOD_HLT_Mu17Mu8_isPS_          = 0;
  AOD_HLT_Mu17TkMu8_isPS_        = 0;

  AOD_HLT_Mu17Mu8_noDZ_isPS_     = 0;
  AOD_HLT_Mu17TkMu8_noDZ_isPS_   = 0;
  AOD_HLT_IsoMu24_isPS_          = 0;
  AOD_HLT_IsoTkMu24_isPS_        = 0;

  AOD_HLT_Mu8Ele23_isPS_          = 0;
  AOD_HLT_Mu23Ele12_isPS_         = 0;
  AOD_HLT_Mu12Ele23_DZ_isPS_      = 0;
  AOD_HLT_Mu23Ele12_DZ_isPS_      = 0;

  AOD_HLT_Photon90_isPS_          = 0;
  AOD_HLT_Photon120_isPS_         = 0;
  AOD_HLT_Photon175_isPS_         = 0;
  AOD_HLT_Photon165_HE10_isPS_    = 0;  


  edm::Handle<pat::TriggerObjectStandAloneCollection> miniAODtriggerHandle;
  e.getByToken(miniAODtriggerObjectsLabel_, miniAODtriggerHandle);

  edm::Handle<edm::TriggerResults> miniAODtrgResultsHandle;
  e.getByToken(miniAODtrgResultsLabel_, miniAODtrgResultsHandle);

  edm::Handle<std::string> miniAODfilterLabels_;
  e.getByLabel("slimmedPatTrigger:filterLabels", miniAODfilterLabels_);

//  if( !( e.getByToken(AODTriggerToken_, AODTriggerHandle_) && e.getByToken(AODTriggerEventToken_, AODTriggerEventHandle_) ) ){ return; }
//  
//  for(int i = 0; i < (int)hltConfig_.size(); i++){
//    string name = hltConfig_.triggerName(i);
//    //std::cout<<name<<"\n";
//    
//    int bitEle23Loose = -1;
//    if      (name.find("HLT_Ele23_WPLoose_Gsf_v1")  != string::npos) bitEle23Loose = 0  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v2")  != string::npos) bitEle23Loose = 1  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v3")  != string::npos) bitEle23Loose = 2  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v4")  != string::npos) bitEle23Loose = 3  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v5")  != string::npos) bitEle23Loose = 4  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v6")  != string::npos) bitEle23Loose = 5  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v7")  != string::npos) bitEle23Loose = 6  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v8")  != string::npos) bitEle23Loose = 7  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v9")  != string::npos) bitEle23Loose = 8  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v10") != string::npos) bitEle23Loose = 9  ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v11") != string::npos) bitEle23Loose = 10 ;
//    else if (name.find("HLT_Ele23_WPLoose_Gsf_v12") != string::npos) bitEle23Loose = 11 ;
//    
//    int bitEle27Tight = -1;
//    if      (name.find("HLT_Ele27_WPTight_Gsf_v1")  != string::npos) bitEle27Tight = 0 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v2")  != string::npos) bitEle27Tight = 1 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v3")  != string::npos) bitEle27Tight = 2 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v4")  != string::npos) bitEle27Tight = 3 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v5")  != string::npos) bitEle27Tight = 4 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v6")  != string::npos) bitEle27Tight = 5 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v7")  != string::npos) bitEle27Tight = 6 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v8")  != string::npos) bitEle27Tight = 7 ;
//    else if (name.find("HLT_Ele27_WPTight_Gsf_v9")  != string::npos) bitEle27Tight = 8 ;
//    
//    // Double Electron
//    int bitEle17Ele12 = -1;
//    if      (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1")  != string::npos) bitEle17Ele12 = 0 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2")  != string::npos) bitEle17Ele12 = 1 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3")  != string::npos) bitEle17Ele12 = 2 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4")  != string::npos) bitEle17Ele12 = 3 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5")  != string::npos) bitEle17Ele12 = 4 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6")  != string::npos) bitEle17Ele12 = 5 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7")  != string::npos) bitEle17Ele12 = 6 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8")  != string::npos) bitEle17Ele12 = 7 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9")  != string::npos) bitEle17Ele12 = 8 ;
//    else if (name.find("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10") != string::npos) bitEle17Ele12 = 9 ;
//    
//    int bitEle23Ele12 = -1;
//    if      (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1")  != string::npos) bitEle23Ele12 = 0 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2")  != string::npos) bitEle23Ele12 = 1 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3")  != string::npos) bitEle23Ele12 = 2 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4")  != string::npos) bitEle23Ele12 = 3 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5")  != string::npos) bitEle23Ele12 = 4 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6")  != string::npos) bitEle23Ele12 = 5 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7")  != string::npos) bitEle23Ele12 = 6 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8")  != string::npos) bitEle23Ele12 = 7 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v9")  != string::npos) bitEle23Ele12 = 8 ;
//    else if (name.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v10") != string::npos) bitEle23Ele12 = 9 ;
//    
//    // Single Iso Mu
//    int bitIsoMu22   = -1;
//    if      (name.find("HLT_IsoMu22_v1") != string::npos) bitIsoMu22 = 0 ;
//    else if (name.find("HLT_IsoMu22_v2") != string::npos) bitIsoMu22 = 1 ;
//    else if (name.find("HLT_IsoMu22_v3") != string::npos) bitIsoMu22 = 2 ;
//    else if (name.find("HLT_IsoMu22_v4") != string::npos) bitIsoMu22 = 3 ;
//    else if (name.find("HLT_IsoMu22_v5") != string::npos) bitIsoMu22 = 4 ;
//    else if (name.find("HLT_IsoMu22_v6") != string::npos) bitIsoMu22 = 5 ;
//    else if (name.find("HLT_IsoMu22_v7") != string::npos) bitIsoMu22 = 6 ;
//    
//    int bitIsoTkMu22 = -1;
//    if      (name.find("HLT_IsoTkMu22_v1")  != string::npos) bitIsoTkMu22 = 0 ;
//    else if (name.find("HLT_IsoTkMu22_v2")  != string::npos) bitIsoTkMu22 = 1 ;
//    else if (name.find("HLT_IsoTkMu22_v3")  != string::npos) bitIsoTkMu22 = 2 ;
//    else if (name.find("HLT_IsoTkMu22_v4")  != string::npos) bitIsoTkMu22 = 3 ;
//    else if (name.find("HLT_IsoTkMu22_v5")  != string::npos) bitIsoTkMu22 = 4 ;
//    else if (name.find("HLT_IsoTkMu22_v6")  != string::npos) bitIsoTkMu22 = 5 ;
//    else if (name.find("HLT_IsoTkMu22_v7")  != string::npos) bitIsoTkMu22 = 6 ;
//    else if (name.find("HLT_IsoTkMu22_v8")  != string::npos) bitIsoTkMu22 = 7 ;
//    else if (name.find("HLT_IsoTkMu22_v9")  != string::npos) bitIsoTkMu22 = 8 ;
//    else if (name.find("HLT_IsoTkMu22_v10") != string::npos) bitIsoTkMu22 = 9 ;
//
//    int bitIsoMu24   = -1;
//    if      (name.find("HLT_IsoMu24_v1") != string::npos) bitIsoMu24 = 0 ;
//    else if (name.find("HLT_IsoMu24_v2") != string::npos) bitIsoMu24 = 1 ;
//    else if (name.find("HLT_IsoMu24_v3") != string::npos) bitIsoMu24 = 2 ;
//    else if (name.find("HLT_IsoMu24_v4") != string::npos) bitIsoMu24 = 3 ;
//    
//    int bitIsoTkMu24 = -1;
//    if      (name.find("HLT_IsoTkMu24_v1")  != string::npos) bitIsoTkMu24 = 0 ;
//    else if (name.find("HLT_IsoTkMu24_v2")  != string::npos) bitIsoTkMu24 = 1 ;
//    else if (name.find("HLT_IsoTkMu24_v3")  != string::npos) bitIsoTkMu24 = 2 ;
//    else if (name.find("HLT_IsoTkMu24_v4")  != string::npos) bitIsoTkMu24 = 3 ;
//    
//    // Double Iso Mu
//    int bitMu17Mu8   = -1;
//    if      (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v1") != string::npos) bitMu17Mu8 = 0 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2") != string::npos) bitMu17Mu8 = 1 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v3") != string::npos) bitMu17Mu8 = 2 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4") != string::npos) bitMu17Mu8 = 3 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v5") != string::npos) bitMu17Mu8 = 4 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v6") != string::npos) bitMu17Mu8 = 5 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v7") != string::npos) bitMu17Mu8 = 6 ;
//    
//    int bitMu17TkMu8 = -1;
//    if      (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v1") != string::npos) bitMu17TkMu8 = 0 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2") != string::npos) bitMu17TkMu8 = 1 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3") != string::npos) bitMu17TkMu8 = 2 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v4") != string::npos) bitMu17TkMu8 = 3 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v5") != string::npos) bitMu17TkMu8 = 4 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v6") != string::npos) bitMu17TkMu8 = 5 ;                   
//
//    int bitMu17Mu8_noDZ   = -1;
//    if      (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1") != string::npos) bitMu17Mu8_noDZ = 0 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v2") != string::npos) bitMu17Mu8_noDZ = 1 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v3") != string::npos) bitMu17Mu8_noDZ = 2 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v4") != string::npos) bitMu17Mu8_noDZ = 3 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v5") != string::npos) bitMu17Mu8_noDZ = 4 ;
//    else if (name.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v6") != string::npos) bitMu17Mu8_noDZ = 5 ;
//    
//    int bitMu17TkMu8_noDZ = -1;
//    if      (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1") != string::npos) bitMu17TkMu8_noDZ = 0 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v2") != string::npos) bitMu17TkMu8_noDZ = 1 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v3") != string::npos) bitMu17TkMu8_noDZ = 2 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v4") != string::npos) bitMu17TkMu8_noDZ = 3 ;                   
//    else if (name.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v5") != string::npos) bitMu17TkMu8_noDZ = 4 ;                   
//
//    //Mu - EG
//    int bitMu8Ele23 = -1;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1") != string::npos) bitMu8Ele23 = 0 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2") != string::npos) bitMu8Ele23 = 1 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v3") != string::npos) bitMu8Ele23 = 2 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v4") != string::npos) bitMu8Ele23 = 3 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5") != string::npos) bitMu8Ele23 = 4 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v6") != string::npos) bitMu8Ele23 = 5 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v7") != string::npos) bitMu8Ele23 = 6 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v8") != string::npos) bitMu8Ele23 = 7 ;
//    if      (name.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v9") != string::npos) bitMu8Ele23 = 8 ; //check
//
//    int bitMu23Ele12 = -1;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1") != string::npos) bitMu23Ele12 = 0;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2") != string::npos) bitMu23Ele12 = 1;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3") != string::npos) bitMu23Ele12 = 2;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v4") != string::npos) bitMu23Ele12 = 3;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5") != string::npos) bitMu23Ele12 = 4;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v6") != string::npos) bitMu23Ele12 = 5;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v7") != string::npos) bitMu23Ele12 = 6;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v8") != string::npos) bitMu23Ele12 = 7;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v9") != string::npos) bitMu23Ele12 = 8;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v10") != string::npos) bitMu23Ele12 = 9;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v11") != string::npos) bitMu23Ele12 = 10;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v12") != string::npos) bitMu23Ele12 = 11;
//
//    int bitMu12Ele23_DZ = -1;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v1") != string::npos)  bitMu12Ele23_DZ = 0;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v2") != string::npos)  bitMu12Ele23_DZ = 1;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v3") != string::npos)  bitMu12Ele23_DZ = 2;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v4") != string::npos)  bitMu12Ele23_DZ = 3;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v5") != string::npos)  bitMu12Ele23_DZ = 4;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v6") != string::npos)  bitMu12Ele23_DZ = 5;
//    if     (name.find("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v7") != string::npos)  bitMu12Ele23_DZ = 6;
//
//    int bitMu23Ele12_DZ = -1;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v1") != string::npos)  bitMu23Ele12_DZ = 0;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2") != string::npos)  bitMu23Ele12_DZ = 1;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v3") != string::npos)  bitMu23Ele12_DZ = 2;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v4") != string::npos)  bitMu23Ele12_DZ = 3;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v5") != string::npos)  bitMu23Ele12_DZ = 4;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v6") != string::npos)  bitMu23Ele12_DZ = 5;
//    if     (name.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v7") != string::npos)  bitMu23Ele12_DZ = 6;
//
//    // Single Photon
//    int bitPhoton90 = -1;
//    if     (name.find("HLT_Photon90_v1") != string::npos)  bitPhoton90 = 0;
//    if     (name.find("HLT_Photon90_v2") != string::npos)  bitPhoton90 = 1;
//    if     (name.find("HLT_Photon90_v3") != string::npos)  bitPhoton90 = 2;
//    if     (name.find("HLT_Photon90_v4") != string::npos)  bitPhoton90 = 3;
//    if     (name.find("HLT_Photon90_v5") != string::npos)  bitPhoton90 = 4;
//    if     (name.find("HLT_Photon90_v6") != string::npos)  bitPhoton90 = 5;
//    if     (name.find("HLT_Photon90_v7") != string::npos)  bitPhoton90 = 6;
//    if     (name.find("HLT_Photon90_v8") != string::npos)  bitPhoton90 = 7;
//    if     (name.find("HLT_Photon90_v9") != string::npos)  bitPhoton90 = 8;
//    if     (name.find("HLT_Photon90_v10") != string::npos)  bitPhoton90 = 9;
//
//    int bitPhoton120 = -1;
//    if     (name.find("HLT_Photon120_v1") != string::npos)  bitPhoton120 = 0;
//    if     (name.find("HLT_Photon120_v2") != string::npos)  bitPhoton120 = 1;
//    if     (name.find("HLT_Photon120_v3") != string::npos)  bitPhoton120 = 2;
//    if     (name.find("HLT_Photon120_v4") != string::npos)  bitPhoton120 = 3;
//    if     (name.find("HLT_Photon120_v5") != string::npos)  bitPhoton120 = 4;
//    if     (name.find("HLT_Photon120_v6") != string::npos)  bitPhoton120 = 5;
//    if     (name.find("HLT_Photon120_v7") != string::npos)  bitPhoton120 = 6;
//    if     (name.find("HLT_Photon120_v8") != string::npos)  bitPhoton120 = 7;
//    if     (name.find("HLT_Photon120_v9") != string::npos)  bitPhoton120 = 8;
//    if     (name.find("HLT_Photon120_v10") != string::npos)  bitPhoton120 = 9;
//
//    int bitPhoton175 = -1;
//    if     (name.find("HLT_Photon175_v1") != string::npos)  bitPhoton175 = 0;
//    if     (name.find("HLT_Photon175_v2") != string::npos)  bitPhoton175 = 1;
//    if     (name.find("HLT_Photon175_v3") != string::npos)  bitPhoton175 = 2;
//    if     (name.find("HLT_Photon175_v4") != string::npos)  bitPhoton175 = 3;
//    if     (name.find("HLT_Photon175_v5") != string::npos)  bitPhoton175 = 4;
//    if     (name.find("HLT_Photon175_v6") != string::npos)  bitPhoton175 = 5;
//    if     (name.find("HLT_Photon175_v7") != string::npos)  bitPhoton175 = 6;
//    if     (name.find("HLT_Photon175_v8") != string::npos)  bitPhoton175 = 7;
//    if     (name.find("HLT_Photon175_v9") != string::npos)  bitPhoton175 = 8;
//    if     (name.find("HLT_Photon175_v10") != string::npos)  bitPhoton175 = 9;
//    if     (name.find("HLT_Photon175_v11") != string::npos)  bitPhoton175 = 10;
//
//    int bitPhoton165_HE10 = -1;
//    if     (name.find("LT_Photon165_HE10_v1") != string::npos)  bitPhoton165_HE10 = 0;
//    if     (name.find("LT_Photon165_HE10_v2") != string::npos)  bitPhoton165_HE10 = 1;
//    if     (name.find("LT_Photon165_HE10_v3") != string::npos)  bitPhoton165_HE10 = 2;
//    if     (name.find("LT_Photon165_HE10_v4") != string::npos)  bitPhoton165_HE10 = 3;
//    if     (name.find("LT_Photon165_HE10_v5") != string::npos)  bitPhoton165_HE10 = 4;
//    if     (name.find("LT_Photon165_HE10_v6") != string::npos)  bitPhoton165_HE10 = 5;
//    if     (name.find("LT_Photon165_HE10_v7") != string::npos)  bitPhoton165_HE10 = 6;
//    if     (name.find("LT_Photon165_HE10_v8") != string::npos)  bitPhoton165_HE10 = 7;
//    if     (name.find("LT_Photon165_HE10_v9") != string::npos)  bitPhoton165_HE10 = 8;
//    if     (name.find("LT_Photon165_HE10_v10") != string::npos)  bitPhoton165_HE10 = 9;
//    if     (name.find("LT_Photon165_HE10_v11") != string::npos)  bitPhoton165_HE10 = 10;
//    
//
//    //printf(" Reading trigger: %s\n" , name.c_str()  );
//    //printf(" bitEle23Loose      %d \n", bitEle23Loose      ); 
//    //printf(" bitEle27Tight      %d \n", bitEle27Tight      ); 
//    //printf(" bitEle17Ele12      %d \n", bitEle17Ele12      ); 
//    //printf(" bitEle23Ele12      %d \n", bitEle23Ele12      ); 
//    //printf(" bitIsoMu22         %d \n", bitIsoMu22         ); 
//    //printf(" bitIsoTkMu22       %d \n", bitIsoTkMu22       ); 
//    //printf(" bitMu17Mu8         %d \n", bitMu17Mu8         ); 
//    //printf(" bitMu17TkMu8       %d \n", bitMu17TkMu8       ); 
//    
//    // indicates prescaling and whether trigger was fired or not
//    unsigned int isPrescaled = (hltPrescale_.prescaleValue(e,es,name)) != 1;
//    const unsigned int triggerIndex = hltConfig_.triggerIndex(name);
//    unsigned int isFired     = AODTriggerHandle_->accept(triggerIndex);
//    //cout << "was run " << AODTriggerHandle_->wasrun(triggerIndex) << endl;
//    
//    // if statement checks if trigger was present (any version XXX_vN)
//    // set the bit N of our trigger variable to 1 or 0 based on isFired
//    if ( bitEle23Loose      >= 0 ){
//      AOD_HLT_Ele23Loose_            |= (isFired     << bitEle23Loose      );
//      AOD_HLT_Ele23Loose_isPS_       |= (isPrescaled << bitEle23Loose      );
//    }
//    if ( bitEle27Tight      >= 0 ){     
//      AOD_HLT_Ele27Tight_            |= (isFired     << bitEle27Tight      );
//      AOD_HLT_Ele27Tight_isPS_       |= (isPrescaled << bitEle27Tight      );
//    }
//    if ( bitEle17Ele12      >= 0 ){     
//      AOD_HLT_Ele17Ele12_            |= (isFired     << bitEle17Ele12      );
//      AOD_HLT_Ele17Ele12_isPS_       |= (isPrescaled << bitEle17Ele12      );
//    }
//    if ( bitEle23Ele12      >= 0 ){     
//      AOD_HLT_Ele23Ele12_            |= (isFired     << bitEle23Ele12      );
//      AOD_HLT_Ele23Ele12_isPS_       |= (isPrescaled << bitEle23Ele12      );
//    }
//    if ( bitIsoMu22         >= 0 ){   
//      AOD_HLT_IsoMu22_               |= (isFired     << bitIsoMu22         );
//      AOD_HLT_IsoMu22_isPS_          |= (isPrescaled << bitIsoMu22         );
//    }
//    if ( bitIsoTkMu22       >= 0 ){     
//      AOD_HLT_IsoTkMu22_             |= (isFired     << bitIsoTkMu22       );
//      AOD_HLT_IsoTkMu22_isPS_        |= (isPrescaled << bitIsoTkMu22       );
//    }
//    if ( bitIsoMu24         >= 0 ){   
//      AOD_HLT_IsoMu24_               |= (isFired     << bitIsoMu24         );
//      AOD_HLT_IsoMu24_isPS_          |= (isPrescaled << bitIsoMu24         );
//    }
//    if ( bitIsoTkMu24       >= 0 ){     
//      AOD_HLT_IsoTkMu24_             |= (isFired     << bitIsoTkMu24       );
//      AOD_HLT_IsoTkMu24_isPS_        |= (isPrescaled << bitIsoTkMu24       );
//    }
//    if ( bitMu17Mu8         >= 0 ){   
//      AOD_HLT_Mu17Mu8_               |= (isFired     << bitMu17Mu8         );
//      AOD_HLT_Mu17Mu8_isPS_          |= (isPrescaled << bitMu17Mu8         );
//    }
//    if ( bitMu17TkMu8       >= 0 ){     
//      AOD_HLT_Mu17TkMu8_             |= (isFired     << bitMu17TkMu8       );
//      AOD_HLT_Mu17TkMu8_isPS_        |= (isPrescaled << bitMu17TkMu8       );
//    }
//    if ( bitMu17Mu8_noDZ         >= 0 ){   
//      AOD_HLT_Mu17Mu8_noDZ_          |= (isFired     << bitMu17Mu8_noDZ    );
//      AOD_HLT_Mu17Mu8_noDZ_isPS_     |= (isPrescaled << bitMu17Mu8_noDZ    );
//    }
//    if ( bitMu17TkMu8_noDZ       >= 0 ){     
//      AOD_HLT_Mu17TkMu8_noDZ_        |= (isFired     << bitMu17TkMu8_noDZ  );
//      AOD_HLT_Mu17TkMu8_noDZ_isPS_   |= (isPrescaled << bitMu17TkMu8_noDZ  );
//    }
//    if ( bitMu8Ele23 >= 0 ){
//      AOD_HLT_Mu8Ele23_       |= (isFired     << bitMu8Ele23 );
//      AOD_HLT_Mu8Ele23_isPS_  |= (isPrescaled << bitMu8Ele23 );
//    }
//    if ( bitMu23Ele12 >= 0 ){
//      AOD_HLT_Mu23Ele12_       |= (isFired     << bitMu23Ele12 );
//      AOD_HLT_Mu23Ele12_isPS_  |= (isPrescaled << bitMu23Ele12 );
//    }
//    if ( bitMu12Ele23_DZ >= 0 ){
//      AOD_HLT_Mu12Ele23_DZ_       |= (isFired     << bitMu12Ele23_DZ );
//      AOD_HLT_Mu12Ele23_DZ_isPS_  |= (isPrescaled << bitMu12Ele23_DZ );
//    }
//    if ( bitMu23Ele12_DZ >= 0 ){
//      AOD_HLT_Mu23Ele12_DZ_       |= (isFired     << bitMu23Ele12_DZ );
//      AOD_HLT_Mu23Ele12_DZ_isPS_  |= (isPrescaled << bitMu23Ele12_DZ );
//    }
//    if ( bitPhoton90 >= 0 ){
//      AOD_HLT_Photon90_       |= (isFired     << bitPhoton90 );
//      AOD_HLT_Photon90_isPS_  |= (isPrescaled << bitPhoton90 );
//    }
//    if ( bitPhoton120 >= 0 ){
//      AOD_HLT_Photon120_       |= (isFired     << bitPhoton120 );
//      AOD_HLT_Photon120_isPS_  |= (isPrescaled << bitPhoton120 );
//    }
//    if ( bitPhoton175 >= 0 ){
//      AOD_HLT_Photon175_        |= (isFired     << bitPhoton175 );
//      AOD_HLT_Photon175_isPS_   |= (isPrescaled << bitPhoton175 );
//    }
//    if ( bitPhoton165_HE10 >= 0 ){
//      AOD_HLT_Photon165_HE10_       |= (isFired     << bitPhoton165_HE10 );
//      AOD_HLT_Photon165_HE10_isPS_  |= (isPrescaled << bitPhoton165_HE10 );
//    }
//    
//
//    
//  }//hltConfig loop
  
}//fillMiniAODTrigger
