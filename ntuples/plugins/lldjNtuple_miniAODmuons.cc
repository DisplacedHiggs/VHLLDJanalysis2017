#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"

using namespace std;

// (local) variables associated with tree branchesMiniAOD
Int_t            miniAOD_nMu_                              ; 
vector<float>    miniAOD_muPt_                            ; 
vector<float>    miniAOD_muEn_                            ; 
vector<float>    miniAOD_muEta_                           ; 
vector<float>    miniAOD_muPhi_                           ; 
vector<int>      miniAOD_muCharge_                        ; 
vector<int>      miniAOD_muType_                          ; 
vector<bool>     miniAOD_muIsGlobalMuon_                  ; 
vector<bool>     miniAOD_muIsPFMuon_                      ; 
vector<bool>     miniAOD_muPassLooseID_                   ;
vector<bool>     miniAOD_muPassMediumBCDEFID_             ;
vector<bool>     miniAOD_muPassMediumGHID_                ;
vector<bool>     miniAOD_muPassTightID_                   ;
vector<float>    miniAOD_muPFdBetaIsolation_              ; 
vector<float>    miniAOD_muDxy_                           ;
vector<float>    miniAOD_muDxyErr_                        ;
vector<float>    miniAOD_muDB_BS2D_                       ;
vector<float>    miniAOD_muDB_PV2D_                       ;

////https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_Mori
//bool lldjNtuple::isMediumMuonBCDEF(const reco::Muon & recoMu) 
//{
//  bool goodGlob = recoMu.isGlobalMuon() && 
//    recoMu.globalTrack()->normalizedChi2() < 3 && 
//    recoMu.combinedQuality().chi2LocalPosition < 12 && 
//    recoMu.combinedQuality().trkKink < 20; 
//  bool isMedium = muon::isLooseMuon(recoMu) && 
//    recoMu.innerTrack()->validFraction() > 0.49 && 
//    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
//  return isMedium; 
//}
//
//
////https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Short_Term_Instructions_for_Mori
//bool lldjNtuple::isMediumMuonGH(const reco::Muon & recoMu) 
//{
//  bool goodGlob = recoMu.isGlobalMuon() && 
//    recoMu.globalTrack()->normalizedChi2() < 3 && 
//    recoMu.combinedQuality().chi2LocalPosition < 12 && 
//    recoMu.combinedQuality().trkKink < 20; 
//  bool isMedium = muon::isLooseMuon(recoMu) && 
//    recoMu.innerTrack()->validFraction() > 0.8 && 
//    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
//  return isMedium; 
//}


void lldjNtuple::branchesMiniAODMuons(TTree* tree) {
 tree->Branch("miniAOD_nMu",                             &miniAOD_nMu_                             ) ; 
 tree->Branch("miniAOD_muPt",                           &miniAOD_muPt_                           ) ; 
 tree->Branch("miniAOD_muEn",                           &miniAOD_muEn_                           ) ; 
 tree->Branch("miniAOD_muEta",                          &miniAOD_muEta_                          ) ; 
 tree->Branch("miniAOD_muPhi",                          &miniAOD_muPhi_                          ) ; 
 tree->Branch("miniAOD_muCharge",                       &miniAOD_muCharge_                       ) ; 
 tree->Branch("miniAOD_muType",                         &miniAOD_muType_                         ) ; 
 tree->Branch("miniAOD_muIsGlobalMuon",                 &miniAOD_muIsGlobalMuon_                 ) ; 
 tree->Branch("miniAOD_muIsPFMuon",                     &miniAOD_muIsPFMuon_                     ) ; 
 tree->Branch("miniAOD_muPassLooseID",                  &miniAOD_muPassLooseID_                  ) ; 
 tree->Branch("miniAOD_muPassMediumBCDEFID",            &miniAOD_muPassMediumBCDEFID_            ) ; 
 tree->Branch("miniAOD_muPassMediumGHID",               &miniAOD_muPassMediumGHID_               ) ; 
 tree->Branch("miniAOD_muPassTightID",                  &miniAOD_muPassTightID_                  ) ; 
 tree->Branch("miniAOD_muPFdBetaIsolation",             &miniAOD_muPFdBetaIsolation_             ) ; 
 tree->Branch("miniAOD_muDxy",                          &miniAOD_muDxy_                          ) ; 
 tree->Branch("miniAOD_muDxyErr",                       &miniAOD_muDxyErr_                       ) ; 
 tree->Branch("miniAOD_muDB_BS2D",                      &miniAOD_muDB_BS2D_                       ) ; 
 tree->Branch("miniAOD_muDB_PV2D",                      &miniAOD_muDB_PV2D_                       ) ; 
}


//void lldjNtuple::fillMiniAODMuons(const edm::Event& e, reco::Vertex vtx) {
void lldjNtuple::fillMiniAODMuons(const edm::Event& e, const edm::EventSetup& es) {

 // cleanup from previous execution
 miniAOD_nMu_ = 0;
 miniAOD_muPt_                          .clear() ; 
 miniAOD_muEn_                          .clear() ; 
 miniAOD_muEta_                         .clear() ; 
 miniAOD_muPhi_                         .clear() ; 
 miniAOD_muCharge_                      .clear() ; 
 miniAOD_muType_                        .clear() ; 
 miniAOD_muIsGlobalMuon_                .clear() ; 
 miniAOD_muIsPFMuon_                    .clear() ; 
 miniAOD_muPassLooseID_                 .clear() ;
 miniAOD_muPassMediumBCDEFID_           .clear() ;
 miniAOD_muPassMediumGHID_              .clear() ;
 miniAOD_muPassTightID_                 .clear() ;
 miniAOD_muPFdBetaIsolation_            .clear() ; 
 miniAOD_muDxy_                         .clear() ; 
 miniAOD_muDxyErr_                      .clear() ; 
 miniAOD_muDB_BS2D_                     .clear() ; 
 miniAOD_muDB_PV2D_                     .clear() ; 

 edm::Handle<edm::View<pat::Muon> > miniAODmuonHandle;
 e.getByToken(miniAODmuonCollection_, miniAODmuonHandle);
 
 if (!miniAODmuonHandle.isValid()) {
   edm::LogWarning("lldjNtuple") << " missing miniAOD muon collection";
  return;
 }

// //Beamspot for impact parameter
// edm::Handle<reco::BeamSpot> beamspotHandle_;
// e.getByToken(beamspotLabel_, beamspotHandle_);

 for (edm::View<pat::Muon>::const_iterator iMu = miniAODmuonHandle->begin(); iMu != miniAODmuonHandle->end(); ++iMu) {

  Float_t pt = iMu->pt();
  Float_t eta = iMu->eta();

  if (pt < 2) continue;
  if (fabs(eta) > 2.4) continue;
  if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;

  miniAOD_nMu_++;

  //const reco::Muon &recoMu = dynamic_cast<const reco::Muon &>(*iMu);

  //Basic info
  miniAOD_muPt_           .push_back(pt);
  miniAOD_muEn_           .push_back(iMu->energy());
  miniAOD_muEta_          .push_back(iMu->eta());
  miniAOD_muPhi_          .push_back(iMu->phi());
  miniAOD_muCharge_       .push_back(iMu->charge());
  miniAOD_muType_         .push_back(iMu->type());
  miniAOD_muIsGlobalMuon_ .push_back(iMu->isGlobalMuon () ) ; 
  miniAOD_muIsPFMuon_     .push_back(iMu->isPFMuon     () ) ; 

  //ID
  miniAOD_muPassLooseID_       .push_back( iMu->isLooseMuon() );
  //miniAOD_muPassMediumBCDEFID_ .push_back( isMediumMuonBCDEF(recoMu) );
  //miniAOD_muPassMediumGHID_    .push_back( isMediumMuonGH(recoMu) );
  //miniAOD_muPassTightID_       .push_back( iMu->isTightMuon(vtx) );

  //Isolation
  Float_t muPFChIso      = iMu->pfIsolationR04().sumChargedHadronPt ;
  Float_t muPFPhoIso     = iMu->pfIsolationR04().sumPhotonEt        ;
  Float_t muPFNeuIso     = iMu->pfIsolationR04().sumNeutralHadronEt ;
  Float_t muPFPUIso      = iMu->pfIsolationR04().sumPUPt            ;
  Float_t pfdBetaIso     = ( muPFChIso + max(0.0,muPFNeuIso + muPFPhoIso - 0.5*muPFPUIso ) ) / pt ;
  miniAOD_muPFdBetaIsolation_.push_back( pfdBetaIso     ) ;   

  ////Impact parameter
  //miniAOD_muDB_BS2D_.push_back (iMu->dB(pat::Muon::BS2D));
  //miniAOD_muDB_PV2D_.push_back (iMu->dB(pat::Muon::PV2D));
  //if(iMu->innerTrack().isNonnull()){
  //  float dxy = fabs(iMu->innerTrack()->dxy(*beamspotHandle_));
  //  float dxyErr = fabs(iMu->innerTrack()->dxyError());
  //  miniAOD_muDxy_.push_back( dxy );
  //  miniAOD_muDxyErr_.push_back( dxyErr );
  //}
  //else{
  // miniAOD_muDxy_.push_back( -1 ) ;
  // miniAOD_muDxyErr_.push_back( -1 ) ;
  //}
 }//End muon collection loop

}//End fillMiniAODMuons
