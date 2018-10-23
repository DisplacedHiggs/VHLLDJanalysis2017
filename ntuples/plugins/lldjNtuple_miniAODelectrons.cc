#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"


using namespace std;
//using namespace reco;


// (local) variables associated with tree branchesMiniAOD
Int_t            miniAOD_nEle_;
vector<float>    miniAOD_elePt_;
vector<float>    miniAOD_eleEn_;
vector<float>    miniAOD_eleEta_;
vector<float>    miniAOD_elePhi_;

vector<float>    miniAOD_eled0_;
vector<float>    miniAOD_eledz_;

vector<int>      miniAOD_eleCharge_;
vector<int>      miniAOD_eleChargeConsistent_;

vector<UShort_t> miniAOD_eleIDbit_;
vector<int>      miniAOD_elePassConversionVeto_;

void lldjNtuple::branchesMiniAODElectrons(TTree* tree) {

 tree->Branch("miniAOD_nEle",                    &miniAOD_nEle_                    );                           
 tree->Branch("miniAOD_elePt",                   &miniAOD_elePt_                   );     
 tree->Branch("miniAOD_eleEn",                   &miniAOD_eleEn_                   );     
 tree->Branch("miniAOD_eleEta",                  &miniAOD_eleEta_                  );     
 tree->Branch("miniAOD_elePhi",                  &miniAOD_elePhi_                  );     

 tree->Branch("miniAOD_eled0",                   &miniAOD_eled0_                   );     
 tree->Branch("miniAOD_eledz",                   &miniAOD_eledz_                   );     

 tree->Branch("miniAOD_eleCharge",               &miniAOD_eleCharge_               );     
 tree->Branch("miniAOD_eleChargeConsistent",     &miniAOD_eleChargeConsistent_     );     

 tree->Branch("miniAOD_eleIDbit",                &miniAOD_eleIDbit_                );     
 tree->Branch("miniAOD_elePassConversionVeto",   &miniAOD_elePassConversionVeto_   );

}

//void lldjNtuple::fillMiniAODElectrons(const edm::Event &e, const edm::EventSetup &es, reco::Vertex vtx) {
void lldjNtuple::fillMiniAODElectrons(const edm::Event &e, const edm::EventSetup &es) {
    
 // cleanup from previous execution
 miniAOD_nEle_                     = 0;

 miniAOD_elePt_                    .clear();     
 miniAOD_eleEn_                    .clear();     
 miniAOD_eleEta_                   .clear();     
 miniAOD_elePhi_                   .clear();     

 miniAOD_eled0_                   .clear();     
 miniAOD_eledz_                   .clear();     

 miniAOD_eleCharge_                .clear();     
 miniAOD_eleChargeConsistent_      .clear();     

 miniAOD_eleIDbit_                 .clear();     
 miniAOD_elePassConversionVeto_    .clear();


 //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
 //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/plugins/ElectronNtuplerVIDDemo.cc
 //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/ElectronNtupler/test/runElectrons_VID_CutBased_Summer16_HLTSafe_demo.py
 
  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(miniAODelectronCollection_, electronHandle);

// edm::Handle<edm::View<reco::GsfElectron> > electronHandle;
// e.getByToken(electronAODToken_, electronHandle);
 if( !electronHandle.isValid() ){
   std::cout << "Invalid electron handle" << std::endl;
   edm::LogWarning("lldjNtuple") << "no pat::Electrons in event";
   return;
 }

 //Skip PV for now (for dz)

// // Get the conversions collection
// edm::Handle<reco::ConversionCollection> conversions;
// e.getByToken(conversionsAODToken_, conversions);
//
// // Beamspot needed for conversion veto
// edm::Handle<reco::BeamSpot> theBeamSpot;
// e.getByToken(beamspotLabel_, theBeamSpot);

// edm::Handle<edm::ValueMap<bool> > ele_id_decisions_loose;
// e.getByToken(AOD_eleLooseIdMapToken_ ,ele_id_decisions_loose);
// edm::Handle<edm::ValueMap<bool> > ele_id_decisions_medium;
// e.getByToken(AOD_eleMediumIdMapToken_ ,ele_id_decisions_medium);
// edm::Handle<edm::ValueMap<bool> > ele_id_decisions_tight;
// e.getByToken(AOD_eleTightIdMapToken_ ,ele_id_decisions_tight);
 
 for (size_t i = 0; i < electronHandle->size(); ++i){
   const auto el = electronHandle->ptrAt(i);

   if( el->pt() < 5 ) continue;
   if( fabs(el->eta()) > 2.5 ) continue;
   
   miniAOD_eleEn_  .push_back( el->energy() );
   miniAOD_elePt_  .push_back( el->pt() );
   miniAOD_eleEta_ .push_back( el->superCluster()->eta() );
   miniAOD_elePhi_ .push_back( el->superCluster()->phi() );

   //reco::GsfTrackRef theTrack = el->gsfTrack();
   ////AOD_eled0_.push_back( theTrack->d0() );
   //miniAOD_eled0_.push_back( (-1) * theTrack->dxy( vtx.position() ) );
   //miniAOD_eledz_.push_back( theTrack->dz( vtx.position() ) );

   //std::cout << "d0 " << theTrack->d0() << ", " << dz << std::endl; 

   miniAOD_eleCharge_          .push_back(el->charge());
   miniAOD_eleChargeConsistent_.push_back((Int_t)el->isGsfCtfScPixChargeConsistent());
   
   //// Conversion rejection
   //bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, 
   //							      conversions,
   //							      theBeamSpot->position());
   //AOD_elePassConversionVeto_.push_back( (int) passConvVeto );


   //ID
   UShort_t tmpeleIDbit = 0;
   bool isPassVeto   = el->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
   if (isPassVeto)   setbit(tmpeleIDbit, 0);
   bool isPassLoose  = el->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
   if (isPassLoose)  setbit(tmpeleIDbit, 1);
   bool isPassMedium = el->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
   if (isPassMedium) setbit(tmpeleIDbit, 2);
   bool isPassTight  = el->electronID("cutBasedElectronID-Fall17-94X-V1-tight");
   if (isPassTight)  setbit(tmpeleIDbit, 3);
   bool isPassHEEP   = el->electronID("heepElectronID-HEEPV70");
   if (isPassHEEP)   setbit(tmpeleIDbit, 4);
   
   miniAOD_eleIDbit_.push_back(tmpeleIDbit);


   //UShort_t tmpeleIDbit = 0;
   //bool isPassEleLooseId  = (*ele_id_decisions_loose)[el];
   //if(isPassEleLooseId) setbit(tmpeleIDbit, 0);
   //bool isPassEleMediumId  = (*ele_id_decisions_medium)[el];
   //if(isPassEleMediumId) setbit(tmpeleIDbit, 1);
   //bool isPassEleTightId  = (*ele_id_decisions_tight)[el];
   //if(isPassEleTightId) setbit(tmpeleIDbit, 2);
   //AOD_eleIDbit_.push_back(tmpeleIDbit);

   
   miniAOD_nEle_++;
 }   
}
