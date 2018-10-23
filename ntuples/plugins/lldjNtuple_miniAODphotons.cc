#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"

using namespace std;

Int_t            miniAOD_nPho_;
vector<float>    miniAOD_phoPt_;
vector<float>    miniAOD_phoEn_;
vector<float>    miniAOD_phoEta_;
vector<float>    miniAOD_phoPhi_;

vector<float>    miniAOD_phoSCEn_;
vector<float>    miniAOD_phoSCEta_;
vector<float>    miniAOD_phoSCPhi_;

vector<float>    miniAOD_phoPassElectronVeto_;
vector<float>    miniAOD_phoHasPixelSeed_;

vector<UShort_t> miniAOD_phoIDbit_;
//vector<Float_t>  miniAOD_phoIDMVA_;

vector<float>    miniAOD_phoObjPFChIso_;
vector<float>    miniAOD_phoObjPFPhoIso_;
vector<float>    miniAOD_phoObjPFNeuIso_;
vector<float>    miniAOD_phoObjPFChWorstIso_;

vector<float>    miniAOD_phoMapPFChIso_;
vector<float>    miniAOD_phoMapPFPhoIso_;
vector<float>    miniAOD_phoMapPFNeuIso_;
vector<float>    miniAOD_phoMapPFChWorstIso_;


void lldjNtuple::branchesMiniAODPhotons(TTree* tree) {

 tree->Branch("miniAOD_nPho",                &miniAOD_nPho_                  );             
 tree->Branch("miniAOD_phoPt",              &miniAOD_phoPt_                );             
 tree->Branch("miniAOD_phoEn",              &miniAOD_phoEn_                );             
 tree->Branch("miniAOD_phoEta",             &miniAOD_phoEta_               );             
 tree->Branch("miniAOD_phoPhi",             &miniAOD_phoPhi_               );             
                                                   
 tree->Branch("miniAOD_phoSCEn",            &miniAOD_phoSCEn_              );             
 tree->Branch("miniAOD_phoSCEta",           &miniAOD_phoSCEta_             );             
 tree->Branch("miniAOD_phoSCPhi",           &miniAOD_phoSCPhi_             );             

 tree->Branch("miniAOD_phoSCPhi",           &miniAOD_phoSCPhi_             );             

 tree->Branch("miniAOD_phoPassElectronVeto", &miniAOD_phoPassElectronVeto_ );
 tree->Branch("miniAOD_phoHasPixelSeed",     &miniAOD_phoHasPixelSeed_     );

 tree->Branch("miniAOD_phoIDbit",           &miniAOD_phoIDbit_           );             
 //tree->Branch("miniAOD_phoIDMVA",           &miniAOD_phoIDMVA_           );             
  
 tree->Branch("miniAOD_phoObjPFChIso",      &miniAOD_phoObjPFChIso_      );
 tree->Branch("miniAOD_phoObjPFPhoIso",     &miniAOD_phoObjPFPhoIso_     );
 tree->Branch("miniAOD_phoObjPFNeuIso",     &miniAOD_phoObjPFNeuIso_     );
 tree->Branch("miniAOD_phoObjPFChWorstIso", &miniAOD_phoObjPFChWorstIso_ );
  
 tree->Branch("miniAOD_phoMapPFChIso",      &miniAOD_phoMapPFChIso_      );
 tree->Branch("miniAOD_phoMapPFPhoIso",     &miniAOD_phoMapPFPhoIso_     );
 tree->Branch("miniAOD_phoMapPFNeuIso",     &miniAOD_phoMapPFNeuIso_     );
 tree->Branch("miniAOD_phoMapPFChWorstIso", &miniAOD_phoMapPFChWorstIso_ );

}

void lldjNtuple::fillMiniAODPhotons(const edm::Event& e, const edm::EventSetup& es) {
 
 // cleanup from previous execution
 miniAOD_nPho_                 = 0;
 miniAOD_phoPt_              .clear();    
 miniAOD_phoEn_              .clear();    
 miniAOD_phoEta_             .clear();     
 miniAOD_phoPhi_             .clear();     

 miniAOD_phoSCEn_            .clear();      
 miniAOD_phoSCEta_           .clear();       
 miniAOD_phoSCPhi_           .clear();       

 miniAOD_phoPassElectronVeto_ .clear();
 miniAOD_phoHasPixelSeed_     .clear();

 miniAOD_phoIDbit_           .clear();       
 //miniAOD_phoIDMVA_           .clear();       

 miniAOD_phoObjPFChIso_      .clear();
 miniAOD_phoObjPFPhoIso_     .clear();
 miniAOD_phoObjPFNeuIso_     .clear();
 miniAOD_phoObjPFChWorstIso_ .clear();

 miniAOD_phoMapPFChIso_      .clear();
 miniAOD_phoMapPFPhoIso_     .clear();
 miniAOD_phoMapPFNeuIso_     .clear();
 miniAOD_phoMapPFChWorstIso_ .clear();

 //edm::Handle<edm::View<reco::Photon> > miniAODphotonHandle;
 edm::Handle<edm::View<pat::Photon> > miniAODphotonHandle;
 e.getByToken(miniAODphotonCollection_, miniAODphotonHandle);

 if (!miniAODphotonHandle.isValid()) {
   std::cout << "Invalid photon handle" << std::endl;
   edm::LogWarning("lldjNtuple") << "no pat::Photons in event";
   return;
 } 
 
// if(!(e.getByToken(miniAOD_phoLooseIdMapToken_ ,loose_id_decisions) 
//    && e.getByToken(miniAOD_phoMediumIdMapToken_ ,medium_id_decisions) 
//    && e.getByToken(miniAOD_phoTightIdMapToken_ ,tight_id_decisions) 
//    && e.getByToken(miniAOD_phoChargedIsolationMapToken_ ,miniAOD_phoChargedIsolationHandle_) 
//    && e.getByToken(miniAOD_phoNeutralHadronIsolationMapToken_ ,miniAOD_phoNeutralHadronIsolationHandle_) 
//    && e.getByToken(miniAOD_phoPhotonIsolationMapToken_ ,miniAOD_phoPhotonIsolationHandle_) 
//    && e.getByToken(miniAOD_phoWorstChargedIsolationMapToken_ ,miniAOD_phoWorstChargedIsolationHandle_)
//      )){
//   edm::LogWarning("lldjNtuple") << "Failure in ID tokens";
//   return;
// } 

 //https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
 //https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_8.0.3/PhotonNtupler/test/runPhotons_VID_CutBased_Spring16_demo.py
 for (size_t i = 0; i < miniAODphotonHandle->size(); ++i){
   const auto pho = miniAODphotonHandle->ptrAt(i);

   Float_t miniAOD_phoPt  = pho->et();
   Float_t miniAOD_phoEta = pho->eta();
   
   // VID decisions    
   UShort_t tmpphoIDbit = 0;    
   bool isPassLoose  = pho->photonID("cutBasedPhotonID-Fall17-94X-V1-loose");
   if (isPassLoose)  setbit(tmpphoIDbit, 0);   
   bool isPassMedium = pho->photonID("cutBasedPhotonID-Fall17-94X-V1-medium");
   if (isPassMedium) setbit(tmpphoIDbit, 1);    
   bool isPassTight  = pho->photonID("cutBasedPhotonID-Fall17-94X-V1-tight");
   if (isPassTight)  setbit(tmpphoIDbit, 2); 
   
   miniAOD_phoIDbit_.push_back(tmpphoIDbit); 

   if (miniAOD_phoPt < 20 || fabs(miniAOD_phoEta) > 2.1) continue;
   
   miniAOD_phoPt_     .push_back(miniAOD_phoPt);
   miniAOD_phoEn_     .push_back(pho->energy());
   miniAOD_phoEta_    .push_back(miniAOD_phoEta);
   miniAOD_phoPhi_    .push_back(pho->phi());
   
   miniAOD_phoSCEn_   .push_back( (*pho).superCluster()->energy());      
   miniAOD_phoSCEta_  .push_back( (*pho).superCluster()->eta());       
   miniAOD_phoSCPhi_  .push_back( (*pho).superCluster()->phi());       

   //https://cmssdt.cern.ch/lxr/source/DataFormats/EgammaCandidates/interface/Photon.h?%21v=CMSSW_8_0_28
   miniAOD_phoPassElectronVeto_ .push_back( -1 ) ; //pho->passElectronVeto() ); 
   miniAOD_phoHasPixelSeed_     .push_back( pho->hasPixelSeed()     ); 
   
   miniAOD_phoObjPFChIso_       .push_back(pho->chargedHadronIso());
   miniAOD_phoObjPFPhoIso_      .push_back(pho->photonIso());
   miniAOD_phoObjPFNeuIso_      .push_back(pho->neutralHadronIso());
   
//   miniAOD_phoMapPFChIso_     .push_back((*miniAOD_phoChargedIsolationHandle_)[pho]);
//   miniAOD_phoMapPFPhoIso_    .push_back((*miniAOD_phoPhotonIsolationHandle_)[pho]);
//   miniAOD_phoMapPFNeuIso_    .push_back((*miniAOD_phoNeutralHadronIsolationHandle_)[pho]);
//   miniAOD_phoMapPFChWorstIso_.push_back((*miniAOD_phoWorstChargedIsolationHandle_)[pho]);
   
   miniAOD_nPho_++;

 }
}

