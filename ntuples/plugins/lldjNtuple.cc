#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"

using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

lldjNtuple::lldjNtuple(const edm::ParameterSet& ps) :
  hltPrescale_(ps,consumesCollector(),*this) 
{
  
  //lldj_pset_ = ps;

  // LLDJ custom miniAOD
  doLLDJAOD_                     = ps.getParameter<bool>("doLLDJAOD");

  //trgFilterDeltaPtCut_       = ps.getParameter<double>("trgFilterDeltaPtCut");
  //trgFilterDeltaRCut_        = ps.getParameter<double>("trgFilterDeltaRCut");

  // Trigger
  miniAODtriggerObjectsLabel_       = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("miniAODtriggerEvent"));
  miniAODtrgResultsLabel_           = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("miniAODtriggerResults"));
  miniAODpatTrgResultsLabel_        = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("miniAODpatTriggerResults"));
  miniAODtrgResultsProcess_         =                                          ps.getParameter<InputTag>("miniAODtriggerResults").process();

  // Global Event
  miniAODrhoLabel_                  = consumes<double>                        (ps.getParameter<InputTag>("miniAODrhoLabel"));
  miniAODrhoCentralLabel_           = consumes<double>                        (ps.getParameter<InputTag>("miniAODrhoCentralLabel"));
  miniAODpuCollection_              = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("miniAODpileupCollection"));
  miniAODvtxLabel_                  = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("miniAODVtxLabel"));

  // Electrons
  miniAODelectronCollection_        = consumes<View<pat::Electron> >          (ps.getParameter<InputTag>("miniAODelectronSrc"));

  // Photons
  miniAODphotonCollection_          = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("miniAODphotonSrc"));

  // Muons
  miniAODmuonCollection_            = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("miniAODmuonSrc"));

  // MET
  //miniAODpatTrgResultsLabel_        = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("miniAODpatTriggerResults"));
  miniAODpfMETlabel_                = consumes<View<pat::MET> >               (ps.getParameter<InputTag>("miniAODpfMETLabel"));

  //miniAODvtxBSLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("miniAODVtxBSLabel"));
  miniAODtrgEventLabel_             = consumes<trigger::TriggerEvent>         (ps.getParameter<InputTag>("miniAODtriggerEvent"));
  miniAODgeneratorLabel_            = consumes<GenEventInfoProduct>           (ps.getParameter<InputTag>("miniAODgeneratorLabel"));
  miniAODlheEventLabel_             = consumes<LHEEventProduct>               (ps.getParameter<InputTag>("miniAODLHEEventLabel"));
  miniAODgenParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("miniAODgenParticleSrc"));
  miniAODgsfTracks_                 = consumes<View<reco::GsfTrack>>          (ps.getParameter<InputTag>("miniAODgsfTrackSrc"));

  miniAODBadChCandFilterToken_      = consumes<bool>                          (ps.getParameter<InputTag>("miniAODBadChargedCandidateFilter"));
  miniAODBadPFMuonFilterToken_      = consumes<bool>                          (ps.getParameter<edm::InputTag>("miniAODBadPFMuonFilter"));

  miniAODebReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("miniAODebReducedRecHitCollection"));
  miniAODeeReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("miniAODeeReducedRecHitCollection"));
  miniAODesReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("miniAODesReducedRecHitCollection")); 
  miniAODrecophotonCollection_      = consumes<reco::PhotonCollection>        (ps.getParameter<InputTag>("miniAODrecoPhotonSrc"));
  miniAODtracklabel_                = consumes<reco::TrackCollection>         (ps.getParameter<InputTag>("miniAODTrackLabel"));
  miniAODgsfElectronlabel_          = consumes<reco::GsfElectronCollection>   (ps.getParameter<InputTag>("miniAODgsfElectronLabel"));
  miniAODtauCollection_             = consumes<vector<pat::Tau> >             (ps.getParameter<InputTag>("miniAODtauSrc"));
  miniAODpfAllParticles_            = consumes<reco::PFCandidateCollection>   (ps.getParameter<InputTag>("miniAODPFAllCandidates"));
  miniAODpckPFCandidateCollection_  = consumes<pat::PackedCandidateCollection>(ps.getParameter<InputTag>("miniAODpackedPFCands"));
  miniAODpckPFCdsLabel_             = consumes<vector<pat::PackedCandidate>>  (ps.getParameter<InputTag>("miniAODpackedPFCands"));
  miniAODrecoCdsLabel_              = consumes<View<reco::Candidate>>         (ps.getParameter<InputTag>("miniAODpackedPFCands"));

  miniAODjetsAK4Label_              = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("miniAODak4JetSrc"));
  miniAODjetsAK8Label_              = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("miniAODak8JetSrc"));


//  // choose AOD or miniAOD or both
//  doAOD_                   = ps.getParameter<bool>("doAOD");
//  doMiniAOD_               = ps.getParameter<bool>("doMiniAOD");
//
//  // electrons
//  electronCollection_      = consumes<View<pat::Electron> > (ps.getParameter<InputTag>("electronSrc"));
//  rhoLabel_                = consumes<double>               (ps.getParameter<InputTag>("rhoLabel"));
//  // electron ID 
//  eleVetoIdMapToken_       = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleVetoIdMap"));
//  eleLooseIdMapToken_      = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleLooseIdMap"));
//  eleMediumIdMapToken_     = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleMediumIdMap"));
//  eleTightIdMapToken_      = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleTightIdMap"));
//  eleHLTIdMapToken_        = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleHLTIdMap"));
//  //eleHEEPIdMapToken_       = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleHEEPIdMap"));
//  //eleMVAValuesMapToken_    = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleMVAValuesMap"));
//  //eleMVAHZZValuesMapToken_ = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleMVAHZZValuesMap"));
//  //elePFClusEcalIsoToken_   = mayConsume<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("elePFClusEcalIsoProducer"));
//  //elePFClusHcalIsoToken_   = mayConsume<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("elePFClusHcalIsoProducer"));
//
//  // AOD electron
//  electronAODToken_    = mayConsume<edm::View<reco::GsfElectron> >(ps.getParameter<edm::InputTag>("electronAODSrc"));
//  conversionsAODToken_ = mayConsume< reco::ConversionCollection >(ps.getParameter<edm::InputTag>("conversions"));
//  AOD_eleLooseIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("AOD_eleLooseIdMap"));
//  AOD_eleMediumIdMapToken_   = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("AOD_eleMediumIdMap"));
//  AOD_eleTightIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("AOD_eleTightIdMap"));
//
//
//  // global event
//  rhoCentralLabel_         = consumes<double>                        (ps.getParameter<InputTag>("rhoCentralLabel"));
//  puCollection_            = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
//  AODpuCollection_         = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("AODpileupCollection"));
//  vtxLabel_                = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));
//  //AODVertexLabel_          = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("AODVertexSrc"));
//  AODVertexLabel_          = consumes<edm::View<reco::Vertex> >      (ps.getParameter<InputTag>("AODVertexSrc"));
//  trgResultsLabel_         = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("triggerResults"));
//  trgResultsProcess_       =                                          ps.getParameter<InputTag>("triggerResults").process();
//
//  // beamspot 
//  beamspotLabel_           = consumes<reco::BeamSpot>                (ps.getParameter<InputTag>("beamspotLabel_"));
//
//  // jets
//  jetsAK4Label_            = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak4JetSrc"));
//  //AODjetsAK4Label_         = consumes<View<pat::Jet> >               (ps.getParameter<InputTag>("ak4JetSrc"));
//  AODak4CaloJetsLabel_     = consumes<View<reco::CaloJet> >          (ps.getParameter<InputTag>("AODak4CaloJetsSrc"));  
//  AODak4PFJetsLabel_       = consumes<View<reco::PFJet>   >          (ps.getParameter<InputTag>("AODak4PFJetsSrc"));    
//  AODak4PFJetsCHSLabel_    = consumes<View<reco::PFJet>   >          (ps.getParameter<InputTag>("AODak4PFJetsCHSSrc")); 
//  selectedPatJetsLabel_    = consumes<edm::View<pat::Jet> >          (ps.getParameter<InputTag>("selectedPatJetsSrc"));
//  AODTrackLabel_           = consumes<edm::View<reco::Track> >       (ps.getParameter<InputTag>("AODTrackSrc"));
//
//  // met
//  patTrgResultsLabel_      = consumes<edm::TriggerResults>           (ps.getParameter<InputTag>("patTriggerResults"));
//  BadChCandFilterToken_    = consumes<bool>                          (ps.getParameter<InputTag>("BadChargedCandidateFilter"));
//  BadPFMuonFilterToken_    = consumes<bool>                          (ps.getParameter<edm::InputTag>("BadPFMuonFilter"));
//  pfMETlabel_              = consumes<View<pat::MET> >               (ps.getParameter<InputTag>("pfMETLabel"));
//  AODCaloMETlabel_         = consumes<edm::View<reco::CaloMET> >     (ps.getParameter<InputTag>("AODCaloMETlabel"));
//  AODpfChMETlabel_         = consumes<edm::View<reco::PFMET> >       (ps.getParameter<InputTag>("AODpfChMETlabel"));
//  AODpfMETlabel_           = consumes<edm::View<reco::PFMET> >       (ps.getParameter<InputTag>("AODpfMETlabel"));
//
//  // muons
//  muonCollection_          = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("muonSrc"));
//  muonAODCollection_       = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("muonAODSrc"));
//
//  // photons
//  photonCollection_        = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonSrc"));
//  //photonAODCollection_     = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonAODSrc"));
//  photonAODCollection_     =  mayConsume<edm::View<reco::Photon> >           (ps.getParameter<InputTag>("photonAODSrc"));
//
//  // Photon ID in VID framwork 
//  phoLooseIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoLooseIdMap"));
//  phoMediumIdMapToken_            = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoMediumIdMap"));
//  phoTightIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoTightIdMap"));
//  //phoMVAValuesMapToken_           = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoMVAValuesMap")); 
//  phoChargedIsolationToken_       = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoChargedIsolation"));
//  phoNeutralHadronIsolationToken_ = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoNeutralHadronIsolation"));
//  phoPhotonIsolationToken_        = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoPhotonIsolation"));
//  phoWorstChargedIsolationToken_  = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoWorstChargedIsolation"));
//
//  //AOD Photon ID
//  AOD_phoLooseIdLabel_ = ps.getParameter<edm::InputTag>("AOD_phoLooseIdMap");
//  AOD_phoMediumIdLabel_ = ps.getParameter<edm::InputTag>("AOD_phoMediumIdMap");
//  AOD_phoTightIdLabel_ = ps.getParameter<edm::InputTag>("AOD_phoTightIdMap");
//  AOD_phoChargedIsolationLabel_ = ps.getParameter<edm::InputTag>("AOD_phoChargedIsolationMap");
//  AOD_phoNeutralHadronIsolationLabel_ = ps.getParameter<edm::InputTag>("AOD_phoNeutralHadronIsolationMap");
//  AOD_phoPhotonIsolationLabel_ = ps.getParameter<edm::InputTag>("AOD_phoPhotonIsolationMap");
//  AOD_phoWorstChargedIsolationLabel_ = ps.getParameter<edm::InputTag>("AOD_phoWorstChargedIsolationMap");
//  //
//  AOD_phoLooseIdMapToken_ = consumes<edm::ValueMap<bool> >(AOD_phoLooseIdLabel_);
//  AOD_phoMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(AOD_phoMediumIdLabel_);
//  AOD_phoTightIdMapToken_ = consumes<edm::ValueMap<bool> >(AOD_phoTightIdLabel_);
//  AOD_phoChargedIsolationMapToken_ = consumes<edm::ValueMap<float> >(AOD_phoChargedIsolationLabel_);
//  AOD_phoNeutralHadronIsolationMapToken_ = consumes<edm::ValueMap<float> >(AOD_phoNeutralHadronIsolationLabel_);
//  AOD_phoPhotonIsolationMapToken_ = consumes<edm::ValueMap<float> >(AOD_phoPhotonIsolationLabel_);
//  AOD_phoWorstChargedIsolationMapToken_ = consumes<edm::ValueMap<float> >(AOD_phoWorstChargedIsolationLabel_);
//
//
//
//  // trigger
//  triggerBits_                    = consumes <edm::TriggerResults>                     (ps.getParameter<edm::InputTag>("bits"));
//  triggerObjects_                 = consumes <edm::View<pat::TriggerObjectStandAlone>> (ps.getParameter<edm::InputTag>("objects"));
//  triggerPrescales_               = consumes <pat::PackedTriggerPrescales>             (ps.getParameter<edm::InputTag>("prescales"));
//  //AOD
//  AODTriggerLabel_                =                                                     ps.getParameter<InputTag>("AODTriggerInputTag");
//  AODTriggerToken_                = consumes<edm::TriggerResults>(AODTriggerLabel_);
//  AODTriggerEventLabel_           =                                                     ps.getParameter<InputTag>("AODTriggerEventInputTag");
//  AODTriggerEventToken_           = consumes<trigger::TriggerEvent>(AODTriggerEventLabel_);
//
//  // gen
//  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));

  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  hTTSF_   = fs->make<TH1F>("hTTSF",      "TTbar scalefactors",   200,  0,   2);
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed events",   1,  0,   2);

 if(doLLDJAOD_ || doMiniAOD_){
  branchesMiniAODTrigger     (tree_);
  branchesMiniAODGlobalEvent (tree_);
  branchesMiniAODElectrons   (tree_);
  branchesMiniAODPhotons     (tree_);
  branchesMiniAODMuons       (tree_);
  branchesMiniAODMET         (tree_);
  //branchesMiniAODGenPart     (tree_);

 }

// if(doLLDJAOD_){
//  void branchesLLDJAODJets        (TTree*);
// }

//  void branchesAODJets            (TTree*);
//  void branchesMiniAODJets        (TTree*);



// if(doMiniAOD_){
//  // make branches for tree
//  branchesGlobalEvent(tree_);
//  branchesTrigger(tree_);
//  branchesPhotons(tree_);
//  branchesElectrons(tree_);
//  branchesMuons(tree_);
//
//  branchesJets(tree_);   
//  branchesMET(tree_);
// }
// if(doAOD_){
//  branchesAODEvent(tree_);
//  branchesGenPart(tree_);
//  branchesAODTrigger(tree_);
//  branchesAODJets(tree_);
//  branchesAODMuons(tree_);
//  branchesAODPhotons(tree_);
//  branchesAODElectrons(tree_);
//  branchesAODMET(tree_);
// }

}

lldjNtuple::~lldjNtuple() {
}


void lldjNtuple::beginRun(edm::Run const& run, edm::EventSetup const& eventsetup) {

  bool changed(true);
  if(hltConfig_.init(run,eventsetup,"HLT",changed)){
  }
  hltPrescale_.init(run,eventsetup,"HLT",changed);

}

void lldjNtuple::analyze(const edm::Event& e, const edm::EventSetup& es) {

 if(doLLDJAOD_ || doMiniAOD_){
  //# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution
  slimmedJetResolution_   = JME::JetResolution::get(es, "AK4PFchs_pt");
  slimmedJetResolutionSF_ = JME::JetResolutionScaleFactor::get(es, "AK4PFchs");

  fillMiniAODTrigger     (e, es);
  fillMiniAODGlobalEvent (e, es);
  fillMiniAODElectrons   (e, es);
  fillMiniAODPhotons     (e, es);
  fillMiniAODMuons       (e, es);
  fillMiniAODMET         (e, es);
  //fillMiniAODGenPart     (e, es);

 }

// if(doLLDJAOD_){
//  void branchesLLDJAODJets        (TTree*);
// }

//  void branchesAODJets            (TTree*);
//  void branchesMiniAODJets        (TTree*);

// 
//  fillGlobalEvent(e, es);
//  fillTrigger(e, es);
//  fillPhotons(e, es);
//  fillElectrons(e, es);
// 
//  // muons use vtx for isolation
//  edm::Handle<reco::VertexCollection> vtxHandle;
//  e.getByToken(vtxLabel_, vtxHandle);
//  reco::Vertex vtx;
//  // best-known primary vertex coordinates
//  math::XYZPoint pv(0, 0, 0); 
//  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
//    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
//    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
//    bool isFake = -(v->chi2() == 0 && v->ndof() == 0); 
// 
//    if (!isFake) {
//      pv.SetXYZ(v->x(), v->y(), v->z());
//      vtx = *v; 
//      break;
//    }   
//  }
//  fillMuons(e, vtx); //muons use vtx for isolation
// 
//  fillJets(e,es);
//  fillMET(e, es);
// }
//
// if(doAOD_){
//  fillAODEvent(e, es);
//  if (!e.isRealData()) fillGenPart(e);
//  fillAODTrigger(e, es);
//  fillAODJets(e, es);
//  fillAODPhotons(e, es);
//
//  //Vertex for Muon 
//  edm::Handle<edm::View<reco::Vertex> > vtxHandle;
//  e.getByToken(AODVertexLabel_, vtxHandle);
//  reco::Vertex vtx;
//  for (edm::View<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
//    if (!v->isFake()) {
//      vtx = *v; 
//      break;
//    }   
//  }
//  fillAODElectrons(e, es, vtx);
//  fillAODMuons(e, vtx); 
//  fillAODMET(e, es);
//
// }

 hEvents_->Fill(1.);
 tree_->Fill();
}

DEFINE_FWK_MODULE(lldjNtuple);
