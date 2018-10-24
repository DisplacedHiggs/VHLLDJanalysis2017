#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "VHLLDJanalysis2017/ntuples/interface/lldjNtuple.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "RecoTracker/DebugTools/interface/GetTrackTrajInfo.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "PhysicsTools/RecoUtils/interface/CheckHitPattern.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
#include "RecoVertex/VertexTools/interface/VertexCompatibleWithBeam.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertError.h"

#include "DataFormats/PatCandidates/interface/Jet.h"



using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// AOD Collection Handles
edm::Handle<edm::View<reco::CaloJet> >  AODak4CaloJetsHandle;   
edm::Handle<edm::View<reco::PFJet>   >  AODak4PFJetsHandle;     
edm::Handle<edm::View<reco::PFJet>   >  AODak4PFJetsCHSHandle;  
//edm::Handle<edm::View<pat::Jet>      >  AODselectedPatJetsHandle;  
edm::Handle<edm::View<reco::Vertex>  >  AODVertexHandle;
edm::Handle<edm::View<reco::Track>   >  AODTrackHandle;
edm::Handle<reco::BeamSpot> beamspotHandle_;
edm::ESHandle<MagneticField> magneticField;
// transient tracks
map<reco::TransientTrack,reco::TrackBaseRef> refMap;
vector<TrajectoryStateOnSurface> tsosList;

// Global stuff for displaced vertices
ConfigurableVertexReconstructor* vtxfitter_; 
VertexDistanceXY vertexDistanceXY_;
VertexCompatibleWithBeam* vertexBeam_ = new VertexCompatibleWithBeam(vertexDistanceXY_,100);

// miniAOD Collection Handles

// AOD ---------------------------------------------
// Calo Jets
Int_t          AOD_nCaloJet_;
vector<float>  AOD_CaloJetPt_;
vector<float>  AOD_CaloJetEta_;
vector<float>  AOD_CaloJetPhi_;

vector<float>  AOD_CaloJetAlphaMax_;
vector<float>  AOD_CaloJetAlphaMax2_;
vector<float>  AOD_CaloJetAlphaMaxPrime_;
vector<float>  AOD_CaloJetAlphaMaxPrime2_;
vector<float>  AOD_CaloJetBeta_;
vector<float>  AOD_CaloJetBeta2_;

vector<float>  AOD_CaloJetSumIP_;
vector<float>  AOD_CaloJetSumIPSig_;
vector<float>  AOD_CaloJetLog10IPSig_;
vector<float>  AOD_CaloJetMedianIP_;
vector<float>  AOD_CaloJetMedianLog10IPSig_;
vector<float>  AOD_CaloJetTrackAngle_;
vector<float>  AOD_CaloJetLogTrackAngle_;
vector<float>  AOD_CaloJetMedianLog10TrackAngle_;
vector<float>  AOD_CaloJetTotalTrackAngle_;

vector<float>  AOD_CaloJetAvfVx_;
vector<float>  AOD_CaloJetAvfVy_;
vector<float>  AOD_CaloJetAvfVz_;
vector<float>  AOD_CaloJetAvfVertexTotalChiSquared_;
vector<float>  AOD_CaloJetAvfVertexDegreesOfFreedom_;
vector<float>  AOD_CaloJetAvfVertexChi2NDoF_;
vector<float>  AOD_CaloJetAvfVertexDistanceToBeam_;
vector<float>  AOD_CaloJetAvfVertexTransverseError_;
vector<float>  AOD_CaloJetAvfVertexTransverseSig_;
vector<float>  AOD_CaloJetAvfVertexDeltaEta_;
vector<float>  AOD_CaloJetAvfVertexDeltaPhi_;
vector<float>  AOD_CaloJetAvfVertexRecoilPt_;
vector<float>  AOD_CaloJetAvfVertexTrackMass_;
vector<float>  AOD_CaloJetAvfVertexTrackEnergy_;
vector<float>  AOD_CaloJetAvfBeamSpotDeltaPhi_;
vector<float>  AOD_CaloJetAvfBeamSpotRecoilPt_;
vector<float>  AOD_CaloJetAvfBeamSpotMedianDeltaPhi_;
vector<float>  AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_;
vector<int>    AOD_CaloJetNMatchedTracks_;
vector<int>    AOD_CaloJetNCleanMatchedTracks_;
vector<int>    AOD_CaloJetSumHitsInFrontOfVert_;
vector<int>    AOD_CaloJetSumMissHitsAfterVert_;
vector<int>    AOD_CaloJetHitsInFrontOfVertPerTrack_;
vector<int>    AOD_CaloJetMissHitsAfterVertPerTrack_;
vector<float>  AOD_CaloJetAvfDistToPV_;
vector<float>  AOD_CaloJetAvfVertexDeltaZtoPV_;
vector<float>  AOD_CaloJetAvfVertexDeltaZtoPV2_;

// PAT Jets
Int_t          AOD_nPATJet_;
vector<int>    AOD_PATJetPartonFlavour_;
vector<float>  AOD_PATJetPt_;
vector<float>  AOD_PATJetEta_;
vector<float>  AOD_PATJetPhi_;
vector<float>  AOD_PATJetCSV_;
vector<float>  AOD_PATJetMVA_;

// PF Jets
Int_t          AOD_nPFJet_;
vector<int>    AOD_PFJetID_; 
vector<float>  AOD_PFJetPt_;
vector<float>  AOD_PFJetEta_;
vector<float>  AOD_PFJetPhi_;
vector<float>  AOD_PFJetAlphaMax_;
vector<float>  AOD_PFJetAlphaMax2_;
vector<float>  AOD_PFJetAlphaMaxPrime_;
vector<float>  AOD_PFJetAlphaMaxPrime2_;
vector<float>  AOD_PFJetBeta_;
vector<float>  AOD_PFJetBeta2_;
vector<float>  AOD_PFJetSumIP_;
vector<float>  AOD_PFJetSumIPSig_;
vector<float>  AOD_PFJetMedianIP_;
vector<float>  AOD_PFJetMedianLog10IPSig_;
vector<float>  AOD_PFJetTrackAngle_;
vector<float>  AOD_PFJetLogTrackAngle_;
vector<float>  AOD_PFJetMedianLog10TrackAngle_;
vector<float>  AOD_PFJetTotalTrackAngle_;

// PFchs Jets
Int_t          AOD_nPFchsJet_;
vector<int>    AOD_PFchsJetID_; 
vector<float>  AOD_PFchsJetPt_;
vector<float>  AOD_PFchsJetEta_;
vector<float>  AOD_PFchsJetPhi_;
vector<float>  AOD_PFchsJetAlphaMax_;
vector<float>  AOD_PFchsJetAlphaMax2_;
vector<float>  AOD_PFchsJetAlphaMaxPrime_;
vector<float>  AOD_PFchsJetAlphaMaxPrime2_;
vector<float>  AOD_PFchsJetBeta_;
vector<float>  AOD_PFchsJetBeta2_;
vector<float>  AOD_PFchsJetSumIP_;
vector<float>  AOD_PFchsJetSumIPSig_;
vector<float>  AOD_PFchsJetMedianIP_;
vector<float>  AOD_PFchsJetMedianLog10IPSig_;
vector<float>  AOD_PFchsJetTrackAngle_;
vector<float>  AOD_PFchsJetLogTrackAngle_;
vector<float>  AOD_PFchsJetMedianLog10TrackAngle_;
vector<float>  AOD_PFchsJetTotalTrackAngle_;

// temporary variables

vector<TVector3> AOD_allTrackPositions; // x,y,z 

vector<float>    AOD_allTrackPt;
vector<float>    AOD_allTrackEta;
vector<float>    AOD_allTrackPhi;
vector<float>    AOD_allTrackIFSPt;

vector<int>      AOD_allTracknMissingInner;
vector<int>      AOD_allTracknMissingOuter;
vector<float>    AOD_allTrackAngle;

vector<int>      AOD_whichVertexByTrack;

vector<float>    AOD_allTrackdxy;
vector<float>    AOD_allTrackdxyerr;

// set parameters for tracks to be accepted
const float minTrackPt_ = 1.0;
const float maxDRtrackJet_ = 0.4;


// start functions ----------------------------------------

// initialize branchesMiniAOD
void lldjNtuple::branchesLLDJAODJets(TTree* tree) {

  tree->Branch("AOD_nCaloJet"                   , &AOD_nCaloJet_);
  tree->Branch("AOD_CaloJetPt"                  , &AOD_CaloJetPt_);
  tree->Branch("AOD_CaloJetEta"                 , &AOD_CaloJetEta_);
  tree->Branch("AOD_CaloJetPhi"                 , &AOD_CaloJetPhi_);

  tree->Branch("AOD_CaloJetAlphaMax"            , &AOD_CaloJetAlphaMax_);                               
  tree->Branch("AOD_CaloJetAlphaMax2"           , &AOD_CaloJetAlphaMax2_);                               
  tree->Branch("AOD_CaloJetAlphaMaxPrime"       , &AOD_CaloJetAlphaMaxPrime_);                               
  tree->Branch("AOD_CaloJetAlphaMaxPrime2"      , &AOD_CaloJetAlphaMaxPrime2_);                               
  tree->Branch("AOD_CaloJetBeta"                , &AOD_CaloJetBeta_);                               
  tree->Branch("AOD_CaloJetBeta2"               , &AOD_CaloJetBeta2_);                               

  tree->Branch("AOD_CaloJetSumIP"               , &AOD_CaloJetSumIP_);
  tree->Branch("AOD_CaloJetSumIPSig"            , &AOD_CaloJetSumIPSig_);
  tree->Branch("AOD_CaloJetMedianIP"            , &AOD_CaloJetMedianIP_);
  tree->Branch("AOD_CaloJetMedianLog10IPSig"    , &AOD_CaloJetMedianLog10IPSig_);
  tree->Branch("AOD_CaloJetTrackAngle"          , &AOD_CaloJetTrackAngle_);
  tree->Branch("AOD_CaloJetLogTrackAngle"       , &AOD_CaloJetLogTrackAngle_);
  tree->Branch("AOD_CaloJetMedianLog10TrackAngle" , &AOD_CaloJetMedianLog10TrackAngle_);
  tree->Branch("AOD_CaloJetTotalTrackAngle"     , &AOD_CaloJetTotalTrackAngle_);
                       
  tree->Branch("AOD_CaloJetAvfVx", &AOD_CaloJetAvfVx_);
  tree->Branch("AOD_CaloJetAvfVy", &AOD_CaloJetAvfVy_);
  tree->Branch("AOD_CaloJetAvfVz", &AOD_CaloJetAvfVz_);
  tree->Branch("AOD_CaloJetAvfVertexTotalChiSquared",     &AOD_CaloJetAvfVertexTotalChiSquared_);
  tree->Branch("AOD_CaloJetAvfVertexDegreesOfFreedom",    &AOD_CaloJetAvfVertexDegreesOfFreedom_);
  tree->Branch("AOD_CaloJetAvfVertexChi2NDoF",            &AOD_CaloJetAvfVertexChi2NDoF_);
  tree->Branch("AOD_CaloJetAvfVertexDistanceToBeam",      &AOD_CaloJetAvfVertexDistanceToBeam_);
  tree->Branch("AOD_CaloJetAvfVertexTransverseError",     &AOD_CaloJetAvfVertexTransverseError_);
  tree->Branch("AOD_CaloJetAvfVertexTransverseSig",       &AOD_CaloJetAvfVertexTransverseSig_);
  tree->Branch("AOD_CaloJetAvfVertexDeltaEta",            &AOD_CaloJetAvfVertexDeltaEta_);
  tree->Branch("AOD_CaloJetAvfVertexDeltaPhi",            &AOD_CaloJetAvfVertexDeltaPhi_);
  tree->Branch("AOD_CaloJetAvfVertexRecoilPt",            &AOD_CaloJetAvfVertexRecoilPt_);
  tree->Branch("AOD_CaloJetAvfVertexTrackMass",           &AOD_CaloJetAvfVertexTrackMass_);
  tree->Branch("AOD_CaloJetAvfVertexTrackEnergy",         &AOD_CaloJetAvfVertexTrackEnergy_);
  tree->Branch("AOD_CaloJetAvfBeamSpotDeltaPhi",          &AOD_CaloJetAvfBeamSpotDeltaPhi_);
  tree->Branch("AOD_CaloJetAvfBeamSpotRecoilPt",          &AOD_CaloJetAvfBeamSpotRecoilPt_);
  tree->Branch("AOD_CaloJetAvfBeamSpotMedianDeltaPhi",    &AOD_CaloJetAvfBeamSpotMedianDeltaPhi_);
  tree->Branch("AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi", &AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_);
  tree->Branch("AOD_CaloJetNCleanMatchedTracks",            &AOD_CaloJetNCleanMatchedTracks_);
  tree->Branch("AOD_CaloJetNMatchedTracks"      ,           &AOD_CaloJetNMatchedTracks_);
  tree->Branch("AOD_CaloJetSumHitsInFrontOfVert",           &AOD_CaloJetSumHitsInFrontOfVert_);
  tree->Branch("AOD_CaloJetSumMissHitsAfterVert",           &AOD_CaloJetSumMissHitsAfterVert_);
  tree->Branch("AOD_CaloJetHitsInFrontOfVertPerTrack",      &AOD_CaloJetHitsInFrontOfVertPerTrack_);
  tree->Branch("AOD_CaloJetMissHitsAfterVertPerTrack",      &AOD_CaloJetMissHitsAfterVertPerTrack_);
  tree->Branch("AOD_CaloJetAvfDistToPV",                    &AOD_CaloJetAvfDistToPV_);
  tree->Branch("AOD_CaloJetAvfVertexDeltaZtoPV",            &AOD_CaloJetAvfVertexDeltaZtoPV_);
  tree->Branch("AOD_CaloJetAvfVertexDeltaZtoPV2",           &AOD_CaloJetAvfVertexDeltaZtoPV2_);
                                                                
  tree->Branch("AOD_nPATJet",              &AOD_nPATJet_);
  tree->Branch("AOD_PATJetPartonFlavour",  &AOD_PATJetPartonFlavour_);
  tree->Branch("AOD_PATJetPt",             &AOD_PATJetPt_);
  tree->Branch("AOD_PATJetEta",            &AOD_PATJetEta_);
  tree->Branch("AOD_PATJetPhi",            &AOD_PATJetPhi_);
  tree->Branch("AOD_PATJetCSV",            &AOD_PATJetCSV_);
  tree->Branch("AOD_PATJetMVA",            &AOD_PATJetMVA_);
                          
  // PF Jets
  tree->Branch("AOD_nPFJet"                        , &AOD_nPFJet_);                                 
  tree->Branch("AOD_PFJetID"                       , &AOD_PFJetID_);                                 
  tree->Branch("AOD_PFJetPt"                       , &AOD_PFJetPt_);                                 
  tree->Branch("AOD_PFJetEta"                      , &AOD_PFJetEta_);                                 
  tree->Branch("AOD_PFJetPhi"                      , &AOD_PFJetPhi_);                                 
  tree->Branch("AOD_PFJetAlphaMax"                 , &AOD_PFJetAlphaMax_);                                 
  tree->Branch("AOD_PFJetAlphaMax2"                , &AOD_PFJetAlphaMax2_);                                 
  tree->Branch("AOD_PFJetAlphaMaxPrime"            , &AOD_PFJetAlphaMaxPrime_);                                 
  tree->Branch("AOD_PFJetAlphaMaxPrime2"           , &AOD_PFJetAlphaMaxPrime2_);                                 
  tree->Branch("AOD_PFJetBeta"                     , &AOD_PFJetBeta_);                                 
  tree->Branch("AOD_PFJetBeta2"                    , &AOD_PFJetBeta2_);                                 
  tree->Branch("AOD_PFJetSumIP"                    , &AOD_PFJetSumIP_);                                 
  tree->Branch("AOD_PFJetSumIPSig"                 , &AOD_PFJetSumIPSig_);                                 
  tree->Branch("AOD_PFJetMedianIP"                 , &AOD_PFJetMedianIP_);                                                                           
  tree->Branch("AOD_PFJetMedianLog10IPSig"         , &AOD_PFJetMedianLog10IPSig_);                           
  tree->Branch("AOD_PFJetTrackAngle"               , &AOD_PFJetTrackAngle_);                                    
  tree->Branch("AOD_PFJetLogTrackAngle"            , &AOD_PFJetLogTrackAngle_);                                         
  tree->Branch("AOD_PFJetMedianLog10TrackAngle"    , &AOD_PFJetMedianLog10TrackAngle_);                           
  tree->Branch("AOD_PFJetTotalTrackAngle"          , &AOD_PFJetTotalTrackAngle_);
  
  // PFchs Jets 
  tree->Branch("AOD_nPFchsJet"                     , &AOD_nPFchsJet_);                                  
  tree->Branch("AOD_PFchsJetID"                    , &AOD_PFchsJetID_);                                  
  tree->Branch("AOD_PFchsJetPt"                    , &AOD_PFchsJetPt_);                                  
  tree->Branch("AOD_PFchsJetEta"                   , &AOD_PFchsJetEta_);                                 
  tree->Branch("AOD_PFchsJetPhi"                   , &AOD_PFchsJetPhi_);                                      
  tree->Branch("AOD_PFchsJetAlphaMax"              , &AOD_PFchsJetAlphaMax_);                                  
  tree->Branch("AOD_PFchsJetAlphaMax2"             , &AOD_PFchsJetAlphaMax2_);                                     
  tree->Branch("AOD_PFchsJetAlphaMaxPrime"         , &AOD_PFchsJetAlphaMaxPrime_);                                  
  tree->Branch("AOD_PFchsJetAlphaMaxPrime2"        , &AOD_PFchsJetAlphaMaxPrime2_);                       
  tree->Branch("AOD_PFchsJetBeta"                  , &AOD_PFchsJetBeta_);                                  
  tree->Branch("AOD_PFchsJetBeta2"                 , &AOD_PFchsJetBeta2_);                                 
  tree->Branch("AOD_PFchsJetSumIP"                 , &AOD_PFchsJetSumIP_);                                    
  tree->Branch("AOD_PFchsJetSumIPSig"              , &AOD_PFchsJetSumIPSig_);                                 
  tree->Branch("AOD_PFchsJetMedianIP"              , &AOD_PFchsJetMedianIP_);                                   
  tree->Branch("AOD_PFchsJetMedianLog10IPSig"      , &AOD_PFchsJetMedianLog10IPSig_);                              
  tree->Branch("AOD_PFchsJetTrackAngle"            , &AOD_PFchsJetTrackAngle_);                                            
  tree->Branch("AOD_PFchsJetLogTrackAngle"         , &AOD_PFchsJetLogTrackAngle_);                                   
  tree->Branch("AOD_PFchsJetMedianLog10TrackAngle" , &AOD_PFchsJetMedianLog10TrackAngle_);
  tree->Branch("AOD_PFchsJetTotalTrackAngle"       , &AOD_PFchsJetTotalTrackAngle_);

}

//fills slimmedJets .clear() to empty vector of old data
void lldjNtuple::fillLLDJAODJets(const edm::Event& e, const edm::EventSetup& es) {

 // bool dodebug = false;
 // cleanup from previous execution

 AOD_nCaloJet_=0;
 AOD_CaloJetPt_.clear();
 AOD_CaloJetEta_.clear();
 AOD_CaloJetPhi_.clear();
       //AODnTracksToCaloJet_.clear();
 AOD_CaloJetNMatchedTracks_.clear();

 AOD_CaloJetAlphaMax_.clear();                               
 AOD_CaloJetAlphaMax2_.clear();                               
 AOD_CaloJetAlphaMaxPrime_.clear();                               
 AOD_CaloJetAlphaMaxPrime2_.clear();                               
 AOD_CaloJetBeta_.clear();                               
 AOD_CaloJetBeta2_.clear();                               

 AOD_CaloJetSumIP_.clear();
 AOD_CaloJetSumIPSig_.clear();
 AOD_CaloJetMedianIP_.clear();
 AOD_CaloJetMedianLog10IPSig_.clear();
 AOD_CaloJetTrackAngle_.clear();
 AOD_CaloJetLogTrackAngle_.clear();
 AOD_CaloJetMedianLog10TrackAngle_.clear();
 AOD_CaloJetTotalTrackAngle_.clear();

 AOD_CaloJetAvfVx_.clear();
 AOD_CaloJetAvfVy_.clear();
 AOD_CaloJetAvfVz_.clear();
 AOD_CaloJetAvfVertexTotalChiSquared_.clear();
 AOD_CaloJetAvfVertexDegreesOfFreedom_.clear();
 AOD_CaloJetAvfVertexChi2NDoF_.clear();
 AOD_CaloJetAvfVertexDistanceToBeam_.clear();
 AOD_CaloJetAvfVertexTransverseError_.clear();
 AOD_CaloJetAvfVertexTransverseSig_.clear();
 AOD_CaloJetAvfVertexDeltaEta_.clear();
 AOD_CaloJetAvfVertexDeltaPhi_.clear();
 AOD_CaloJetAvfVertexRecoilPt_.clear();
 AOD_CaloJetAvfVertexTrackMass_.clear();
 AOD_CaloJetAvfVertexTrackEnergy_.clear();
 AOD_CaloJetAvfBeamSpotDeltaPhi_.clear();
 AOD_CaloJetAvfBeamSpotRecoilPt_.clear();
 AOD_CaloJetAvfBeamSpotMedianDeltaPhi_.clear();
 AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_.clear();
 AOD_CaloJetNCleanMatchedTracks_.clear();
 AOD_CaloJetSumHitsInFrontOfVert_.clear();
 AOD_CaloJetSumMissHitsAfterVert_.clear();
 AOD_CaloJetHitsInFrontOfVertPerTrack_.clear();
 AOD_CaloJetMissHitsAfterVertPerTrack_.clear();
 AOD_CaloJetAvfDistToPV_.clear();
 AOD_CaloJetAvfVertexDeltaZtoPV_.clear();
 AOD_CaloJetAvfVertexDeltaZtoPV2_.clear();

 // PAT Jets
 AOD_nPATJet_ = 0;
 AOD_PATJetPartonFlavour_.clear();
 AOD_PATJetPt_.clear();
 AOD_PATJetEta_.clear();
 AOD_PATJetPhi_.clear();
 AOD_PATJetCSV_.clear();
 AOD_PATJetMVA_.clear();

 // PF Jets
 AOD_nPFJet_=0;
 AOD_PFJetID_.clear();
 AOD_PFJetPt_.clear();
 AOD_PFJetEta_.clear();
 AOD_PFJetPhi_.clear();
 AOD_PFJetAlphaMax_.clear();
 AOD_PFJetAlphaMax2_.clear();
 AOD_PFJetAlphaMaxPrime_.clear();
 AOD_PFJetAlphaMaxPrime2_.clear();
 AOD_PFJetBeta_.clear();
 AOD_PFJetBeta2_.clear();
 AOD_PFJetSumIP_.clear();
 AOD_PFJetSumIPSig_.clear();
 AOD_PFJetMedianIP_.clear();
 AOD_PFJetMedianLog10IPSig_.clear();
 AOD_PFJetTrackAngle_.clear();
 AOD_PFJetLogTrackAngle_.clear();
 AOD_PFJetMedianLog10TrackAngle_.clear();
 AOD_PFJetTotalTrackAngle_.clear();
 
 // PFchs Jets
 AOD_nPFchsJet_=0;
 AOD_PFchsJetID_.clear();
 AOD_PFchsJetPt_.clear();
 AOD_PFchsJetEta_.clear();
 AOD_PFchsJetPhi_.clear();
 AOD_PFchsJetAlphaMax_.clear();
 AOD_PFchsJetAlphaMax2_.clear();
 AOD_PFchsJetAlphaMaxPrime_.clear();
 AOD_PFchsJetAlphaMaxPrime2_.clear();
 AOD_PFchsJetBeta_.clear();
 AOD_PFchsJetBeta2_.clear();
 AOD_PFchsJetSumIP_.clear();
 AOD_PFchsJetSumIPSig_.clear();
 AOD_PFchsJetMedianIP_.clear();
 AOD_PFchsJetMedianLog10IPSig_.clear();
 AOD_PFchsJetTrackAngle_.clear();
 AOD_PFchsJetLogTrackAngle_.clear();
 AOD_PFchsJetMedianLog10TrackAngle_.clear();
 AOD_PFchsJetTotalTrackAngle_.clear();
 
 //e.getByToken(rhoLabel_, rhoHandle);
 //float rho = *(rhoHandle.product());
 
// e.getByToken(vtxLabel_, vtxHandle);
// if (!vtxHandle.isValid()) edm::LogWarning("lldjNtuple") << "Primary vertices info not unavailable";

// // Accessing the JEC uncertainties 
// //ak4  
// edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
// es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
// JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
// JetCorrectionUncertainty *jecUnc=0;
// jecUnc = new JetCorrectionUncertainty(JetCorPar);
// // //ak8
// // edm::ESHandle<JetCorrectorParametersCollection> AK8JetCorParColl;
// // es.get<JetCorrectionsRecord>().get("AK8PFchs",AK8JetCorParColl);
// // JetCorrectorParameters const & AK8JetCorPar = (*AK8JetCorParColl)["Uncertainty"];
// // JetCorrectionUncertainty *AK8jecUnc=0;
// // AK8jecUnc = new JetCorrectionUncertainty(AK8JetCorPar);
 
 // ----------------------------------------------------------
 // AOD Section ----------------------------------------------
 
 bool verbose_AOD = false;
 //bool verbose_AOD = true;
 
 // AOD Jet Handles
 e.getByToken( AODak4CaloJetsLabel_ ,  AODak4CaloJetsHandle  );   
 e.getByToken( AODak4PFJetsLabel_   ,  AODak4PFJetsHandle    );     
 e.getByToken( AODak4PFJetsCHSLabel_,  AODak4PFJetsCHSHandle );  
 //e.getByToken( selectedPatJetsLabel_,  selectedPatJetsHandle );  
 e.getByToken( AODVertexLabel_      ,  AODVertexHandle );
 e.getByToken( AODTrackLabel_       ,  AODTrackHandle );

 // Magnetic field
 es.get<IdealMagneticFieldRecord>().get(magneticField);
 magneticField_ = &*magneticField;

 // beamspot
 e.getByToken(beamspotLabel_, beamspotHandle_);

 // clear values
 int nMissingInner = 0;
 int nMissingOuter = 0; 

 // propagator
 std::string thePropagatorName_ = "PropagatorWithMaterial";
 es.get<TrackingComponentsRecord>().get(thePropagatorName_,thePropagator_);
 StateOnTrackerBound stateOnTracker(thePropagator_.product());

 // Vertex
 vector<int> whichVertex_;
 whichVertex_.clear();
 whichVertex_ = vector<int>(AODTrackHandle->size(),-1);

 vtxfitter_ = new ConfigurableVertexReconstructor(lldj_pset_.getParameter<edm::ParameterSet>("vertexFitterConfig"));

 // clear master track vectors before starting track loop
 AOD_allTrackPositions       .clear(); 
 AOD_allTrackPt              .clear(); 
 AOD_allTrackEta             .clear(); 
 AOD_allTrackPhi             .clear(); 
 AOD_allTrackIFSPt           .clear(); 
 AOD_allTracknMissingInner   .clear(); 
 AOD_allTracknMissingOuter   .clear(); 
 AOD_allTrackAngle           .clear(); 
 AOD_whichVertexByTrack      .clear(); 
 AOD_allTrackdxy             .clear(); 
 AOD_allTrackdxyerr          .clear(); 

 for(int j = 0; j < (int)AODTrackHandle->size(); j++){


  // get track j using the AOD track handle 
  reco::TrackBaseRef tref(AODTrackHandle,j);
  // make transient track (unfolding effects of B field ?)
  reco::TransientTrack tt(AODTrackHandle->at(j),magneticField_);

  if(!tt.isValid())continue;

  // track pt first
  float trackpt  = tref->pt();

  // make selections on track
  // for alphaMax, track angle we use ttIFSpt, not tref->pt()
  //   ---------!!!!--------------
  if (trackpt < minTrackPt_) continue;  // minimum pT for track
  if (!tref->quality(reco::TrackBase::highPurity)) continue; // track must be highPurity

  // find where track (no B field) would hit outer tracker
  FreeTrajectoryState fts = trajectoryStateTransform::initialFreeState(AODTrackHandle->at(j),magneticField_);
  TrajectoryStateOnSurface outer = stateOnTracker(fts);
  if(!outer.isValid()) continue;
  GlobalPoint outerPos = outer.globalPosition();
  TVector3 trackPos(outerPos.x(),outerPos.y(),outerPos.z());

  // push back track position to save in master vector
  AOD_allTrackPositions.push_back(trackPos);

  // track basics (trackpt above)
  float tracketa = tref->eta();
  float trackphi = tref->phi();
  float ttIFSpt  = tt.initialFreeState().momentum().transverse();
  AOD_allTrackPt .push_back(trackpt );  
  AOD_allTrackEta.push_back(tracketa);
  AOD_allTrackPhi.push_back(trackphi);
  AOD_allTrackIFSPt.push_back(ttIFSpt);

  /// Find best vertex associated with this track
  // we are on track j
  // loop over vertices
  // reassign index bestk if trackWeight is new max
  float maxWeight = 0; 
  int bestk = -1;
  for(int k = 0; k < (int)AODVertexHandle->size();k++){        
   if(AODVertexHandle->at(k).trackWeight(tref) > maxWeight){  
    maxWeight = AODVertexHandle->at(k).trackWeight(tref);    
    bestk = k;                                               
   }                                                          
  }                                                            
  AOD_whichVertexByTrack.push_back(bestk); 
  
  // Number of missing tracker hits 
  nMissingInner = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_INNER_HITS);
  nMissingOuter = tref->hitPattern().numberOfLostTrackerHits(reco::HitPattern::MISSING_OUTER_HITS);
  AOD_allTracknMissingInner.push_back(nMissingInner) ;
  AOD_allTracknMissingOuter.push_back(nMissingOuter) ;

  /// For track angle
  // get track trajectory info
  static GetTrackTrajInfo getTrackTrajInfo; 
  vector<GetTrackTrajInfo::Result> trajInfo = getTrackTrajInfo.analyze(es, (*tref));
  if ( trajInfo.size() > 0 && trajInfo[0].valid) {
   // get inner tracker hit from trajectory state 
   const TrajectoryStateOnSurface& tsosInnerHit = trajInfo[0].detTSOS;

   //  here's the track angle
   // find beamspot x,y coordinates
   const reco::BeamSpot& pat_beamspot = (*beamspotHandle_);
   TVector2 bmspot(pat_beamspot.x0(),pat_beamspot.y0());
   // find track trajectory state on surface inner hit
   GlobalPoint  innerPos = tsosInnerHit.globalPosition();
   GlobalVector innerMom = tsosInnerHit.globalMomentum();
   
   // calculate the difference between inner hit and beamspot
   TVector2 sv(innerPos.x(),innerPos.y());
   TVector2 diff = (sv-bmspot);
   //cout<<"bs x: "<<bmspot.X()<<" y: "<<bmspot.Y()<<endl;
   //cout<<" sv x: "<<sv.X()<<" y: "<<sv.Y()<<endl;
   //cout<<" diff phi: "<<diff.Phi()<<endl;
   TVector2 momentum(innerMom.x(),innerMom.y());
   //cout<<" p x: "<<momentum.X()<<" y: "<<momentum.Y()<<endl;
   //cout<<" p phi: "<<momentum.Phi()<<endl;
   //cout<<" dPhi: "<<diff.DeltaPhi(momentum)<<endl;
   float ta = fabs( diff.DeltaPhi(momentum) ) ;

   AOD_allTrackAngle.push_back(ta);
  }
  else{ AOD_allTrackAngle.push_back(0); }

  // beamspot info, track impact parameter
  float dxy = fabs(tref->dxy(*beamspotHandle_));
  float dxyerr = tref->dxyError();
  if(verbose_AOD) printf(" dxy dxyerr: %0.4f %0.4f\n", dxy, dxyerr);
  AOD_allTrackdxy   .push_back(dxy   ) ;
  AOD_allTrackdxyerr.push_back(dxyerr) ;
   
 }//end track loop

 //Debug printing
 if(verbose_AOD){
  for(int j = 0; j < (int)AODTrackHandle->size(); j++){
   reco::TrackBaseRef tref(AODTrackHandle,j);
   printf("AOD track pt eta phi: %f %f %f\n",tref->pt(),tref->eta(),tref->phi());
  }
 
  printf("  AOD_allTrackPositions      %i \n",  (int)AOD_allTrackPositions    .size() ); 
  printf("  AOD_allTrackPt             %i \n",  (int)AOD_allTrackPt           .size() ); 
  printf("  AOD_allTrackEta            %i \n",  (int)AOD_allTrackEta          .size() ); 
  printf("  AOD_allTrackPhi            %i \n",  (int)AOD_allTrackPhi          .size() ); 
  printf("  AOD_allTrackIFSPt          %i \n",  (int)AOD_allTrackIFSPt        .size() ); 
  printf("  AOD_allTracknMissingInner  %i \n",  (int)AOD_allTracknMissingInner.size() ); 
  printf("  AOD_allTracknMissingOuter  %i \n",  (int)AOD_allTracknMissingOuter.size() ); 
  printf("  AOD_allTrackAngle          %i \n",  (int)AOD_allTrackAngle        .size() ); 
  printf("  AOD_whichVertexByTrack     %i \n",  (int)AOD_whichVertexByTrack   .size() ); 
  printf("  AOD_allTrackdxy            %i \n",  (int)AOD_allTrackdxy          .size() ); 
  printf("  AOD_allTrackdxyerr         %i \n",  (int)AOD_allTrackdxyerr       .size() ); 
 }


 // AOD Calo Jets -------------------------------------------
 for (edm::View<reco::CaloJet>::const_iterator iJet = AODak4CaloJetsHandle->begin(); iJet != AODak4CaloJetsHandle->end(); ++iJet) {

  if(verbose_AOD) printf("Calo %f \n",iJet->pt());
  
  float jetpt  = iJet->pt();
  float jeteta = iJet->eta();
  float jetphi = iJet->phi();

  // ID and jet selections
  bool passID = false;
  if( iJet->emEnergyFraction()>=0.0
   && iJet->emEnergyFraction()<=0.9
   && iJet->energyFractionHadronic()>=0.0
   && iJet->energyFractionHadronic()<=0.9)  passID = true; 

  if(iJet->pt()<20.0 || fabs(iJet->eta())>2.4 || !passID) continue;

  // caloJetTrackIDs is a vector of ints where each int is the 
  // index of a track passing deltaR requirement to this jet
  // out of the master track record of tracks passing basic selections
  vector<int>   caloJetTrackIDs = getJetTrackIndexs( jeteta, jetphi );
  AOD_CaloJetNMatchedTracks_.push_back( caloJetTrackIDs.size() );

  if(caloJetTrackIDs.size()<1) continue;

  if(verbose_AOD){
   printf(" AOD Jet pt eta phi: %0.1f %0.1f %0.1f\n",jetpt,jeteta,jetphi);
   for( int i=0; i<(int)AOD_allTrackPositions.size(); i++){
    printf("  allTrack %i %0.1f %0.1f %0.1f \n",i,
     AOD_allTrackPt [i] ,
     AOD_allTrackEta[i] ,
     AOD_allTrackPhi[i] );

   }
   for( int i=0; i<(int)caloJetTrackIDs.size(); i++){
    printf(" Track %i at %i \n",i,caloJetTrackIDs[i]);
    printf("  caloTrack %i=%i %0.1f %0.1f %0.1f \n",i,
     caloJetTrackIDs[i],
     AOD_allTrackPt [caloJetTrackIDs[i]] ,
     AOD_allTrackEta[caloJetTrackIDs[i]] ,
     AOD_allTrackPhi[caloJetTrackIDs[i]] );
   }
  }

  // initialize variables
  float alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2 = -1.;
  float totalTrackAngle, totalTrackAnglePt = 0.;
  float sumIP, sumIPSig = 0.;
  vector<float> caloJetTrackAngles; 
  caloJetTrackAngles.clear();
  vector<float> caloJetIPs; 
  caloJetIPs.clear();
  vector<float> caloJetIPSigs; 
  caloJetIPSigs.clear();

  // do calculations
  calculateAlphaMax(caloJetTrackIDs,alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2);
  calculateTrackAngle(caloJetTrackIDs, caloJetTrackAngles, totalTrackAngle, totalTrackAnglePt);
  calculateIP(caloJetTrackIDs, caloJetIPs, caloJetIPSigs, sumIP, sumIPSig);
  calculateDisplacedVertices(es, caloJetTrackIDs);

  // find medians
  float medianTrackAngle;
  medianTrackAngle = findMedian(caloJetTrackAngles);
  float medianIP;
  medianIP = findMedian(caloJetIPs);
  float medianIPSig;
  medianIPSig = findMedian(caloJetIPSigs);

  ////////////////////////
  // Fill tree
  /////////////////////////
  AOD_nCaloJet_++;
  
  //Pt, Eta, Phi
  AOD_CaloJetPt_.push_back(jetpt);
  AOD_CaloJetEta_.push_back(jeteta);
  AOD_CaloJetPhi_.push_back(jetphi);
  
  //AlphaMax-type variables
  AOD_CaloJetAlphaMax_       .push_back(alphaMax      ) ; 
  AOD_CaloJetAlphaMax2_      .push_back(alphaMax2     ) ; 
  AOD_CaloJetAlphaMaxPrime_  .push_back(alphaMaxPrime ) ; 
  AOD_CaloJetAlphaMaxPrime2_ .push_back(alphaMaxPrime2) ; 
  AOD_CaloJetBeta_           .push_back(beta          ) ; 
  AOD_CaloJetBeta2_          .push_back(beta2         ) ; 

  //Totals
  AOD_CaloJetSumIP_.push_back(sumIP);
  AOD_CaloJetSumIPSig_.push_back(sumIPSig);

  AOD_CaloJetTotalTrackAngle_.push_back(totalTrackAngle);    

  /////Medians
  AOD_CaloJetMedianIP_             .push_back(medianIP);
  AOD_CaloJetMedianLog10IPSig_     .push_back(log10(medianIPSig));
  AOD_CaloJetMedianLog10TrackAngle_.push_back(log10(medianTrackAngle));

 }


// // PAT jets for heavy flavor studies (ak4 pf chs)
// for (edm::View<pat::Jet>::const_iterator iJet = selectedPatJetsHandle->begin(); iJet != selectedPatJetsHandle->end(); ++iJet) {
//
//   if(iJet->pt()<15.0 || fabs(iJet->eta())>3) continue; 
//   
//   AOD_nPATJet_++;
//   AOD_PATJetPt_.push_back(iJet->pt());
//   AOD_PATJetEta_.push_back(iJet->eta());
//   AOD_PATJetPhi_.push_back(iJet->phi());
//   AOD_PATJetPartonFlavour_.push_back(iJet->partonFlavour()); 
//   AOD_PATJetCSV_.push_back(iJet->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
//   AOD_PATJetMVA_.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
// }
  
	
 // AOD PF Jets -------------------------------------------
 for (edm::View<reco::PFJet>::const_iterator iJet = AODak4PFJetsHandle->begin(); iJet != AODak4PFJetsHandle->end(); ++iJet) {
	 
// Remove code casting reco jet as PAT jet, because these will not have btag or partonFlavour
//
// if TagInfo present
//  if( patJet.hasTagInfo("pfInclusiveSecondaryVertexFinder") ) // need to omit 'TagInfos' from the label since PAT strips it away
//  {
//     std::cout<<" PF jet has Tag info"<<std::endl;
//    //const reco::CandSecondaryVertexTagInfo *candSVTagInfo = jet->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
//    //// if there is at least one reconstructed SV
//    //if( candSVTagInfo->nVertices() >= 1 ) 
//    //{
//    //  std::cout << "Found secondary vertex with a flight distance of " << candSVTagInfo->flightDistance(0).value() << " cm" << std::endl;
//    //}
//  }

  if(verbose_AOD) printf("PF %f \n",iJet->pt());
  
  //test btagging
  //std::cout << "partonFlavour " << iJet->partonFlavour() << std::endl;

  float jetpt  = iJet->pt();
  float jeteta = iJet->eta();
  float jetphi = iJet->phi();

  // ID and jet selections
  bool pfjetPassLooseID = false;
  bool pfjetPassTightID = false;

  // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016
  float pfjetNHF                 = iJet->neutralHadronEnergyFraction();
  float pfjetNEMF                = iJet->neutralEmEnergyFraction();
  float pfjetCHF                 = iJet->chargedHadronEnergyFraction();
  //float pfjetMUF                 = iJet->muonEnergyFraction();
  float pfjetCEMF                = iJet->chargedEmEnergyFraction();
  float pfjetNumConst            = iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
  float pfjetNumNeutralParticle  = iJet->neutralMultiplicity();
  float pfjetCHM                 = iJet->chargedMultiplicity(); 

  if ( fabs(jeteta) <= 2.7 ){
   pfjetPassLooseID = (pfjetNHF<0.99 && pfjetNEMF<0.99 && pfjetNumConst>1) && ((abs(jeteta)<=2.4 && pfjetCHF>0 && pfjetCHM>0 && pfjetCEMF<0.99) || abs(jeteta)>2.4) && abs(jeteta)<=2.7 ;
   pfjetPassTightID = (pfjetNHF<0.90 && pfjetNEMF<0.90 && pfjetNumConst>1) && ((abs(jeteta)<=2.4 && pfjetCHF>0 && pfjetCHM>0 && pfjetCEMF<0.99) || abs(jeteta)>2.4) && abs(jeteta)<=2.7 ;
  }
  if ( fabs(jeteta) > 2.7 && fabs(jeteta) <= 3.0 ){
   pfjetPassLooseID = (pfjetNEMF<0.90 && pfjetNumNeutralParticle>2 ) ;
   pfjetPassTightID = (pfjetNEMF<0.90 && pfjetNumNeutralParticle>2 ) ;
  }
  if ( fabs(jeteta) > 3.0 ){
   pfjetPassLooseID = (pfjetNEMF<0.90 && pfjetNumNeutralParticle>10 ) ;
   pfjetPassTightID = (pfjetNEMF<0.90 && pfjetNumNeutralParticle>10 ) ; 
  }

  Int_t AOD_PFJetIDdecision = 0;
  if (pfjetPassLooseID) AOD_PFJetIDdecision += pow(2, 1);
  if (pfjetPassTightID) AOD_PFJetIDdecision += pow(2, 2);

  // selections
  if(iJet->pt()<20.0 || fabs(iJet->eta())>2.4 || pfjetPassLooseID==0 ) continue;

  // index of a track passing deltaR requirement to this jet
  // out of the master track record of tracks passing basic selections
  vector<int>   pfJetTrackIDs = getJetTrackIndexs( jeteta, jetphi );
  if(pfJetTrackIDs.size()<1) continue;

  if(verbose_AOD){
   printf(" AOD Jet pt eta phi: %0.1f %0.1f %0.1f\n",jetpt,jeteta,jetphi);
   for( int i=0; i<(int)AOD_allTrackPositions.size(); i++){
    printf("  allTrack %i %0.1f %0.1f %0.1f \n",i,
     AOD_allTrackPt [i] ,
     AOD_allTrackEta[i] ,
     AOD_allTrackPhi[i] );

   }
   for( int i=0; i<(int)pfJetTrackIDs.size(); i++){
    printf(" Track %i at %i \n",i,pfJetTrackIDs[i]);
    printf("  pfTrack %i=%i %0.1f %0.1f %0.1f \n",i,
     pfJetTrackIDs[i],
     AOD_allTrackPt [pfJetTrackIDs[i]] ,
     AOD_allTrackEta[pfJetTrackIDs[i]] ,
     AOD_allTrackPhi[pfJetTrackIDs[i]] );
   }
  }

  // initialize variables
  float alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2 = -1.;
  float totalTrackAngle, totalTrackAnglePt = 0.;
  float sumIP, sumIPSig = 0.;
  vector<float> pfJetTrackAngles; 
  pfJetTrackAngles.clear();
  vector<float> pfJetIPs; 
  pfJetIPs.clear();
  vector<float> pfJetIPSigs; 
  pfJetIPSigs.clear();

  // do calculations
  calculateAlphaMax(pfJetTrackIDs,alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2);
  calculateTrackAngle(pfJetTrackIDs, pfJetTrackAngles, totalTrackAngle, totalTrackAnglePt);
  calculateIP(pfJetTrackIDs, pfJetIPs, pfJetIPSigs, sumIP, sumIPSig);
  //calculateDisplacedVertices(es, pfJetTrackIDs);

  // find medians
  float medianTrackAngle;
  medianTrackAngle = findMedian(pfJetTrackAngles);
  float medianIP;
  medianIP = findMedian(pfJetIPs);
  float medianIPSig;
  medianIPSig = findMedian(pfJetIPSigs);

  ////////////////////////
  // Fill tree
  /////////////////////////
  AOD_nPFJet_++;
  
  //Pt, Eta, Phi
  AOD_PFJetPt_.push_back(jetpt);
  AOD_PFJetEta_.push_back(jeteta);
  AOD_PFJetPhi_.push_back(jetphi);
  AOD_PFJetID_.push_back(AOD_PFJetIDdecision);    
  
  //AlphaMax-type variables
  AOD_PFJetAlphaMax_       .push_back(alphaMax      ) ; 
  AOD_PFJetAlphaMax2_      .push_back(alphaMax2     ) ; 
  AOD_PFJetAlphaMaxPrime_  .push_back(alphaMaxPrime ) ; 
  AOD_PFJetAlphaMaxPrime2_ .push_back(alphaMaxPrime2) ; 
  AOD_PFJetBeta_           .push_back(beta          ) ; 
  AOD_PFJetBeta2_          .push_back(beta2         ) ; 

  //Totals
  AOD_PFJetSumIP_.push_back(sumIP);
  AOD_PFJetSumIPSig_.push_back(sumIPSig);

  AOD_PFJetTotalTrackAngle_.push_back(totalTrackAngle);    

  /////Medians
  AOD_PFJetMedianIP_             .push_back(medianIP);
  AOD_PFJetMedianLog10IPSig_     .push_back(log10(medianIPSig));
  AOD_PFJetMedianLog10TrackAngle_.push_back(log10(medianTrackAngle));

 }//end pf loop

 // AOD PFchs Jets -------------------------------------------
 for (edm::View<reco::PFJet>::const_iterator iJet = AODak4PFJetsCHSHandle->begin(); iJet != AODak4PFJetsCHSHandle->end(); ++iJet) {

  if(verbose_AOD) printf("PFchs %f \n",iJet->pt());
  
  float jetpt  = iJet->pt();
  float jeteta = iJet->eta();
  float jetphi = iJet->phi();

  // ID and jet selections
  bool pfchsjetPassLooseID = false;
  bool pfchsjetPassTightID = false;

  // from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016
  float pfchsjetNHF                 = iJet->neutralHadronEnergyFraction();
  float pfchsjetNEMF                = iJet->neutralEmEnergyFraction();
  float pfchsjetCHF                 = iJet->chargedHadronEnergyFraction();
  //float pfchsjetMUF                 = iJet->muonEnergyFraction();
  float pfchsjetCEMF                = iJet->chargedEmEnergyFraction();
  float pfchsjetNumConst            = iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
  float pfchsjetNumNeutralParticle  = iJet->neutralMultiplicity();
  float pfchsjetCHM                 = iJet->chargedMultiplicity(); 

  if ( fabs(jeteta) <= 2.7 ){
   pfchsjetPassLooseID = (pfchsjetNHF<0.99 && pfchsjetNEMF<0.99 && pfchsjetNumConst>1) && ((abs(jeteta)<=2.4 && pfchsjetCHF>0 && pfchsjetCHM>0 && pfchsjetCEMF<0.99) || abs(jeteta)>2.4) && abs(jeteta)<=2.7 ;
   pfchsjetPassTightID = (pfchsjetNHF<0.90 && pfchsjetNEMF<0.90 && pfchsjetNumConst>1) && ((abs(jeteta)<=2.4 && pfchsjetCHF>0 && pfchsjetCHM>0 && pfchsjetCEMF<0.99) || abs(jeteta)>2.4) && abs(jeteta)<=2.7 ;
  }
  if ( fabs(jeteta) > 2.7 && fabs(jeteta) <= 3.0 ){
   pfchsjetPassLooseID = (pfchsjetNEMF<0.90 && pfchsjetNumNeutralParticle>2 ) ;
   pfchsjetPassTightID = (pfchsjetNEMF<0.90 && pfchsjetNumNeutralParticle>2 ) ;
  }
  if ( fabs(jeteta) > 3.0 ){
   pfchsjetPassLooseID = (pfchsjetNEMF<0.90 && pfchsjetNumNeutralParticle>10 ) ;
   pfchsjetPassTightID = (pfchsjetNEMF<0.90 && pfchsjetNumNeutralParticle>10 ) ; 
  }

  Int_t AOD_PFchsJetIDdecision = 0;
  if (pfchsjetPassLooseID) AOD_PFchsJetIDdecision += pow(2, 1);
  if (pfchsjetPassTightID) AOD_PFchsJetIDdecision += pow(2, 2);

  if(iJet->pt()<20.0 || fabs(iJet->eta())>2.4 || AOD_PFchsJetIDdecision==0) continue;

  // pfchsJetTrackIDs is a vector of ints where each int is the 
  // index of a track passing deltaR requirement to this jet
  // out of the master track record of tracks passing basic selections
  vector<int>   pfchsJetTrackIDs = getJetTrackIndexs( jeteta, jetphi );
  if(pfchsJetTrackIDs.size()<1) continue;

  if(verbose_AOD){
   printf(" AOD Jet pt eta phi: %0.1f %0.1f %0.1f\n",jetpt,jeteta,jetphi);
   for( int i=0; i<(int)AOD_allTrackPositions.size(); i++){
    printf("  allTrack %i %0.1f %0.1f %0.1f \n",i,
     AOD_allTrackPt [i] ,
     AOD_allTrackEta[i] ,
     AOD_allTrackPhi[i] );

   }
   for( int i=0; i<(int)pfchsJetTrackIDs.size(); i++){
    printf(" Track %i at %i \n",i,pfchsJetTrackIDs[i]);
    printf("  pfchsTrack %i=%i %0.1f %0.1f %0.1f \n",i,
     pfchsJetTrackIDs[i],
     AOD_allTrackPt [pfchsJetTrackIDs[i]] ,
     AOD_allTrackEta[pfchsJetTrackIDs[i]] ,
     AOD_allTrackPhi[pfchsJetTrackIDs[i]] );
   }
  }

  // initialize variables
  float alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2 = -1.;
  float totalTrackAngle, totalTrackAnglePt = 0.;
  float sumIP, sumIPSig = 0.;
  vector<float> pfchsJetTrackAngles; 
  pfchsJetTrackAngles.clear();
  vector<float> pfchsJetIPs; 
  pfchsJetIPs.clear();
  vector<float> pfchsJetIPSigs; 
  pfchsJetIPSigs.clear();

  // do calculations
  calculateAlphaMax(pfchsJetTrackIDs,alphaMax,alphaMaxPrime,beta,alphaMax2,alphaMaxPrime2,beta2);
  calculateTrackAngle(pfchsJetTrackIDs, pfchsJetTrackAngles, totalTrackAngle, totalTrackAnglePt);
  calculateIP(pfchsJetTrackIDs, pfchsJetIPs, pfchsJetIPSigs, sumIP, sumIPSig);
  //calculateDisplacedVertices(es, pfchsJetTrackIDs);

  // find medians
  float medianTrackAngle;
  medianTrackAngle = findMedian(pfchsJetTrackAngles);
  float medianIP;
  medianIP = findMedian(pfchsJetIPs);
  float medianIPSig;
  medianIPSig = findMedian(pfchsJetIPSigs);

  ////////////////////////
  // Fill tree
  /////////////////////////
  AOD_nPFchsJet_++;
  
  //Pt, Eta, Phi
  AOD_PFchsJetPt_.push_back(jetpt);
  AOD_PFchsJetEta_.push_back(jeteta);
  AOD_PFchsJetPhi_.push_back(jetphi);
  AOD_PFchsJetID_.push_back(AOD_PFchsJetIDdecision);
  
  //AlphaMax-type variables
  AOD_PFchsJetAlphaMax_       .push_back(alphaMax      ) ; 
  AOD_PFchsJetAlphaMax2_      .push_back(alphaMax2     ) ; 
  AOD_PFchsJetAlphaMaxPrime_  .push_back(alphaMaxPrime ) ; 
  AOD_PFchsJetAlphaMaxPrime2_ .push_back(alphaMaxPrime2) ; 
  AOD_PFchsJetBeta_           .push_back(beta          ) ; 
  AOD_PFchsJetBeta2_          .push_back(beta2         ) ; 

  //Totals
  AOD_PFchsJetSumIP_.push_back(sumIP);
  AOD_PFchsJetSumIPSig_.push_back(sumIPSig);

  AOD_PFchsJetTotalTrackAngle_.push_back(totalTrackAngle);    

  /////Medians
  AOD_PFchsJetMedianIP_             .push_back(medianIP);
  AOD_PFchsJetMedianLog10IPSig_     .push_back(log10(medianIPSig));
  AOD_PFchsJetMedianLog10TrackAngle_.push_back(log10(medianTrackAngle));

 }//end pfchs loop
 
}//end fill jets


vector<int> lldjNtuple::getJetTrackIndexs( float jeteta, float jetphi )
{
 vector<int> idvector;
 // loop over all selected tracks, dR match to jet
 for( int i=0; i<(int)AOD_allTrackPositions.size(); i++){
  float tracketa = AOD_allTrackPositions[i].Eta(); 
  float trackphi = AOD_allTrackPositions[i].Phi(); 
  float drt = deltaR( jeteta, jetphi, tracketa, trackphi );
  if(drt < maxDRtrackJet_){ idvector.push_back(i); }
 }
 return idvector;
}


void lldjNtuple::calculateAlphaMax(vector<int> jetTrackIDs, float& aMax, float& aMaxP, float& beta, float& aMax2, float& aMaxP2, float& beta2)
{

  float trackSumPt = 0; 
  float trackSumPt2 = 0; 
  float trackSumPtVtxMatched = 0; 
  float trackSumPtVtxMatched2 = 0; 

  int nrvtxs = AODVertexHandle->size();
  vector<float> trackSumPtByVtx( nrvtxs,0);
  vector<float> trackSumPtByVtx2(nrvtxs,0);

  // printf(" jetTracksIDs size: %i \n",(int)jetTrackIDs.size() );
  // printf(" AOD_whichVertexByTrack size: %i \n",(int)AOD_whichVertexByTrack.size() );

  for(int t=0; t< (int)jetTrackIDs.size(); t++){
   int trackID = jetTrackIDs[t];
  
   // sum pt of all tracks passing dR cut
   float trackpt = AOD_allTrackIFSPt[trackID];
   trackSumPt += trackpt;
   trackSumPt2 += (trackpt*trackpt);
  
   // now only for tracks associated with a vertex
   // the index of best vertex for track j is AOD_whichVertexByTrack[j]
   if(AOD_whichVertexByTrack[trackID] < 0)continue;

   // technically we could get this by summing over trackSumPtByVtx later
   trackSumPtVtxMatched += trackpt;
   trackSumPtVtxMatched2 += trackpt*trackpt;

   //// trackSumPtByVtx are sorted by vertex 
   //printf("  track %i to vtx %i TS %i\n",
   // trackID, AOD_whichVertexByTrack[trackID],
   // (int)trackSumPtByVtx.size());
   trackSumPtByVtx[AOD_whichVertexByTrack[trackID]] += trackpt;
   trackSumPtByVtx2[AOD_whichVertexByTrack[trackID]] += (trackpt*trackpt);
  }
  
  // clear variables from previous execution
  double apMax =0;
  double apMax2 = 0;
  
  // calculate beta
  // beta = ratio of sum of track pt matched to any vertex / sum on all tracks
  beta = 1.0 - trackSumPtVtxMatched/trackSumPt;
  beta2 = 1.0 - trackSumPtVtxMatched2 / trackSumPt2;
  
  // calculate alpha
  // alpha[a] = sum pt tracks matched to vtx [a] / total track sum pt
  // ap[a] = sum pt tracks matched to vtx [a] / ( same sum + beta ) 
  float tmpMaxSumPt = -1.;
  for(int i = 0; i < nrvtxs; i++){
    // find vertex number of max pt sum
    if(trackSumPtByVtx[i] > tmpMaxSumPt){
     tmpMaxSumPt = trackSumPtByVtx[i];
    }
    // calculate and fill apMax 
    double ap = trackSumPtByVtx[i] / (trackSumPtByVtx[i] + beta);
    double ap2 = trackSumPtByVtx2[i] / (trackSumPtByVtx2[i] + beta2);
    if(ap > apMax) apMax = ap;
    if(ap2 > apMax2) apMax2 = ap2;
  }
  aMax   = tmpMaxSumPt  / trackSumPt;
  aMax2  = ( tmpMaxSumPt * tmpMaxSumPt ) / trackSumPt2;
  aMaxP  = apMax;
  aMaxP2 = apMax2;
  return;

}

void lldjNtuple::calculateTrackAngle(vector<int> jetTrackIDs, vector<float> &allTrackAngles,
 float &totalTrackAngle, float &totalTrackAnglePt)
{
  
  for(int t=0; t< (int)jetTrackIDs.size(); t++){
    int trackID = jetTrackIDs[t];
    // sum pt of all tracks passing dR cut
    float trackpt    = AOD_allTrackIFSPt[trackID];
    float trackangle = AOD_allTrackAngle[trackID];

    totalTrackAngle   += trackangle;
    totalTrackAnglePt += trackangle*trackpt;

    allTrackAngles.push_back(trackangle);
  }

  return;
  
}

void lldjNtuple::calculateIP(vector<int> jetTrackIDs, vector<float> &jetIPs, vector<float> &jetIPSigs, float &tsumIP, float &tsumIPSig)
{
  
  for(int t=0; t< (int)jetTrackIDs.size(); t++){

    int trackID = jetTrackIDs[t];
    // sum pt of all tracks passing dR cut
    float trackdxy    = AOD_allTrackdxy   [trackID];
    float trackdxyerr = AOD_allTrackdxyerr[trackID];
    float trackIPSig  = 0;
    if(trackdxyerr>0.) trackIPSig = trackdxy/trackdxyerr;

    tsumIP    += trackdxy;
    tsumIPSig += trackIPSig;

    //printf(" aa trackdyx: %0.4f %0.4f \n", trackdxy, trackIPSig);

    jetIPs.push_back(trackdxy);
    jetIPSigs.push_back(trackIPSig);

  }

  return;
  
}


float lldjNtuple::findMedian( vector<float> thevector){

 float themedian = -999. ;

 //printf(" thevector: ");
 //for ( int i=0; i<(int)thevector.size(); i++ ){
 // printf(" %0.4f ",thevector.at(i));
 //}
 //printf("\n");
 std::sort (thevector.begin(), thevector.end());
 //printf(" sorted thevector: ");
 //for ( int i=0; i<(int)thevector.size(); i++ ){
 // printf(" %0.4f ",thevector.at(i));
 //}
 //printf("\n");

 if(thevector.size() > 0){
  if((thevector.size()%2 == 0)){
   themedian = (
    thevector.at( (thevector.size() / 2) - 1 )
  + thevector.at( (thevector.size() / 2)     ) 
   ) / 2 ;
  } else {
   themedian = thevector.at( ( thevector.size() - 1 ) / 2 ) ;
  }
 }

 return themedian;

}

void lldjNtuple::deltaVertex3D(GlobalPoint secVert, std::vector<reco::TransientTrack> tracks, double& dEta, double& dPhi, double& pt, double& m, double& energy)
{
  TVector3 pv(AODVertexHandle->at(0).x(),AODVertexHandle->at(0).y(),AODVertexHandle->at(0).z());
  TVector3 sv(secVert.x(),secVert.y(),secVert.z());
  TVector3 diff = (sv-pv);
  TVector3 trackPt(0,0,0);
  TLorentzVector trackP4(0,0,0,0);
  for(int i = 0; i < (int)tracks.size(); i++){
    TVector3 tt;
    tt.SetPtEtaPhi(tracks[i].trajectoryStateClosestToPoint(secVert).momentum().transverse(),tracks[i].trajectoryStateClosestToPoint(secVert).momentum().eta(),tracks[i].trajectoryStateClosestToPoint(secVert).momentum().phi());
    trackPt += tt;
    trackP4 += TLorentzVector(tt,tracks[i].trajectoryStateClosestToPoint(secVert).momentum().mag());
  }
  dEta = diff.Eta()-trackPt.Eta();
  dPhi = diff.DeltaPhi(trackPt);
  pt = (trackPt - ((trackPt * diff)/(diff * diff) * diff)).Mag();
  m = trackP4.M();
  energy = trackP4.E();
}

void lldjNtuple::deltaVertex2D(GlobalPoint secVert, std::vector<reco::TransientTrack> tracks, double& dPhi, double& pt, double& mediandPhi)
{

  //edm::Handle<reco::BeamSpot> beamspotHandle_;//make this global??? BENTODO
  //e.getByToken(beamspotLabel_, beamspotHandle_);

  const reco::BeamSpot& pat_beamspot = (*beamspotHandle_);
  TVector2 bmspot(pat_beamspot.x0(),pat_beamspot.y0());
  TVector2 sv(secVert.x(),secVert.y());
  TVector2 diff = (sv-bmspot);
  TVector2 trackPt(0,0);
  vector<double> trackAngles;
  for(int i = 0; i < (int)tracks.size(); i++){
    TVector2 tt;
    tt.SetMagPhi(tracks[i].trajectoryStateClosestToPoint(secVert).momentum().transverse(),tracks[i].trajectoryStateClosestToPoint(secVert).momentum().phi());
    trackPt += tt;
    trackAngles.push_back(fabs(diff.DeltaPhi(tt)));
  }
  sort(trackAngles.begin(), trackAngles.end());

  if(trackAngles.size() == 0){
    //do nothing
  }else if((trackAngles.size()%2 == 0)){
    mediandPhi = trackAngles.at(trackAngles.size()/2-1);
  }else{
    mediandPhi = trackAngles.at((trackAngles.size()-1)/2);
  }

  dPhi = diff.DeltaPhi(trackPt);
  pt = (trackPt - ((trackPt * diff)/(diff * diff) * diff)).Mod();

}


vector<reco::TransientTrack> lldjNtuple::cleanTracks(vector<reco::TransientTrack> tracks, GlobalPoint vertPos)
{
  vector<reco::TransientTrack> cleanTracks;
  for(int i = 0; i < (int)tracks.size(); i++){
    if(tracks[i].trajectoryStateClosestToPoint(vertPos).perigeeError().transverseImpactParameterError() > 0 && tracks[i].trajectoryStateClosestToPoint(vertPos).perigeeParameters().transverseImpactParameter() / tracks[i].trajectoryStateClosestToPoint(vertPos).perigeeError().transverseImpactParameterError() > 3.0)continue;
    cleanTracks.push_back(tracks[i]);
  }
  return cleanTracks;
}



void lldjNtuple::calculateDisplacedVertices(const edm::EventSetup& es, vector<int> jetTrackIDs){
  
  //Select transient tracks for this jet
  vector<reco::TransientTrack> transientTracks;
  for(unsigned int j = 0; j < jetTrackIDs.size(); j++){
    reco::TransientTrack tt(AODTrackHandle->at( jetTrackIDs.at(j) ),magneticField_);
    transientTracks.push_back(tt);
  }
  
  if(transientTracks.size() > 1){
    vector<TransientVertex> avfVerts = vtxfitter_->vertices(transientTracks);
    if(avfVerts.size() > 0 && avfVerts[0].isValid()){
      GlobalPoint vertPos = avfVerts[0].position();
      GlobalError vertErr = avfVerts[0].positionError();
      
      AOD_CaloJetAvfVx_.push_back( vertPos.x() );
      AOD_CaloJetAvfVy_.push_back( vertPos.y() );
      AOD_CaloJetAvfVz_.push_back( vertPos.z() );
      
      float vertBeamXY = vertexBeam_->distanceToBeam(reco::Vertex(RecoVertex::convertPos(vertPos),RecoVertex::convertError(vertErr)));
      
      AOD_CaloJetAvfVertexTotalChiSquared_.push_back( avfVerts[0].totalChiSquared() );
      AOD_CaloJetAvfVertexDegreesOfFreedom_.push_back( avfVerts[0].degreesOfFreedom() );
      if(avfVerts[0].degreesOfFreedom() > 0) AOD_CaloJetAvfVertexChi2NDoF_.push_back( avfVerts[0].totalChiSquared()/avfVerts[0].degreesOfFreedom() );
       else AOD_CaloJetAvfVertexChi2NDoF_.push_back( 0. );
      AOD_CaloJetAvfVertexDistanceToBeam_.push_back( vertBeamXY );
      double rerr = vertErr.rerr(vertPos);
      AOD_CaloJetAvfVertexTransverseError_.push_back( rerr );
      if(rerr > 0) AOD_CaloJetAvfVertexTransverseSig_.push_back( vertBeamXY/rerr );
       else AOD_CaloJetAvfVertexTransverseSig_.push_back( 0. );
      
      vector<reco::TransientTrack> cleanTrackColl = cleanTracks(avfVerts[0].refittedTracks(),vertPos);
      
      double d3deta = 0, d3dphi = 0, d3dpt = 0, d3m = 0, d3en;
      deltaVertex3D(vertPos, cleanTrackColl,d3deta,d3dphi,d3dpt,d3m,d3en);
      AOD_CaloJetAvfVertexDeltaEta_.push_back( d3deta );
      AOD_CaloJetAvfVertexDeltaPhi_.push_back( d3dphi );
      AOD_CaloJetAvfVertexRecoilPt_.push_back( d3dpt );
      AOD_CaloJetAvfVertexTrackMass_.push_back( d3m );
      AOD_CaloJetAvfVertexTrackEnergy_.push_back( d3en );
      double d2dphi = 0,d2dpt = 0,d2dphiMed=1e-6;
      deltaVertex2D(vertPos,cleanTrackColl,d2dphi,d2dpt,d2dphiMed);
      AOD_CaloJetAvfBeamSpotDeltaPhi_.push_back( d2dphi );
      AOD_CaloJetAvfBeamSpotRecoilPt_.push_back( d2dpt );
      AOD_CaloJetAvfBeamSpotMedianDeltaPhi_.push_back( d2dphiMed );
      AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_.push_back( log10(d2dphiMed) );
      AOD_CaloJetNCleanMatchedTracks_.push_back( (int)cleanTrackColl.size() );
      
      int nHitsInFront = 0;
      int nMissingAfter = 0;
      CheckHitPattern checkHitPattern_;
      for(int j = 0; j < (int)cleanTrackColl.size(); j++){
	CheckHitPattern::Result res = checkHitPattern_.analyze(es,cleanTrackColl[j].track(),avfVerts[0].vertexState(),false);
	nHitsInFront += res.hitsInFrontOfVert;
	nMissingAfter += res.missHitsAfterVert;
      }
      AOD_CaloJetSumHitsInFrontOfVert_.push_back( nHitsInFront );
      AOD_CaloJetSumMissHitsAfterVert_.push_back( nMissingAfter );
      AOD_CaloJetHitsInFrontOfVertPerTrack_.push_back( double(nHitsInFront)/double(transientTracks.size()) );
      AOD_CaloJetMissHitsAfterVertPerTrack_.push_back( double(nMissingAfter)/double(transientTracks.size()) );
      
      AOD_CaloJetAvfDistToPV_.push_back( 
				       sqrt(pow((AODVertexHandle->at(0).x() - avfVerts[0].position().x())/AODVertexHandle->at(0).x(),2)
					    +pow((AODVertexHandle->at(0).y() - avfVerts[0].position().y())/AODVertexHandle->at(0).y(),2)
					    +pow((AODVertexHandle->at(0).z() - avfVerts[0].position().z())/AODVertexHandle->at(0).z(),2)) );
      if(AODVertexHandle->size() > 0)AOD_CaloJetAvfVertexDeltaZtoPV_.push_back( AODVertexHandle->at(0).z() - avfVerts[0].position().z() );
       else AOD_CaloJetAvfVertexDeltaZtoPV_.push_back( 0. );
      if(AODVertexHandle->size() > 1)AOD_CaloJetAvfVertexDeltaZtoPV2_.push_back( AODVertexHandle->at(1).z() - avfVerts[0].position().z() );
       else AOD_CaloJetAvfVertexDeltaZtoPV2_.push_back( 0. );
      
      
    }//end valid avf vertex
    else{
     AOD_CaloJetAvfVx_.push_back( 0. );
     AOD_CaloJetAvfVy_.push_back( 0. );
     AOD_CaloJetAvfVz_.push_back( 0. );
     AOD_CaloJetAvfVertexTotalChiSquared_.push_back( 0. );
     AOD_CaloJetAvfVertexDegreesOfFreedom_.push_back( 0. );
     AOD_CaloJetAvfVertexChi2NDoF_.push_back( 0. );
     AOD_CaloJetAvfVertexDistanceToBeam_.push_back( 0. );
     AOD_CaloJetAvfVertexTransverseError_.push_back( 0. );
     AOD_CaloJetAvfVertexTransverseSig_.push_back( 0. );
     AOD_CaloJetAvfVertexDeltaEta_.push_back( 0. );
     AOD_CaloJetAvfVertexDeltaPhi_.push_back( 0. );
     AOD_CaloJetAvfVertexRecoilPt_.push_back( 0. );
     AOD_CaloJetAvfVertexTrackMass_.push_back( 0. );
     AOD_CaloJetAvfVertexTrackEnergy_.push_back( 0. );
     AOD_CaloJetAvfBeamSpotDeltaPhi_.push_back( 0. );
     AOD_CaloJetAvfBeamSpotRecoilPt_.push_back( 0. );
     AOD_CaloJetAvfBeamSpotMedianDeltaPhi_.push_back( 0. );
     AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_.push_back( 0. );
     AOD_CaloJetNCleanMatchedTracks_.push_back( 0. );
     AOD_CaloJetSumHitsInFrontOfVert_.push_back( 0. );
     AOD_CaloJetSumMissHitsAfterVert_.push_back( 0. );
     AOD_CaloJetHitsInFrontOfVertPerTrack_.push_back( 0. );
     AOD_CaloJetMissHitsAfterVertPerTrack_.push_back( 0. );
     AOD_CaloJetAvfDistToPV_.push_back( 0. );
     AOD_CaloJetAvfVertexDeltaZtoPV_.push_back( 0. );
     AOD_CaloJetAvfVertexDeltaZtoPV2_.push_back( 0. );
    }

  }//end if transientTracks
  else{
   AOD_CaloJetAvfVx_.push_back( 0. );
   AOD_CaloJetAvfVy_.push_back( 0. );
   AOD_CaloJetAvfVz_.push_back( 0. );
   AOD_CaloJetAvfVertexTotalChiSquared_.push_back( 0. );
   AOD_CaloJetAvfVertexDegreesOfFreedom_.push_back( 0. );
   AOD_CaloJetAvfVertexChi2NDoF_.push_back( 0. );
   AOD_CaloJetAvfVertexDistanceToBeam_.push_back( 0. );
   AOD_CaloJetAvfVertexTransverseError_.push_back( 0. );
   AOD_CaloJetAvfVertexTransverseSig_.push_back( 0. );
   AOD_CaloJetAvfVertexDeltaEta_.push_back( 0. );
   AOD_CaloJetAvfVertexDeltaPhi_.push_back( 0. );
   AOD_CaloJetAvfVertexRecoilPt_.push_back( 0. );
   AOD_CaloJetAvfVertexTrackMass_.push_back( 0. );
   AOD_CaloJetAvfVertexTrackEnergy_.push_back( 0. );
   AOD_CaloJetAvfBeamSpotDeltaPhi_.push_back( 0. );
   AOD_CaloJetAvfBeamSpotRecoilPt_.push_back( 0. );
   AOD_CaloJetAvfBeamSpotMedianDeltaPhi_.push_back( 0. );
   AOD_CaloJetAvfBeamSpotLog10MedianDeltaPhi_.push_back( 0. );
   AOD_CaloJetNCleanMatchedTracks_.push_back( 0. );
   AOD_CaloJetSumHitsInFrontOfVert_.push_back( 0. );
   AOD_CaloJetSumMissHitsAfterVert_.push_back( 0. );
   AOD_CaloJetHitsInFrontOfVertPerTrack_.push_back( 0. );
   AOD_CaloJetMissHitsAfterVertPerTrack_.push_back( 0. );
   AOD_CaloJetAvfDistToPV_.push_back( 0. );
   AOD_CaloJetAvfVertexDeltaZtoPV_.push_back( 0. );
   AOD_CaloJetAvfVertexDeltaZtoPV2_.push_back( 0. );
  }

}//end calculateDisplacedVertices
