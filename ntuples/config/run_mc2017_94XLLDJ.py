import FWCore.ParameterSet.Config as cms

##########################################################################################
# Setup

# this is the process run by cmsRun
process = cms.Process('LLDJ')
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
process.load("FWCore.MessageLogger.MessageLogger_cfi")

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

process.load("RecoTracker.TkNavigation.NavigationSchoolESProducer_cfi")


process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v12')

#process.Tracer = cms.Service("Tracer")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'file:/uscms_data/d3/tmperry/roots/miniAOD2017_TTJets_C2D182B7-D454-E811-93EE-3417EBE70663.root'
        'file:/uscms_data/d3/tmperry/aaLLDJ_slc6_630_CMSSW_9_4_9_cand2/src/VHLLDJanalysis2017/ntuples/config/TT_reMiniAOD.root',
        #'file:/uscms_data/d3/tmperry/ggAnalysis_slc6_630_CMSSW_9_4_9_cand2/src/ggAnalysis/ggNtuplizer/test/TT_reMiniAOD_1000.root'
        ))

#process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=False, 
                       era='2017-Nov17ReReco')  

process.TFileService = cms.Service("TFileService", fileName = cms.string('LLDJntuple_mc.root'))

# MET correction and uncertainties
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=False
                           )


# For AOD Track variables
process.MaterialPropagator = cms.ESProducer('PropagatorWithMaterialESProducer',
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)

process.TransientTrackBuilderESProducer = cms.ESProducer('TransientTrackBuilderESProducer',
    ComponentName = cms.string('TransientTrackBuilder')
)

#NTuplizer
process.lldjNtuple = cms.EDAnalyzer('lldjNtuple',

 doLLDJAOD                 = cms.bool(True),


 #trgFilterDeltaPtCut  = cms.double(0.5),
 #trgFilterDeltaRCut   = cms.double(0.3),

 miniAODtriggerEvent         = cms.InputTag("slimmedPatTrigger", "", ""),
 miniAODtriggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
 miniAODpatTriggerResults    = cms.InputTag("TriggerResults", "", "PAT"),
 #patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
 miniAODgenParticleSrc       = cms.InputTag("prunedGenParticles"),
 miniAODgeneratorLabel       = cms.InputTag("generator"),
 miniAODLHEEventLabel        = cms.InputTag("externalLHEProducer"),
 miniAODnewParticles         = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
 miniAODpileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
 miniAODVtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
 miniAODVtxBSLabel           = cms.InputTag("offlinePrimaryVerticesWithBS"),
 miniAODrhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
 miniAODrhoCentralLabel      = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
 miniAODpfMETLabel           = cms.InputTag("slimmedMETs"),
 miniAODelectronSrc          = cms.InputTag("slimmedElectrons"),
 #calibelectronSrc     = cms.InputTag("calibratedPatElectrons"),
 miniAODcalibelectronSrc     = cms.InputTag("slimmedElectrons"),
 miniAODphotonSrc            = cms.InputTag("slimmedPhotons"),
 #calibphotonSrc       = cms.InputTag("calibratedPatPhotons"),
 miniAODcalibphotonSrc       = cms.InputTag("slimmedPhotons"),
 miniAODmuonSrc              = cms.InputTag("slimmedMuons"),
 miniAODgsfTrackSrc          = cms.InputTag("reducedEgamma", "reducedGsfTracks"),
 miniAODebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
 miniAODeeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
 miniAODesReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits"),
 miniAODrecoPhotonSrc             = cms.InputTag("reducedEgamma", "reducedGedPhotonCores"),
 miniAODTrackLabel                = cms.InputTag("generalTracks"),
 miniAODgsfElectronLabel          = cms.InputTag("gsfElectrons"),
 miniAODPFAllCandidates           = cms.InputTag("particleFlow"),
 #ak4JetSrc                 = cms.InputTag("updatedJets"),
 miniAODak4JetSrc                 = cms.InputTag("slimmedJets"),
 #ak4JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
 miniAODak8JetSrc                 = cms.InputTag("slimmedJetsAK8"),
 #ak8JetSrc                 = cms.InputTag("selectedUpdatedPatJetsUpdatedJECAK8"),
 #boostedDoubleSVLabel      = cms.InputTag("pfBoostedDoubleSecondaryVertexAK8BJetTags"),
 miniAODtauSrc                    = cms.InputTag("slimmedTaus"),
 #pfLooseId                 = pfJetIDSelector.clone(),

 miniAODpackedPFCands             = cms.InputTag("packedPFCandidates"),
 miniAODelePFClusEcalIsoProducer  = cms.InputTag("electronEcalPFClusterIsolationProducer"),
 miniAODelePFClusHcalIsoProducer  = cms.InputTag("electronHcalPFClusterIsolationProducer"),
 miniAODBadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
 miniAODBadPFMuonFilter           = cms.InputTag("BadPFMuonFilter")




#
#
#
# miniAODelectronSrc        = cms.InputTag('slimmedElectrons'),
#
# doAOD                     = cms.bool(False),
# doMiniAOD                 = cms.bool(False),
#
# electronSrc               = cms.InputTag('selectedElectrons','','LLDJ'),
# rhoLabel                  = cms.InputTag('fixedGridRhoFastjetAll'),
# eleVetoIdMap              = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto'),
# eleLooseIdMap             = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose'),
# eleMediumIdMap            = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium'),
# eleTightIdMap             = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight'),
# eleHLTIdMap               = cms.InputTag('egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1'),
# #eleMVAValuesMap           = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values'),
# #eleMVAHZZValuesMap        = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values'),
# #elePFClusEcalIsoProducer  = cms.InputTag('electronEcalPFClusterIsolationProducer'),
# #elePFClusHcalIsoProducer  = cms.InputTag('electronHcalPFClusterIsolationProducer'),
#
# rhoCentralLabel           = cms.InputTag('fixedGridRhoFastjetCentralNeutral'),
# pileupCollection          = cms.InputTag('slimmedAddPileupInfo'),
# AODpileupCollection       = cms.InputTag('addPileupInfo', '', 'HLT'),
# VtxLabel                  = cms.InputTag('offlineSlimmedPrimaryVertices'),
# triggerResults            = cms.InputTag('TriggerResults', '', 'HLT'),
#
# AODTriggerInputTag           = cms.InputTag("TriggerResults","","HLT"),
# AODTriggerEventInputTag      = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
#
# beamspotLabel_            = cms.InputTag('offlineBeamSpot'),
#
# ak4JetSrc                 = cms.InputTag('slimmedJets'),
# AODak4CaloJetsSrc         = cms.InputTag('ak4CaloJets' , '', 'RECO'),
# AODak4PFJetsSrc           = cms.InputTag('ak4PFJets'   , '', 'RECO'),
# AODak4PFJetsCHSSrc        = cms.InputTag('ak4PFJetsCHS', '', 'RECO'),
# selectedPatJetsSrc        = cms.InputTag('selectedPatJets'),                                   
# AODVertexSrc              = cms.InputTag('offlinePrimaryVertices', '', 'RECO'),
# AODTrackSrc               = cms.InputTag('generalTracks', '', 'RECO'),
# vertexFitterConfig = cms.PSet(
#        finder = cms.string('avf'),
#        sigmacut = cms.double(10.),
#        Tini = cms.double(256.),
#        ratio = cms.double(0.25),
#        ),
#
# patTriggerResults         = cms.InputTag('TriggerResults', '', 'PAT'),
# BadChargedCandidateFilter = cms.InputTag('BadChargedCandidateFilter'),
# BadPFMuonFilter           = cms.InputTag('BadPFMuonFilter'),
# pfMETLabel                = cms.InputTag('slimmedMETs'),
# AODCaloMETlabel           = cms.InputTag('caloMet','','RECO'),    
# AODpfChMETlabel           = cms.InputTag('pfChMet','','RECO'),    
# AODpfMETlabel             = cms.InputTag('pfMet','','RECO'),  
#
# muonSrc                   = cms.InputTag('slimmedMuons'),
# muonAODSrc                = cms.InputTag('selectedPatMuons'),
#
# photonSrc                 = cms.InputTag('selectedPhotons','','LLDJ'),
# phoLooseIdMap             = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose'),
# phoMediumIdMap            = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium'),
# phoTightIdMap             = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight'),
# phoChargedIsolation       = cms.InputTag('photonIDValueMapProducer:phoChargedIsolation'),
# phoNeutralHadronIsolation = cms.InputTag('photonIDValueMapProducer:phoNeutralHadronIsolation'),
# phoPhotonIsolation        = cms.InputTag('photonIDValueMapProducer:phoPhotonIsolation'),
# phoWorstChargedIsolation  = cms.InputTag('photonIDValueMapProducer:phoWorstChargedIsolation'),
#
# #photonAODSrc              = cms.InputTag('selectedPatPhotons'),
# photonAODSrc              = cms.InputTag('gedPhotons'),
# AOD_phoLooseIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
# AOD_phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
# AOD_phoTightIdMap  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
# AOD_phoChargedIsolationMap       = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation"),
# AOD_phoNeutralHadronIsolationMap = cms.InputTag("photonIDValueMapProducer", "phoNeutralHadronIsolation"),
# AOD_phoPhotonIsolationMap        = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
# AOD_phoWorstChargedIsolationMap  = cms.InputTag("photonIDValueMapProducer", "phoWorstChargedIsolation"),
#
# electronAODSrc = cms.InputTag("gedGsfElectrons"),
# #AOD_eleIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),#doesn't work with AOD
# AOD_eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
# AOD_eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
# AOD_eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
# conversions  = cms.InputTag('allConversions'),                                    
#
# genParticleSrc    = cms.InputTag("genParticles"),
#
# bits = cms.InputTag("TriggerResults","","HLT"),
# prescales = cms.InputTag("patTrigger"),
# objects = cms.InputTag("selectedPatTrigger"),

)


#builds Ntuple
process.p = cms.Path(
    #process.egmPhotonIDSequence *
    process.lldjNtuple
    )


# process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")
# process.ggNtuplizer.year=cms.int32(2017)
# process.ggNtuplizer.doGenParticles=cms.bool(True)
# process.ggNtuplizer.dumpPFPhotons=cms.bool(False)
# process.ggNtuplizer.dumpHFElectrons=cms.bool(False)
# process.ggNtuplizer.dumpJets=cms.bool(True)
# process.ggNtuplizer.dumpAK8Jets=cms.bool(False)
# process.ggNtuplizer.dumpSoftDrop= cms.bool(True)
# process.ggNtuplizer.dumpTaus=cms.bool(False)
# process.ggNtuplizer.triggerEvent=cms.InputTag("slimmedPatTrigger", "", "PAT")
# 
# process.p = cms.Path(
#     process.ggNtuplizer
#     )

#print process.dumpPython()
