import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process("Demo")
sourceFile = 'root://polgrid4.in2p3.fr//store/data/Run2016H/ZeroBiasBunchTrains0/AOD/09Nov2016-v1/110000/005EAF54-21A8-E611-9D30-0CC47A4D7640.root'
outdir = '/data_CMS/cms/amendola/RateStudiesL1Ntuples/Ntuples23Gen2016/'
outFile = 'Data_L1ntuple_ZeroBias_test1000entries.root'
MAX = 1000

readme = open(outdir+"readme.txt","w")
readme.write('Source '+ sourceFile+'\nEvents {0}'.format(MAX))

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAX) )
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(sourceFile)
)


process.L1Tree = cms.EDAnalyzer("L1ntuples",
                                stage2TauCollection = cms.InputTag("caloStage2Digis","Tau"),
                                stage2EGammaCollection = cms.InputTag("caloStage2Digis","EGamma"),
                                stage2MuonCollection = cms.InputTag("gmtStage2Digis","Muon"),
                                stage2JetCollection = cms.InputTag("caloStage2Digis","Jet")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outdir+outFile),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



process.p = cms.Path(process.L1Tree)

