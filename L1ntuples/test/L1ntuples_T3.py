import FWCore.ParameterSet.Config as cms
import sys,os

PyFilePath = os.environ['CMSSW_BASE']+"/src/RateStudies/L1ntuples/"
process = cms.Process("Demo")


process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(XXX_MAXEVENTS_XXX) )

execfile(PyFilePath+"test/XXX_SAMPLEFILENAME_XXX")
process.source = cms.Source("PoolSource", 

                            fileNames = FILELIST
)


#process.demo = cms.EDAnalyzer('L1ntuples'
#     , tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
#)
process.L1Tree = cms.EDAnalyzer("L1ntuples",
                                stage2TauCollection = cms.InputTag("caloStage2Digis","Tau"),
                                stage2EGammaCollection = cms.InputTag("caloStage2Digis","EGamma"),
                                stage2MuonCollection = cms.InputTag("gmtStage2Digis","Muon"),
                                stage2JetCollection = cms.InputTag("caloStage2Digis","Jet")
                                )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("XXX_OUTPUTFILENTUPLE_XXX"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



process.p = cms.Path(process.L1Tree)
#process.ep = cms.EndPath(process.TFileService)
