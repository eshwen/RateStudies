import FWCore.ParameterSet.Config as cms
import sys
process = cms.Process("Demo")
outdir = '/data_CMS/cms/amendola/RateStudiesL1Ntuples/Ntuples22Gen2016/'
readme = open(outdir+"readme.txt","w")
sourceFile = 'root://polgrid4.in2p3.fr//store/data/Run2016D/ZeroBias/MINIAOD/23Sep2016-v1/100000/047F4BED-BD84-E611-9045-44A842CFD5FF.root'
readme.write('Source '+ sourceFile)

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
process.source = cms.Source("PoolSource", 
#    fileNames = cms.untracked.vstring("root://polgrid4.in2p3.fr//store/data/Run2016D/ZeroBias/MINIAOD/23Sep2016-v1/100000/047F4BED-BD84-E611-9045-44A842CFD5FF.root")
                            fileNames = cms.untracked.vstring(sourceFile)
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
                                   fileName = cms.string(outdir+'Data_L1ntuple_ZeroBias.root'),
                                   closeFileFast = cms.untracked.bool(True)
                                   )



process.p = cms.Path(process.L1Tree)
#process.ep = cms.EndPath(process.TFileService)
