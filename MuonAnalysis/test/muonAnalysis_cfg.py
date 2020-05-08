import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'muonAnalyzer.root'
options.inputFiles = 'file:step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT.root'
options.maxEvents = -1
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Input file list")

options.parseArguments()

if options.inputFileList:
    with open(options.inputFileList) as f:
        inputFiles = f.readlines()
    options.inputFiles[:] = []
    options.inputFiles = [f if f.startswith('/store') else 'file:'+f for f in inputFiles]

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process = cms.Process("MuonAnalysis", Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

# output definition
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputFile),
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
)

process.schedule = cms.Schedule()

process.muonAnalyzer = cms.EDAnalyzer('MuonAnalyzer',
    muonsSrc = cms.InputTag('muons'),
    genSrc = cms.InputTag('genParticles'),
    verticesSrc = cms.InputTag('offlinePrimaryVertices4D'),
    pfSrc = cms.InputTag('particleFlow'),
    MuonSimInfoConfiguration = cms.VPSet(
        cms.PSet(
            muonSimInfo = cms.InputTag("muonSimClassifier"),
            muonSimInfoGenPrimary = cms.InputTag("muonSimClassifier","toPrimaries"),
            label = cms.string('global'),
        ),
        cms.PSet(
            muonSimInfo = cms.InputTag("muonSimClassifierInner"),
            muonSimInfoGenPrimary = cms.InputTag("muonSimClassifierInner","toPrimaries"),
            label = cms.string('inner'),
        ),
        cms.PSet(
            muonSimInfo = cms.InputTag("muonSimClassifierOuter"),
            muonSimInfoGenPrimary = cms.InputTag("muonSimClassifierOuter","toPrimaries"),
            label = cms.string('outer'),
        ),
    ),
)
process.muonAnalyzerPath = cms.Path()
process.muonAnalyzerPath += process.muonAnalyzer
process.schedule.append(process.muonAnalyzerPath)
