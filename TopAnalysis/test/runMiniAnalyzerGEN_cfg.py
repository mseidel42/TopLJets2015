import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('outFilename', 'MiniEvents.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('inputFile', None,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "input file to process"
                 )
options.register('saveTree', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "save summary tree"
                 )
options.register('savePF', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'save PF candidates'
                 )
options.register('scale', 1., VarParsing.multiplicity.singleton, VarParsing.varType.float, "factor for fsr ren scale")
options.register('asfsr', 0.1365, VarParsing.multiplicity.singleton, VarParsing.varType.float, "alpha_s for fsr")
options.register('me', 'on', VarParsing.multiplicity.singleton, VarParsing.varType.string, "ME corrections")
options.register('generator', 'pythia8', VarParsing.multiplicity.singleton, VarParsing.varType.string, "PS generator")
options.register('cr', 'default', VarParsing.multiplicity.singleton, VarParsing.varType.string, "color reconnection mode")
options.setDefault('maxEvents', 100)
options.parseArguments()

process = cms.Process("MiniAnalysisGEN")

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper
randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randSvc.populate()

# set input to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.options   = cms.untracked.PSet(
   # wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)
# Input source
#process.source = cms.Source("PoolSource",
#    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
#    fileNames = cms.untracked.vstring( 
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E0E60B0A-EFD9-E411-944B-002590494C44.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E216B380-EFD9-E411-A286-003048FEC15C.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E41672D8-EFD9-E411-96A4-02163E00F2F9.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E41DFBBB-EFD9-E411-B7E3-02163E00F319.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E4F6FEEE-EED9-E411-ABA1-02163E00F4EF.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E612966B-EFD9-E411-A4DB-02163E00EAD0.root',
#        '/store/mc/RunIIWinter15wmLHE/TT_13TeV-powheg/LHE/MCRUN2_71_V1_ext1-v1/40000/E85BB4DA-EED9-E411-8FE2-02163E00E9BC.root',
#    ),
#    inputCommands = cms.untracked.vstring('keep *', 
#        'drop LHEXMLStringProduct_*_*_*'),
#    secondaryFileNames = cms.untracked.vstring()
#)

process.source = cms.Source("EmptySource")

if options.generator == 'sherpa':
    process.generator = cms.EDFilter("SherpaGeneratorFilter",
        FetchSherpack = cms.bool(True),
        SherpaDefaultWeight = cms.double(1.0),
        SherpaParameters = cms.PSet(
            MPI_Cross_Sections = cms.vstring(' MPIs in Sherpa, Model = Amisic:', 
                ' semihard xsec = 39.8027 mb,', 
                ' non-diffractive xsec = 17.0318 mb with nd factor = 0.3142.'),
            Run = cms.vstring(' (run){', 
                ' EVENTS 1M; ERROR 0.99;', 
                ' EVENT_OUTPUT=HepMC_GenEvent[sample]', 
                ' FSF:=1.; RSF:=1.; QSF:=1.;', 
                ' SCALES METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2};', 
                ' CORE_SCALE QCD;', 
                ' METS_BBAR_MODE 5;', 
                ' NJET:=0; LJET:=2; QCUT:=20.;', 
                ' ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;', 
                ' OL_PREFIX=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/openloops/1.3.1-mlhled2', 
                ' EVENT_GENERATION_MODE Weighted;', 
                ' LOOPGEN:=OpenLoops;', 
                ' BEAM_1 2212; BEAM_ENERGY_1 = 6500.;', 
                ' BEAM_2 2212; BEAM_ENERGY_2 = 6500.;', 
                ' HARD_DECAYS On;', 
                ' STABLE[24] 0; STABLE[6] 0; WIDTH[6] 0;', 
                ' NLO_SMEAR_THRESHOLD 1;', 
                ' NLO_SMEAR_POWER 2;', 
                '}(run)', 
                ' (processes){', 
                ' Process : 93 93 ->  6 -6 93{NJET};', 
                ' Order (*,0); CKKW sqr(QCUT/E_CMS);', 
                ' NLO_QCD_Mode MC@NLO {LJET};', 
                ' ME_Generator Amegic {LJET};', 
                ' RS_ME_Generator Comix {LJET};', 
                ' Loop_Generator LOOPGEN {LJET};', 
                ' Max_N_Quarks 6 {5,6,7,8};', 
                ' Integration_Error 0.05 {5,6,7,8};', 
                ' Scales LOOSE_METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2} {5,6,7,8};', 
                ' End process', 
                '}(processes)'),
            parameterSets = cms.vstring('MPI_Cross_Sections', 
                'Run')
        ),
        SherpaPath = cms.string('./'),
        SherpaPathPiece = cms.string('./'),
        SherpaProcess = cms.string('Tops'),
        SherpaResultDir = cms.string('Result'),
        SherpackChecksum = cms.string('a43f22edc29e0bbbf1716e438e999ff0'),
        SherpackLocation = cms.string('/afs/cern.ch/work/m/mseidel/TopAnalysis/CMSSW_9_3_1/src/Sherpack/Tops/test/'),
        crossSection = cms.untracked.double(-1),
        filterEfficiency = cms.untracked.double(1.0),
        maxEventsToPrint = cms.int32(0)
    )

#pseudo-top
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
process.particleLevel.src = 'generator:unsmeared'

#tfile service
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outFilename)
                                   )


#analysis
#process.analysis = cms.EDAnalyzer("MiniAnalyzerGEN",
#                          saveTree               = cms.bool(True),
#                          savePF                 = cms.bool(True),
#                          )

process.load('TopLJets2015.TopAnalysis.miniAnalyzer_cfi')
process.analysis.saveTree = options.saveTree
process.analysis.savePF   = options.savePF
process.analysis.runOnGEN = True
process.analysis.prunedGenParticles = 'genParticles'
if not process.analysis.saveTree :
    print '\t Summary tree won\'t be saved'
if not process.analysis.savePF :
    print 'Summary PF info won\'t be saved'

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.analysis_step = cms.Path(process.particleLevel*process.analysis)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# FastGenParticleCandidateProducer broken, replace by GenParticleProducer
# (actually no idea where this is from)
process.genParticleCandidates = cms.EDProducer("GenParticleProducer",
    abortOnUnknownPDGCode = cms.untracked.bool(False),
    saveBarCodes = cms.untracked.bool(True),
    src = cms.InputTag("generatorSmeared")
)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.analysis_step,process.genfiltersummary_step,process.endjob_step)
# filter all path with the production filter sequence
for path in process.paths:
	  if path in ['lhe_step']: continue
	  getattr(process,path)._seq = process.generator * getattr(process,path)._seq 
