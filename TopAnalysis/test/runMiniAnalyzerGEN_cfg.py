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

if options.generator == 'herwig7':
    process.externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
        nEvents = cms.untracked.uint32(options.maxEvents),
        outputFile = cms.string('cmsgrid_final.lhe'),
        scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
        numberOfParameters = cms.uint32(1),
        args = cms.vstring('/afs/cern.ch/work/m/mseidel/generator/CMSSW_7_1_20/src/TT_hdamp_TuneT4_noweights_NNPDF30_13TeV_powheg_hvq.tgz')
    )
    
    process.generator = cms.EDFilter("Herwig7GeneratorFilter",
        run = cms.string('InterfaceTest'),
       # dumpConfig = cms.untracked.string('HerwigConfig.in'),
        repository = cms.string('${HERWIGPATH}/HerwigDefaults.rpo'),
        dataLocation = cms.string('${HERWIGPATH}'),
        generatorModule = cms.string('/Herwig/Generators/EventGenerator'),
        eventHandlers = cms.string('/Herwig/EventHandlers'),
        hwpp_cmsDefaults = cms.vstring('+hwpp_basicSetup', 
            '+hwpp_setParticlesStableForDetector'),
        hwpp_basicSetup = cms.vstring('#create ThePEG::RandomEngineGlue /Herwig/RandomGlue', 
            '#set /Herwig/Generators/LHCGenerator:RandomNumberGenerator /Herwig/RandomGlue', 
            'set /Herwig/Generators/EventGenerator:NumberOfEvents -1', 
            'set /Herwig/Generators/EventGenerator:DebugLevel 2', 
            'set /Herwig/Generators/EventGenerator:PrintEvent 1', 
            'set /Herwig/Generators/EventGenerator:MaxErrors 10000'),
        hwpp_setParticlesStableForDetector = cms.vstring('set /Herwig/Particles/mu-:Stable Stable', 
            'set /Herwig/Particles/mu+:Stable Stable', 
            'set /Herwig/Particles/Sigma-:Stable Stable', 
            'set /Herwig/Particles/Sigmabar+:Stable Stable', 
            'set /Herwig/Particles/Lambda0:Stable Stable', 
            'set /Herwig/Particles/Lambdabar0:Stable Stable', 
            'set /Herwig/Particles/Sigma+:Stable Stable', 
            'set /Herwig/Particles/Sigmabar-:Stable Stable', 
            'set /Herwig/Particles/Xi-:Stable Stable', 
            'set /Herwig/Particles/Xibar+:Stable Stable', 
            'set /Herwig/Particles/Xi0:Stable Stable', 
            'set /Herwig/Particles/Xibar0:Stable Stable', 
            'set /Herwig/Particles/Omega-:Stable Stable', 
            'set /Herwig/Particles/Omegabar+:Stable Stable', 
            'set /Herwig/Particles/pi+:Stable Stable', 
            'set /Herwig/Particles/pi-:Stable Stable', 
            'set /Herwig/Particles/K+:Stable Stable', 
            'set /Herwig/Particles/K-:Stable Stable', 
            'set /Herwig/Particles/K_S0:Stable Stable', 
            'set /Herwig/Particles/K_L0:Stable Stable'),
        hwpp_lhe = cms.vstring('library LesHouches.so', 
            'cd /Herwig/EventHandlers',
            'create ThePEG::LesHouchesFileReader myReader',
            'set myReader:FileName ttbar.lhe',
            'create ThePEG::Cuts NoCuts',
            'set myReader:Cuts NoCuts',
            'create ThePEG::LesHouchesEventHandler myLesHouchesHandler',
            'set myLesHouchesHandler:CascadeHandler /Herwig/Shower/ShowerHandler',
            'set myLesHouchesHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
            'set myLesHouchesHandler:DecayHandler /Herwig/Decays/DecayHandler',
            'set myLesHouchesHandler:PartonExtractor /Herwig/MatrixElements/SubProcess:PartonExtractor',
            'insert myLesHouchesHandler:PostSubProcessHandlers 0 /Herwig/QEDRadiation/QEDRadiationHandler',
            'insert myLesHouchesHandler:LesHouchesReaders 0 myReader',
            'cd /Herwig/Generators',
            'cp EventGenerator myLesHouchesGenerator',
            'set myLesHouchesGenerator:EventHandler /Herwig/EventHandlers/myLesHouchesHandler',
            'saverun myLesHouches myLesHouchesGenerator',
        ),
        h7_PPCollider = cms.vstring(
            '# -*- ThePEG-repository -*-',
            'cd /Herwig/EventHandlers',
            'create ThePEG::FixedCMSLuminosity Luminosity FixedCMSLuminosity.so',
            'set Luminosity:Energy 13000.0',
            'set EventHandler:LuminosityFunction Luminosity',
            'set EventHandler:BeamA /Herwig/Particles/p+',
            'set EventHandler:BeamB /Herwig/Particles/p+',
            'cd /Herwig/Cuts',
            '# create the cuts object for hadron collisions',
            'set Cuts:ScaleMin 2.0*GeV2',
            'set Cuts:X1Min 1.0e-5',
            'set Cuts:X2Min 1.0e-5',
            '# Matchbox settings',
            'cd /Herwig/MatrixElements/Matchbox',
            'set Factory:FirstPerturbativePDF Yes',
            'set Factory:SecondPerturbativePDF Yes',
            'set Factory:PartonExtractor /Herwig/Partons/PPExtractor',
            'cd /Herwig/Merging',
            'set MergingFactory:FirstPerturbativePDF Yes',
            'set MergingFactory:SecondPerturbativePDF Yes',
            'set MergingFactory:PartonExtractor /Herwig/Partons/PPExtractor',
            'cd /Herwig/Generators/',
            'set /Herwig/MatrixElements/SubProcess:PartonExtractor /Herwig/Partons/PPExtractor',
            '# Read in parameters to use the soft model',
            '#read snippets/SoftModel.in',
            '+h7_SoftModel',
        ),
        h7_SoftModel = cms.vstring(
            '# Parameters for soft interactions',
            'set /Herwig/Partons/RemnantDecayer:colourDisrupt 0.0',
            'set /Herwig/Hadronization/ColourReconnector:ColourReconnection Yes',
            '# Use multiperipheral kinematics',
            'set /Herwig/Partons/RemnantDecayer:MultiPeriph Yes',
            '# Set gaussian width of longitudinal momentum fraction fluctuation',
            'set /Herwig/Partons/RemnantDecayer:gaussWidth 0.03',
            '# Tuned to min-bias data (5B tune)',
            'set /Herwig/Hadronization/ColourReconnector:ReconnectionProbability 0.652710',
            'set /Herwig/UnderlyingEvent/MPIHandler:pTmin0 3.568157',
            'set /Herwig/UnderlyingEvent/MPIHandler:InvRadius 1.489997',
            'set /Herwig/UnderlyingEvent/MPIHandler:Power 0.420445',
            'set /Herwig/Partons/RemnantDecayer:ladderPower -0.088983',
            'set /Herwig/Partons/RemnantDecayer:ladderNorm 1.086029',
        ),
        h7_LHE = cms.vstring(
            '##################################################',
            '# Example generator based on LHC parameters',
            '# usage: Herwig read LHE.in',
            '##################################################',
            '#read snippets/PPCollider.in',
            '+h7_PPCollider',
            '##################################################',
            '# Technical parameters for this run',
            '##################################################',
            'cd /Herwig/Generators',
            '#set EventGenerator:NumberOfEvents 10000000',
            '#set EventGenerator:RandomNumberGenerator:Seed 31122001',
            '#set EventGenerator:DebugLevel 0',
            '#set EventGenerator:PrintEvent 10',
            '#set EventGenerator:MaxErrors 10000',
            '##################################################',
            '#   Create the Les Houches file handler and reader',
            '##################################################',
            'cd /Herwig/EventHandlers',
            'library LesHouches.so',
            '# create the event handler',
            'create ThePEG::LesHouchesEventHandler LesHouchesHandler',
            '# set the various step handlers',
            'set LesHouchesHandler:PartonExtractor /Herwig/Partons/PPExtractor',
            'set LesHouchesHandler:CascadeHandler /Herwig/Shower/ShowerHandler',
            'set LesHouchesHandler:DecayHandler /Herwig/Decays/DecayHandler',
            'set LesHouchesHandler:HadronizationHandler /Herwig/Hadronization/ClusterHadHandler',
            '# set the weight option (e.g. for MC@NLO)',
            'set LesHouchesHandler:WeightOption VarNegWeight',
            '# set event hander as one to be used',
            'set /Herwig/Generators/EventGenerator:EventHandler /Herwig/EventHandlers/LesHouchesHandler',
            '# Set up an EMPTY CUTS object',
            '# Normally you will have imposed any cuts you want',
            '# when generating the event file and do not want any more',
            '# in particular for POWHEG and MC@NLO you must not apply cuts on the',
            '# the extra jet',
            'create ThePEG::Cuts /Herwig/Cuts/NoCuts',
            '# You may wish to use the same PDF as the events were generated with',
            'create ThePEG::LHAPDF /Herwig/Partons/LHAPDF ThePEGLHAPDF.so',
            'set /Herwig/Partons/LHAPDF:PDFName NNPDF30_nlo_as_0118',
            'set /Herwig/Partons/LHAPDF:RemnantHandler /Herwig/Partons/HadronRemnants',
            'set /Herwig/Particles/p+:PDF /Herwig/Partons/LHAPDF',
            'set /Herwig/Particles/pbar-:PDF /Herwig/Partons/LHAPDF',
            'set /Herwig/Partons/PPExtractor:FirstPDF  /Herwig/Partons/LHAPDF',
            'set /Herwig/Partons/PPExtractor:SecondPDF /Herwig/Partons/LHAPDF',
            '# We would recommend the shower uses the default PDFs with which it was tuned.',
            '# However it can be argued that the same set as for the sample should be used for',
            '# matched samples, i.e. MC@NLO (and less so POWHEG)',
            '#set /Herwig/Shower/ShowerHandler:PDFA /Herwig/Partons/LHAPDF',
            '#set /Herwig/Shower/ShowerHandler:PDFB /Herwig/Partons/LHAPDF',
            '# You can in principle also change the PDFs for the remnant extraction and',
            '# multiple scattering. As the generator was tuned with the default values',
            '# this is STRONGLY DISCOURAGED without retuning the MPI parameters',
            '# create the reader and set cuts',
            'create ThePEG::LesHouchesFileReader LesHouchesReader',
            'set LesHouchesReader:FileName cmsgrid_final.lhe',
            'set LesHouchesReader:AllowedToReOpen No',
            'set LesHouchesReader:InitPDFs 0',
            'set LesHouchesReader:Cuts /Herwig/Cuts/NoCuts',
            '# option to ensure momentum conservation is O.K. due rounding errors (recommended)',
            'set LesHouchesReader:MomentumTreatment RescaleEnergy',
            '# set the pdfs',
            'set LesHouchesReader:PDFA /Herwig/Partons/LHAPDF',
            'set LesHouchesReader:PDFB /Herwig/Partons/LHAPDF',
            '# if using BSM models with QNUMBER info',
            '#set LesHouchesReader:QNumbers Yes',
            '#set LesHouchesReader:Decayer /Herwig/Decays/Mambo',
            '# and add to handler',
            'insert LesHouchesHandler:LesHouchesReaders 0 LesHouchesReader',
            '##################################################',
            '#  Shower parameters',
            '##################################################',
            '# normally, especially for POWHEG, you want',
            '# the scale supplied in the event files (SCALUP)',
            '# to be used as a pT veto scale in the parton shower',
            'set /Herwig/Shower/ShowerHandler:MaxPtIsMuF Yes',
            'set /Herwig/Shower/ShowerHandler:RestrictPhasespace Yes',
            '# Shower parameters',
            '# treatment of wide angle radiation',
            'set /Herwig/Shower/PartnerFinder:PartnerMethod Random',
            'set /Herwig/Shower/PartnerFinder:ScaleChoice Partner',
            '# fix issue before 7.0.5 (not needed after this)',
            'set /Herwig/Shower/GtoQQbarSplitFn:AngularOrdered Yes',
            'set /Herwig/Shower/GammatoQQbarSplitFn:AngularOrdered Yes',
            '# with MC@NLO these parameters are required for consistency of the subtraction terms',
            '# suggested parameters (give worse physics results with POWHEG)',
            '#set /Herwig/Shower/KinematicsReconstructor:InitialInitialBoostOption LongTransBoost',
            '#set /Herwig/Shower/KinematicsReconstructor:ReconstructionOption General',
            '#set /Herwig/Shower/KinematicsReconstructor:FinalStateReconOption Default',
            '#set /Herwig/Shower/KinematicsReconstructor:InitialStateReconOption Rapidity',
            '#set /Herwig/Shower/ShowerHandler:SpinCorrelations No',
            '##################################################',
            '# LHC physics parameters (override defaults here) ',
            '##################################################',
            '# e.g if different top mass used',
            'set /Herwig/Particles/t:NominalMass 172.5',
            '##################################################',
            '# Save run for later usage with `Herwig run`',
            '##################################################',
            'cd /Herwig/Generators',
            'saverun LHE EventGenerator',
        ),
        configFiles = cms.vstring(),
        crossSection = cms.untracked.double(-1),
        parameterSets = cms.vstring('hwpp_cmsDefaults', 'h7_LHE'),
        filterEfficiency = cms.untracked.double(1.0),
        dummyprocess = cms.vstring('insert /Herwig/MatrixElements/SimpleQCD:MatrixElements[0] /Herwig/MatrixElements/MEPP2ttbarH')
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

process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")

process.rivetAnalyzer.AnalysisNames = cms.vstring(
    #'CMS_2016_I1434354', # diff xs lepton+jets
    #'MC_TTBAR', # MC analysis for lepton+jets
    'MC_TOPMASS_LJETS', # MC analysis for lepton+jets top mass
    #'CMS_LesHouches2015', # MC analysis for dilepton
    #'MC_GENERIC', # MC generic analysis
    #'MC_XS', # MC xs analysis
)
process.rivetAnalyzer.OutputFile      = 'run.yoda'
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator:unsmeared")
process.rivetAnalyzer.CrossSection    = 831.76 # NNLO (arXiv:1303.6254)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.analysis_step = cms.Path(process.particleLevel*process.analysis*process.rivetAnalyzer)
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
