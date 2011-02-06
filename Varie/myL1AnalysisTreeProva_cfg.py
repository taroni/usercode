
#######################################################
# CONFIGURATION FILE TO ANALYSE OUTPUT OF SIMULATIONS #
#######################################################

#####################
# preliminary setup #
#####################
#
# definition of CMS interface to framework
import FWCore.ParameterSet.Config as cms               ### REMEMBER THIS DEFINITION

######################################################
######################################################
###                                                ###
###                                                ###
###    DEFINITION OF ANALYSIS PLUGIN PARAMETERS    ###
###                                                ###
###                                                ###
######################################################
######################################################

########################################
# which data are going to be analysed? #
########################################
#
# definition of process
# in this case, Physics Analysis Toolkit template
process = cms.Process("PATtemplate")
#
# source files to feed the analyser
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_0.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_1.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_2.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_3.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_4.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_5.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_6.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_7.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_8.root",
"rfio:/castor/cern.ch/cms/store/caf/user/taroni/testVari/myHWW2mu2nu_9.root"

)
)
#
###################################
##################
# geometry setup #
##################
#
# magnetic field
process.load("Configuration.StandardSequences.MagneticField_40T_cff")
process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = True
process.TrackerDigiGeometryESModule = cms.ESProducer("TrackerDigiGeometryESModule",
    fromDDD = cms.bool(True),
    applyAlignment = cms.bool(False),
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string('')
)
#
# longbarrel geometry
from Geometry.CMSCommonData.cmsIdealGeometryXML_cfi import *
process.load("SLHCUpgradeSimulations.Geometry.longbarrel_cmsIdealGeometryXML_cff")
#
process.TrackerGeometricDetESModule = cms.ESProducer("TrackerGeometricDetESModule",
    fromDDD = cms.bool(True)
)


process.load("SLHCUpgradeSimulations.Utilities.StackedTrackerGeometry_cfi")

###################################
#
# define number of events to be processed (-1 means all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#
# define some input tags for the analyzer used in this module
process.PATtemplate = cms.EDAnalyzer("myL1TrackAnalyzerTree",
###process.PATtemplate = cms.EDAnalyzer("TestL1Track",
  trackletVTX = cms.string("center"),
##  trackletVTX = cms.string("offcenter"),
  seedSuperLayer = cms.int32(0),
  numberStubs = cms.int32(6),
  windowSize = cms.double(99.0),
  tightStubsL1Trk = cms.int32(0),
#
   src = cms.InputTag("siPixelRecHits"),
#
   stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
   rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
   matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
   simtrackHits = cms.InputTag("famosSimHits"),
#
   locStubsPixelDigis = cms.InputTag("LocalStubsFromPixelDigis"),
   locStubsSimHits = cms.InputTag("LocalStubsFromSimHits"),
   globStubsPixelDigis = cms.InputTag("GlobalStubsFromPixelDigis"),
   globStubsSimHits = cms.InputTag("GlobalStubsFromSimHits"),
#
   trackletsPixelDigis = cms.InputTag("TrackletsFromPixelDigis","ShortTracklets"),
   trackletsSimHits = cms.InputTag("TrackletsFromSimHits","ShortTracklets"),
#
   l1tracksPixelDigis = cms.InputTag("L1TracksFromPixelDigis","L1Tracks"),
   l1tracksSimHits = cms.InputTag("L1TracksFromSimHits","L1Tracks"),
#
   trackProducer = cms.InputTag("ctfWithMaterialTracks"),
#
    L1EGamma       = cms.InputTag("L1ExtraMaker","EGamma"),
    L1IsoEGamma    = cms.InputTag("L1ExtraMaker","IsoEGamma"),
    L1Tau          = cms.InputTag("L1ExtraMaker","Taus"),
    L1IsoTau       = cms.InputTag("L1ExtraMaker","IsoTaus"),
    Jets           = cms.InputTag("L1ExtraMaker","Jets"),
#
    genParticles   = cms.InputTag("genParticles"),


   # for using track hit association
   associatePixel = cms.bool(True),
   associateStrip = cms.bool(False),
   associateRecoTracks = cms.bool(False),
   ROUList = cms.vstring('famosSimHitsTrackerHits')
)

################
# define output #
################
# 
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('outfileL1.root' )
                                   )
#
# define the global process for cmsRun
process.p = cms.Path( process.PATtemplate )

#############################
# END OF CONFIGURATION FILE #
#############################

