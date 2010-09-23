import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
import FWCore.Framework.test.cmsExceptionsFatal_cff

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageService.test.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

from GluonFusion.GenQuantities.datafiles2_cfi import *

process.source = cms.Source("PoolSource",
                            fileNames           = readFiles, 
                             secondaryFileNames = secFiles
##     fileNames = cms.untracked.vstring(
##     'file:/data06/users/taroni/trigger/HtoWW/H140GGfusiontoWW_cfi_py_GEN-1000.root'
##    ),
)
process.MessageLogger = cms.Service("MessageLogger",
         destinations = cms.untracked.vstring(
                        "cout",
                        #"detailedINFO100000"
                        ),
         cout = cms.untracked.PSet(
         #threshold = cms.untracked.string('DEBUG')
         #),
         #detailedINFO100000 = cms.untracked.PSet(
         threshold = cms.untracked.string('INFO')
         ),
                                    
)

process.demo = cms.EDProducer('GenQuantities',
##      filename = cms.untracked.string('H140toWWto2Mu2NuHistos-100000.root'),
##      outfile  = cms.untracked.string('MuonsTree-100000.root')
## ##     decayChain = cms.string('HtoWWto2Mu2Nu')
     filename = cms.untracked.string('HtoWWGenHisto.root'),
     outfile  = cms.untracked.string('HtoWWGenTree.root')
)

## process.output = cms.OutputModule("PoolOutputModule",
##            fileName = cms.untracked.string('prova.root'),
##             outputCommands = cms.untracked.vstring(
##             'drop *',
##             'keep *_demo_*_*')
## ##             SelectEvents = cms.untracked.PSet(
## ##                  SelectEvents = cms.vstring('Demo')
## ##              ) 
##          )
 
## process.p = cms.Path(process.demo+process.output)
process.p = cms.Path(process.demo)
