import FWCore.ParameterSet.Config as cms
import sys
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()

process = cms.Process('ALCAHARVEST')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.AlCaHarvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

                     
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('ALCAPROMPT_files.txt')
readFiles = cms.untracked.vstring( *mylist)

process.source = cms.Source("PoolSource",
    fileNames = readFiles,
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


process.PoolDBOutputService.toPut.append(process.ALCAHARVESTSiPixelLA_dbOutput)

process.pclMetadataWriter.recordsToMap.append(process.ALCAHARVESTSiPixelLA_metadata)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "121X_dataRun2_v10"


process.SiPixelLA      = cms.Path(process.ALCAHARVESTSiPixelLorentzAngle)


process.ALCAHARVESTDQMSaveAndMetadataWriter = cms.Path(process.dqmSaver+process.pclMetadataWriter)


process.schedule = cms.Schedule(process.SiPixelLA,
                                process.ALCAHARVESTDQMSaveAndMetadataWriter)

