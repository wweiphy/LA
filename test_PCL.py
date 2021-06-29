# cmsRun test_PCL.py maxEvents=1000 treeout=Tree.root inputFiles=file:/uscms/home/wwei/nobackup/LA/CMSSW_11_3_0/src/SiPixelCalSingleMuon.root


import FWCore.ParameterSet.Config as cms
import sys
import os

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register( "treeout", "Tree.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "name and path of the output tree files (without extension)" )
options.parseArguments()

#if options.maxEvents is -1: # maxEvents is set in VarParsing class by default to -1
#    options.maxEvents = 100

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("LA", Run2_2018)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# check for the correct tag on https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions

process.load("DQMServices.Core.DQM_cfg")
process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("CondCore.CondDB.CondDB_cfi")

process.GlobalTag.globaltag = "113X_dataRun2_v6"
#process.GlobalTag.globaltag = "112X_dataRun2_v7" CMSSW 11_2_0_pre10
# https://github.com/cms-sw/cmssw/blob/CMSSW_11_3_X/Configuration/AlCa/python/autoCond.py

process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *
process.offlineBeamSpot = offlineBeamSpot

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

process.MeasurementTrackerEvent.pixelClusterProducer = 'ALCARECOSiPixelCalSingleMuon'
process.MeasurementTrackerEvent.stripClusterProducer = 'ALCARECOSiPixelCalSingleMuon'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag()
process.MeasurementTrackerEvent.inactiveStripDetectorLabels = cms.VInputTag()


process.load("RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilderWithoutRefit_cfi")
process.TrackRefitter.src = 'ALCARECOSiPixelCalSingleMuon'
process.TrackRefitter.TrajectoryInEvent = True

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('simul', 
        'cout'),
    simul = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    ),
)

process.SiPixelLorentzAnglePCLWorker = cms.EDProducer(
    "SiPixelLorentzAnglePCLWorker",
    folder = cms.string('CalibTracker/SiPixelLorentzAngle'),
    fileName = cms.string(options.treeout),
    src = cms.InputTag("TrackRefitter"),
    binsDepth    = cms.int32(50),
    binsDrift =    cms.int32(200),
    ptMin = cms.double(3),
    normChi2Max = cms.double(2),
    clustSizeYMin = cms.int32(4),
    residualMax = cms.double(0.005),
    clustChargeMax = cms.double(120000)
)

process.CondDB.connect = 'sqlite_file:Harvester.db'

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDB,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string("SiPixelLorentzAngleRcd"),
        tag = cms.string('Harvester')
    ))
)

process.PixelLorentzAnglePCLHarvester = cms.EDProducer(
    "SiPixelLorentzAnglePCLHarvester",
    dqmDir = cms.string('CalibTracker/SiPixelLorentzAngle'),
    record = cms.string("SiPixelLorentzAngleRcd"),
    fitProbCut = cms.double(0.5)
)
                                     
process.p = cms.Path(process.offlineBeamSpot*
                     process.MeasurementTrackerEvent*
                     process.TrackRefitter*
                     process.SiPixelLorentzAnglePCLWorker*
                     process.PixelLorentzAnglePCLHarvester)

process.dqmsave_step = cms.Path(process.dqmSaver)

process.schedule = cms.Schedule(
    process.p,
    process.dqmsave_step
    )
    
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(int(options.maxEvents))
)

#process.options.numberOfThreads=cms.untracked.uint32(8)


process.source = cms.Source("PoolSource",
	#put here the sample you want to use
#    fileNames = cms.untracked.vstring('file:/uscms/home/wwei/nobackup/LA/CMSSW_11_2_0_pre10/src/trial2/SiPixelCalSingleMuon_1.root'),
    fileNames = cms.untracked.vstring(options.inputFiles),
#   skipEvents = cms.untracked.uint32(100)
#
)
process.dqmSaver.workflow = '/CalibTracker/SiPixelLorentzAngle/test'
