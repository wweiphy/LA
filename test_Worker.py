# cmsRun test_Worker.py maxEvents=5000 notInPCL=False PCLoutName=PCL_worker.root inputFiles=/store/data/Run2018D/SingleMuon/ALCARECO/SiPixelCalSingleMuon-ForPixelALCARECO_UL2018-v1/230000/F33B7CA6-256A-B34F-9536-7594FBC6F75B.root

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.register( "outName", "Tree.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "name and path of the Tree output files (without extension)" )
options.register( "PCLoutName", "DQM_Worker.root", VarParsing.multiplicity.singleton, VarParsing.varType.string, "name and path of the PCL Worker output files (without extension)" )

options.register( "notInPCL", False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, "Is it in PCL or not ?" )
options.parseArguments()

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("LA", Run2_2018)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = "120X_dataRun2_v2" # for CMSSW 12_0_0_pre

process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFindingBeamSpot
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


process.load("CalibTracker.SiPixelLorentzAngle.SiPixelLorentzAnglePCLWorker_cfi")
process.SiPixelLorentzAnglePCLWorker.folder = cms.string('AlCaReco/SiPixelLorentzAngleHarvesting/')
process.SiPixelLorentzAnglePCLWorker.fileName = cms.string(options.outName)
process.SiPixelLorentzAnglePCLWorker.notInPCL = cms.bool(options.notInPCL)


process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
                                     fileName = cms.untracked.string(options.PCLoutName))
                                     
process.p = cms.Path(process.offlineBeamSpot*
                     process.MeasurementTrackerEvent*
                     process.TrackRefitter*
                     process.SiPixelLorentzAnglePCLWorker)

process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(
    process.p,
    process.DQMoutput_step
    )
                                     
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
    
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

#process.options.numberOfThreads=cms.untracked.uint32(4)
