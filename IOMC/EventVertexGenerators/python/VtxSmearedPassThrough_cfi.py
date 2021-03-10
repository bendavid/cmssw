import FWCore.ParameterSet.Config as cms

from IOMC.EventVertexGenerators.VtxSmearedParameters_cfi import VtxSmearedCommon
VtxSmeared = cms.EDProducer("PassThroughEvtVtxGenerator",
    VtxSmearedCommon
)



