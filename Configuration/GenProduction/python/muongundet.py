import FWCore.ParameterSet.Config as cms


generator = cms.EDProducer("TrackerDetGunProducer",
    PGunParameters = cms.PSet(
        #MaxPt = cms.double(150.01),
        #MinPt = cms.double(2.99),
        MaxPt = cms.double(150.00),
        MinPt = cms.double(150.00),
        PartID = cms.vint32([-13,13]),
        ProbParts = cms.vdouble([0.5,0.5]),
        MinEta = cms.double(-2.4),
        MaxEta = cms.double(-2.3),
        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),
        detId = cms.uint32(470098246),
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single mu plus and minus pt 3 to 150'),
    AddAntiParticle = cms.bool(False),
    firstRun = cms.untracked.uint32(1)
)
    
