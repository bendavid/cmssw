#ifndef IOMC_ParticleGuns_FlatRandomMultiParticlePtGunProducer_H
#define IOMC_ParticleGuns_FlatRandomMultiParticlePtGunProducer_H

#include <vector>
#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"

namespace edm {

  class FlatRandomMultiParticlePtGunProducer : public BaseFlatGunProducer {
  
  public:
    FlatRandomMultiParticlePtGunProducer(const ParameterSet & pset);
    ~FlatRandomMultiParticlePtGunProducer() override;
    
    void produce(Event &e, const EventSetup& es) override;

  private:
    
    // data members
    std::vector<double> fProbParticle_;
    double              fMinPt;
    double              fMaxPt;
    
  };
} 
#endif
