#ifndef TrackerDetGunProducer_H
#define TrackerDetGunProducer_H


#include "IOMC/ParticleGuns/interface/BaseFlatGunProducer.h"
#include "DataFormats/DetId/interface/DetId.h"


namespace edm
{
  
  class TrackerDetGunProducer : public BaseFlatGunProducer
  {
  
  public:
    TrackerDetGunProducer(const ParameterSet & pset);
    ~TrackerDetGunProducer() override;
   
    void produce(Event & e, const EventSetup& es) override;

  private:
    
    // data members
    std::vector<double> fProbParticle_;
    double            fMinPt   ;
    double            fMaxPt   ;
    double            fMinX;
    double            fMaxX;
    double            fMinY;
    double            fMaxY;
    DetId detId_;

  };
} 

#endif
