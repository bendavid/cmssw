#ifndef RecoLuminosity_LumiProducer_LumiInfoRunHeaderProducer_h
#define RecoLuminosity_LumiProducer_LumiInfoRunHeaderProducer_h


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <vector>


class LumiInfoRunHeaderProducer : public edm::one::EDProducer<> {
  public:
  
  LumiInfoRunHeaderProducer(const edm::ParameterSet&);  
  ~LumiInfoRunHeaderProducer() {}
    
  //virtual void produce(edm::Event &, edm::EventSetup const&) override;
  virtual void endRunProduce(edm::Run &, edm::EventSetup const&) override;
    
  private:
    bool fillSchemeFromConfig_;
    int bunchSpacingFromConfig_;
    
};
#endif


