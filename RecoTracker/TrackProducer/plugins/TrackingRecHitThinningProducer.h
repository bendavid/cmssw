#ifndef MuonIdentification_TrackingRecHitThinningProducer_h
#define MuonIdentification_TrackingRecHitThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

namespace edm {
  class ParameterSetDescription;
}

class TrackingRecHitSelector {

public:
  TrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);
  static void fillDescription(edm::ParameterSetDescription & desc);
  void preChoose(edm::Handle<TrackingRecHitCollection> hits, edm::Event const& event, edm::EventSetup const& es);
  bool choose(unsigned int iIndex, TrackingRecHit const& iItem);
  
private:
  edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;
  std::vector<bool> keepMask_;
};

typedef edm::ThinningProducer<TrackingRecHitCollection, TrackingRecHitSelector> TrackingRecHitThinningProducer;

#endif
