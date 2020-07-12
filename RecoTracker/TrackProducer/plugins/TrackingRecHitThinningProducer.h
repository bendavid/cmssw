#ifndef TrackProducer_TrackingRecHitThinningProducer_h
#define TrackProducer_TrackingRecHitThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "DataFormats/Common/interface/ThinnedRefSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"

namespace edm {
  class ParameterSetDescription;
}

class TrackingRecHitSelector {
public:
  TrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  void preChoose(edm::Handle<TrackingRecHitCollection> hits, edm::Event const& event, edm::EventSetup const& es);

  bool choose(unsigned int iIndex, TrackingRecHit const& iItem) const;

  void reset() { keysToSave_.clear(); }

private:
  edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;
  edm::ThinnedRefSet<TrackingRecHitCollection> keysToSave_;
};

using TrackingRecHitThinningProducer = edm::ThinningProducer<TrackingRecHitCollection, TrackingRecHitSelector>;

#endif
