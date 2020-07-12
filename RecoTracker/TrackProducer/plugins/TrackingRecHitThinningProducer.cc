#include "RecoTracker/TrackProducer/plugins/TrackingRecHitThinningProducer.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

TrackingRecHitSelector::TrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc)
    : trackExtraToken_(cc.consumes<reco::TrackExtraCollection>(pset.getParameter<edm::InputTag>("trackExtraTag"))) {}

void TrackingRecHitSelector::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<edm::InputTag>("trackExtraTag");
}

void TrackingRecHitSelector::preChoose(edm::Handle<TrackingRecHitCollection> hits,
                                       edm::Event const& event,
                                       edm::EventSetup const& es) {
  auto trackExtras = event.getHandle(trackExtraToken_);

  auto filler = keysToSave_.fill(edm::RefProd(hits), event.productGetter());
  for (const auto& trackExtra : *trackExtras) {
    for (unsigned int i = 0; i < trackExtra.recHitsSize(); ++i) {
      filler.insert(trackExtra.recHit(i));
    }
  }
}

bool TrackingRecHitSelector::choose(unsigned int iIndex, TrackingRecHit const& iItem) const {
  return keysToSave_.contains(iIndex);
}
