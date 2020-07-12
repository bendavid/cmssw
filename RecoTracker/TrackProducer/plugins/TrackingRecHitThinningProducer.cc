#include "RecoTracker/TrackProducer/plugins/TrackingRecHitThinningProducer.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/Common/interface/RefItemGet.h"

TrackingRecHitSelector::TrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc)
     : trackExtraToken_(cc.consumes<reco::TrackExtraCollection>(pset.getParameter<edm::InputTag>("trackExtraTag"))) {}

void TrackingRecHitSelector::fillDescription(edm::ParameterSetDescription & desc) {
  desc.add<edm::InputTag>("trackExtraTag");  
}

void TrackingRecHitSelector::preChoose(edm::Handle<TrackingRecHitCollection> hits, edm::Event const& event, edm::EventSetup const& es) {
  auto trackExtras = event.getHandle(trackExtraToken_);
  
  keepMask_.clear();
  keepMask_.resize(hits->size(), false);
  
  for (const auto &trackExtra : *trackExtras) {
    for (unsigned int i=0; i<trackExtra.recHitsSize(); ++i) {
      //make sure reference points to desired product id in case of (multiple) thinned collections
      //TODO this could be more efficient with a RefVector version of thinnedRefFrom
      const TrackingRecHitRef& hitRef = edm::thinnedRefFrom(trackExtra.recHit(i), hits.id(), event.productGetter());
      if (hitRef.isNonnull()) {
        keepMask_[hitRef.key()] = true;
      }
    }
  }
}

bool TrackingRecHitSelector::choose(unsigned int iIndex, TrackingRecHit const& iItem) {
  return keepMask_[iIndex];
}
