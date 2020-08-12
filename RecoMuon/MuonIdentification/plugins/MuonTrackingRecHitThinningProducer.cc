#include "RecoMuon/MuonIdentification/plugins/MuonTrackingRecHitThinningProducer.h"
#include "DataFormats/MuonReco/interface/Muon.h"

MuonTrackingRecHitSelector::MuonTrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc)
    : cut_(pset.getParameter<std::string>("cut")),
      muonToken_(cc.consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"))),
      selector_(cut_) {}

void MuonTrackingRecHitSelector::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<std::string>("cut", "");
  desc.add<edm::InputTag>("muonTag", edm::InputTag("muons"));
}

void MuonTrackingRecHitSelector::preChoose(edm::Handle<TrackingRecHitCollection> hits,
                                           edm::Event const& event,
                                           edm::EventSetup const& es) {
  auto muons = event.getHandle(muonToken_);

  auto filler = keysToSave_.fill(edm::RefProd(hits), event.productGetter());
  for (const auto& muon : *muons) {
    if (!selector_(muon)) {
      continue;
    }
    reco::Track const& track = *muon.bestTrack();
    for (unsigned int i = 0; i < track.recHitsSize(); ++i) {
      filler.insert(track.recHit(i));
    }
  }
}

bool MuonTrackingRecHitSelector::choose(unsigned int iIndex, TrackingRecHit const& iItem) const {
  return keysToSave_.contains(iIndex);
}
