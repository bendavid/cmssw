#include "RecoMuon/MuonIdentification/plugins/MuonTrackExtraThinningProducer.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/RefItemGet.h"

MuonTrackExtraSelector::MuonTrackExtraSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc)
    : cut_(pset.getParameter<std::string>("cut")),
      muonToken_(cc.consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"))),
      selector_(cut_) {}

void MuonTrackExtraSelector::fillDescription(edm::ParameterSetDescription & desc) {
  desc.add<std::string>("cut");
  desc.add<edm::InputTag>("muonTag");
}

void MuonTrackExtraSelector::preChoose(edm::Handle<reco::TrackExtraCollection> trackExtras, edm::Event const& event, edm::EventSetup const& es) {
  auto muons = event.getHandle(muonToken_);
  
  keepMask_.clear();
  keepMask_.resize(trackExtras->size(), false);
  
  for (const auto &muon : *muons) {
    //make sure reference points to desired product id in case of (multiple) thinned collections
    printf("getting thinnedref\n");
    const edm::Ref<reco::TrackExtraCollection>& trackExtraRef = edm::thinnedRefFrom(muon.bestTrack()->extra(), trackExtras.id(), event.productGetter());
    printf("got thinnedref\n");
    if (trackExtraRef.isNonnull()) {
      keepMask_[trackExtraRef.key()] = true;
    }
  }
}

bool MuonTrackExtraSelector::choose(unsigned int iIndex, reco::TrackExtra const& iItem) {
  return keepMask_[iIndex];
}
