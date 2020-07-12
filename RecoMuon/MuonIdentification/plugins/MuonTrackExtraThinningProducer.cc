#include "RecoMuon/MuonIdentification/plugins/MuonTrackExtraThinningProducer.h"
#include "DataFormats/MuonReco/interface/Muon.h"

MuonTrackExtraSelector::MuonTrackExtraSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc)
    : cut_(pset.getParameter<std::string>("cut")),
      slimTrajParams_(pset.getParameter<bool>("slimTrajParams")),
      slimResiduals_(pset.getParameter<bool>("slimResiduals")),
      muonToken_(cc.consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"))),
      selector_(cut_) {}

void MuonTrackExtraSelector::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<std::string>("cut", "");
  desc.add<bool>("slimTrajParams", false);
  desc.add<bool>("slimResiduals", false);
  desc.add<edm::InputTag>("muonTag", edm::InputTag("muons"));
}

void MuonTrackExtraSelector::preChoose(edm::Handle<reco::TrackExtraCollection> trackExtras,
                                       edm::Event const& event,
                                       edm::EventSetup const& es) {
  auto muons = event.getHandle(muonToken_);

  auto filler = keysToSave_.fill(edm::RefProd(trackExtras), event.productGetter());
  for (const auto& muon : *muons) {
    if (!selector_(muon)) {
      continue;
    }
    filler.insert(muon.bestTrack()->extra());
  }
}

std::optional<reco::TrackExtra> MuonTrackExtraSelector::choose(unsigned int iIndex,
                                                               reco::TrackExtra const& iItem) const {
  if (!keysToSave_.contains(iIndex)) {
    return std::nullopt;
  }

  reco::TrackExtra trackExtra = iItem;
  if (slimTrajParams_) {
    trackExtra.setTrajParams(reco::TrackExtraBase::TrajParams(), reco::TrackExtraBase::Chi2sFive());
  }
  if (slimResiduals_) {
    trackExtra.setResiduals(reco::TrackResiduals());
  }
  return trackExtra;
}
