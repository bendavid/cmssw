#ifndef MuonIdentification_MuonTrackExtraThinningProducer_h
#define MuonIdentification_MuonTrackExtraThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "DataFormats/Common/interface/ThinnedRefSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace edm {
  class ParameterSetDescription;
}

class MuonTrackExtraSelector {
public:
  MuonTrackExtraSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  void preChoose(edm::Handle<reco::TrackExtraCollection> trackExtras,
                 edm::Event const& event,
                 edm::EventSetup const& es);

  std::optional<reco::TrackExtra> choose(unsigned int iIndex, reco::TrackExtra const& iItem) const;

  void reset() { keysToSave_.clear(); }

private:
  std::string cut_;
  bool slimTrajParams_;
  bool slimResiduals_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  StringCutObjectSelector<reco::Muon> selector_;
  edm::ThinnedRefSet<reco::TrackExtraCollection> keysToSave_;
};

using MuonTrackExtraThinningProducer = edm::ThinningProducer<reco::TrackExtraCollection, MuonTrackExtraSelector>;

#endif
