#ifndef MuonIdentification_MuonTrackingRecHitThinningProducer_h
#define MuonIdentification_MuonTrackingRecHitThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "DataFormats/Common/interface/ThinnedRefSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace edm {
  class ParameterSetDescription;
}

class MuonTrackingRecHitSelector {
public:
  MuonTrackingRecHitSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  void preChoose(edm::Handle<TrackingRecHitCollection> hits, edm::Event const& event, edm::EventSetup const& es);

  bool choose(unsigned int iIndex, TrackingRecHit const& iItem) const;

  void reset() { keysToSave_.clear(); }

private:
  std::string cut_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  StringCutObjectSelector<reco::Muon> selector_;
  edm::ThinnedRefSet<TrackingRecHitCollection> keysToSave_;
};

using MuonTrackingRecHitThinningProducer = edm::ThinningProducer<TrackingRecHitCollection, MuonTrackingRecHitSelector>;

#endif
