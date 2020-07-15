#ifndef TrackProducer_SiPixelClusterThinningProducer_h
#define TrackProducer_SiPixelClusterThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "DataFormats/Common/interface/ThinnedRefSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

class SiPixelCluster;

namespace edm {
  class ParameterSetDescription;
}

class SiPixelClusterSelector {
public:
  SiPixelClusterSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  void preChoose(edm::Handle<edmNew::DetSetVector<SiPixelCluster>> clusters,
                 edm::Event const& event,
                 edm::EventSetup const& es);

  bool choose(unsigned int iIndex, SiPixelCluster const& iItem) const;

  void reset() { keysToSave_.clear(); }

private:
  std::vector<edm::EDGetTokenT<TrackingRecHitCollection>> trackingRecHitsTokens_;
  edm::ThinnedRefSet<edmNew::DetSetVector<SiPixelCluster>> keysToSave_;
};

using SiPixelClusterThinningProducer =
    edm::ThinningProducer<edmNew::DetSetVector<SiPixelCluster>, SiPixelClusterSelector>;

#endif
