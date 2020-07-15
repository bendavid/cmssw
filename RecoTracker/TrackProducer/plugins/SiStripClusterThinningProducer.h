#ifndef TrackProducer_SiStripClusterThinningProducer_h
#define TrackProducer_SiStripClusterThinningProducer_h

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/stream/ThinningProducer.h"
#include "DataFormats/Common/interface/ThinnedRefSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterfwd.h"

namespace edm {
  class ParameterSetDescription;
}

class SiStripClusterSelector {
public:
  SiStripClusterSelector(edm::ParameterSet const& pset, edm::ConsumesCollector&& cc);

  static void fillPSetDescription(edm::ParameterSetDescription& desc);

  void preChoose(edm::Handle<edmNew::DetSetVector<SiStripCluster>> clusters,
                 edm::Event const& event,
                 edm::EventSetup const& es);

  bool choose(unsigned int iIndex, SiStripCluster const& iItem) const;

  void reset() { keysToSave_.clear(); }

private:
  std::vector<edm::EDGetTokenT<TrackingRecHitCollection>> trackingRecHitsTokens_;
  edm::ThinnedRefSet<edmNew::DetSetVector<SiStripCluster>> keysToSave_;
};

using SiStripClusterThinningProducer =
    edm::ThinningProducer<edmNew::DetSetVector<SiStripCluster>, SiStripClusterSelector>;

#endif
