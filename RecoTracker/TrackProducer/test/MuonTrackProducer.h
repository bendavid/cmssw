#ifndef RecoTracker_TrackProducer_MuonTrackProducer_H
#define RecoTracker_TrackProducer_MuonTrackProducer_H

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/Common/interface/Association.h"

class SiPixelCluster;
class SiStripCluster;

namespace pat {
  class Muon;
}

class MuonTrackProducer : public edm::stream::EDProducer<> {
public:
  MuonTrackProducer(const edm::ParameterSet&);

  ~MuonTrackProducer() override {}
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
  edm::EDPutTokenT<reco::TrackCollection> trackOutToken_;
};
#endif
