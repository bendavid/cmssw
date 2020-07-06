#ifndef RecoMuon_MuonIdentification_ThinnedTrackExtraProducer_H
#define RecoMuon_MuonIdentification_ThinnedTrackExtraProducer_H

/** \class ThinnedTrackExtraProducer
 *
 *  \author J. Bendavid
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "FWCore/Framework/interface/stream/ThinningProducerBase.h"

typedef edm::ThinningProducerBase<reco::TrackExtraCollection> ThinnedTrackExtraProducerBase;

class ThinnedTrackExtraProducer : public ThinnedTrackExtraProducerBase {
public:
  /// Constructor
  ThinnedTrackExtraProducer(const edm::ParameterSet&);

  /// Destructor
  ~ThinnedTrackExtraProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  std::string cut_;
  bool storeHits_;
  bool isInputThinned_;
  bool skipMissingInput_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
  edm::EDPutTokenT<TrackingRecHitCollection> hitsToken_;
  StringCutObjectSelector<reco::Muon> selector_;
};
#endif
