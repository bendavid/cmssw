/** \class ThinnedTrackExtraProducer
 *  See header file.
 *
 *  \author J. Bendavid
 */

#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoMuon/MuonIdentification/plugins/ThinnedTrackExtraProducer.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/ThinnedAssociation.h"
#include "DataFormats/Provenance/interface/ProductID.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

/// Constructor
ThinnedTrackExtraProducer::ThinnedTrackExtraProducer(const edm::ParameterSet& pset)
    : ThinnedTrackExtraProducerBase(pset),
      cut_(pset.getParameter<std::string>("cut")),
      storeHits_(pset.getParameter<bool>("storeHits")),
      isInputThinned_(pset.getParameter<bool>("isInputThinned")),
      skipMissingInput_(pset.getParameter<bool>("skipMissingInput")),
      muonToken_(consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonTag"))),
      selector_(cut_) {
  if (storeHits_) {
    hitsToken_ = produces<TrackingRecHitCollection>();
  }
}

/// Destructor
ThinnedTrackExtraProducer::~ThinnedTrackExtraProducer() {}

void ThinnedTrackExtraProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Produces thinned collection of TrackExtras and TrackingRecHits");
  desc.add<edm::InputTag>("inputTag", edm::InputTag("generalTracks"));
  desc.add<std::string>("cut", "pt > 3. || isPFMuon");
  desc.add<bool>("storeHits", true);
  desc.add<bool>("isInputThinned", false);
  desc.add<bool>("skipMissingInput", false);
  desc.add<edm::InputTag>("muonTag", edm::InputTag("muons"));
  descriptions.add("thinnedGeneralTracks", desc);
}

void ThinnedTrackExtraProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  auto inputCollection = event.getHandle(ThinnedTrackExtraProducerBase::inputToken_);
  auto muonCollection = event.getHandle(muonToken_);

  reco::TrackExtraCollection thinnedCollection;
  edm::ThinnedAssociation thinnedAssociation;
  TrackingRecHitCollection hitCollection;

  // build set of TrackExtras to keep (the ones associated with bestTrack for muons passing the
  // selection criteria)
  std::set<std::pair<edm::ProductID, edm::Ref<reco::TrackExtraCollection>::key_type> > trackExtraRefSet;
  std::set<std::array<double, 6> > trackExtraValSet;
  for (const auto& muon : *muonCollection) {
    if (!selector_(muon)) {
      continue;
    }
    const edm::Ref<reco::TrackExtraCollection>& trackExtraRef = muon.bestTrack()->extra();
    if (isInputThinned_) {
      // if the input collection is already thinned, references won't match, so need to match by value
      // The track extra is not guaranteed to be available, so check
      if (trackExtraRef.isAvailable()) {
        const reco::TrackExtra* trackExtra = muon.bestTrack()->extra().get();
        std::array<double, 6> parms = {{trackExtra->innerPosition().x(),
                                        trackExtra->innerPosition().y(),
                                        trackExtra->innerPosition().z(),
                                        trackExtra->innerMomentum().x(),
                                        trackExtra->innerMomentum().y(),
                                        trackExtra->innerMomentum().z()}};
        trackExtraValSet.insert(parms);
      }
    } else {
      //if the input collection is the original, match by reference since this should be faster
      trackExtraRefSet.emplace(trackExtraRef.id(), trackExtraRef.key());
    }
  }

  // Build thinned collection for selected TrackExtras and fill additional output collection
  // of TrackingRecHits if necessary
  if (!skipMissingInput_ || inputCollection.isValid()) {
    TrackingRecHitRefProd rHits = storeHits_ ? event.getRefBeforePut(hitsToken_) : TrackingRecHitRefProd();
    for (auto it = inputCollection->begin(); it != inputCollection->end(); ++it) {
      unsigned int idx = it - inputCollection->begin();
      bool keep = false;
      if (isInputThinned_) {
        std::array<double, 6> parms = {{it->innerPosition().x(),
                                        it->innerPosition().y(),
                                        it->innerPosition().z(),
                                        it->innerMomentum().x(),
                                        it->innerMomentum().y(),
                                        it->innerMomentum().z()}};
        if (trackExtraValSet.count(parms)) {
          keep = true;
        }
      } else {
        if (trackExtraRefSet.count(std::make_pair(inputCollection.id(), idx))) {
          keep = true;
        }
      }
      if (keep) {
        // use explicit TrackExtra constructor here because it is otherwise
        // not possible to reset the hits after the fact
        thinnedCollection.emplace_back(it->outerPosition(),
                                       it->outerMomentum(),
                                       it->outerOk(),
                                       it->innerPosition(),
                                       it->innerMomentum(),
                                       it->innerOk(),
                                       it->outerStateCovariance(),
                                       it->outerDetId(),
                                       it->innerStateCovariance(),
                                       it->innerDetId(),
                                       it->seedDirection(),
                                       it->seedRef());
        thinnedAssociation.push_back(idx);

        if (storeHits_) {
          thinnedCollection.back().setHits(rHits, hitCollection.size(), it->recHitsSize());
          for (const auto& hit : it->recHits()) {
            hitCollection.push_back(hit->clone());
          }
        }
      }
    }
    thinnedAssociation.setParentCollectionID(inputCollection.id());
  }

  if (storeHits_) {
    event.emplace(hitsToken_, std::move(hitCollection));
  }

  edm::OrphanHandle<reco::TrackExtraCollection> orphanHandle =
      event.emplace(ThinnedTrackExtraProducerBase::outputToken_, std::move(thinnedCollection));

  thinnedAssociation.setThinnedCollectionID(orphanHandle.id());
  event.emplace(ThinnedTrackExtraProducerBase::thinnedOutToken_, std::move(thinnedAssociation));
}
