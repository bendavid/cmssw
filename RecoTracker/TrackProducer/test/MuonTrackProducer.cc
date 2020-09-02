#include "RecoTracker/TrackProducer/test/MuonTrackProducer.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "FWCore/Framework/interface/MakerMacros.h"

MuonTrackProducer::MuonTrackProducer(const edm::ParameterSet& pset)
 :  muonToken_(consumes<std::vector<pat::Muon> >(pset.getParameter<edm::InputTag>("muonTag"))),
    trackOutToken_(produces<reco::TrackCollection>()) {}

void MuonTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Reproduces tracks from those embedded in PatMuon to test track refitter.");
  desc.add<edm::InputTag>("muonTag", edm::InputTag("slimmedMuons"));
  descriptions.addWithDefaultLabel(desc);
}

void MuonTrackProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  
  auto muons = event.getHandle(muonToken_);
  std::vector<reco::Track> tracksOut;

  //loop over muons and store tracks
  for (auto const& muon : *muons) {
//     if (muon.pt()>4.5 && !muon.bestTrack()->extra().isAvailable()) {
//       std::cout << "missing trackextra which is supposed to be here." << std::endl;
//       assert(false);
//     }
//     
//     if (!muon.innerTrack().isAvailable()) {
//       continue;
//     }
//     const reco::Track& track = *muon.innerTrack();
    const reco::Track& track = *muon.bestTrack();
    std::cout << "found muon track" << std::endl;
    if (track.extra().isAvailable()) {
      std::cout << "trackExtra available" << std::endl;
      std::cout << track.p() << " " << std::sqrt(track.extra()->innerMomentum().mag2()) << std::endl;
      tracksOut.emplace_back(track);
    }
    else {
      std::cout << "trackExtra not available" << std::endl;
    }
  }
  
  event.emplace(trackOutToken_, std::move(tracksOut));  
}

DEFINE_FWK_MODULE(MuonTrackProducer);
