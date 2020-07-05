#ifndef FWCore_Framework_ThinningProducer_h
#define FWCore_Framework_ThinningProducer_h

/** \class edm::ThinningProducer
\author W. David Dagenhart, created 11 June 2014
*/

#include "FWCore/Framework/interface/stream/ThinningProducerBase.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/Common/interface/ThinnedAssociation.h"
#include "DataFormats/Provenance/interface/ProductRegistry.h"
#include "DataFormats/Provenance/interface/ThinnedAssociationsHelper.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/propagate_const.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <memory>

namespace edm {

  class EventSetup;

  template <typename Collection, typename Selector>
  class ThinningProducer : public ThinningProducerBase<Collection> {
  public:
    explicit ThinningProducer(ParameterSet const& pset);
    ~ThinningProducer() override;

    static void fillDescriptions(ConfigurationDescriptions& descriptions);

    void produce(Event& event, EventSetup const& eventSetup) override;

  private:
    edm::propagate_const<std::unique_ptr<Selector>> selector_;
  };

  template <typename Collection, typename Selector>
  ThinningProducer<Collection, Selector>::ThinningProducer(ParameterSet const& pset)
      : ThinningProducerBase<Collection>(pset),
        selector_(new Selector(pset, ThinningProducerBase<Collection>::consumesCollector())) {}

  template <typename Collection, typename Selector>
  ThinningProducer<Collection, Selector>::~ThinningProducer() {}

  template <typename Collection, typename Selector>
  void ThinningProducer<Collection, Selector>::fillDescriptions(ConfigurationDescriptions& descriptions) {
    ParameterSetDescription desc;
    desc.setComment("Produces thinned collections and associations to them");
    desc.add<edm::InputTag>("inputTag");
    Selector::fillDescription(desc);
    descriptions.addDefault(desc);
  }

  template <typename Collection, typename Selector>
  void ThinningProducer<Collection, Selector>::produce(Event& event, EventSetup const& eventSetup) {
    auto inputCollection = event.getHandle(ThinningProducerBase<Collection>::inputToken_);

    edm::Event const& constEvent = event;
    selector_->preChoose(inputCollection, constEvent, eventSetup);

    Collection thinnedCollection;
    ThinnedAssociation thinnedAssociation;

    unsigned int iIndex = 0;
    for (auto iter = inputCollection->begin(), iterEnd = inputCollection->end(); iter != iterEnd; ++iter, ++iIndex) {
      if (selector_->choose(iIndex, *iter)) {
        thinnedCollection.push_back(*iter);
        thinnedAssociation.push_back(iIndex);
      }
    }
    OrphanHandle<Collection> orphanHandle =
        event.emplace(ThinningProducerBase<Collection>::outputToken_, std::move(thinnedCollection));

    thinnedAssociation.setParentCollectionID(inputCollection.id());
    thinnedAssociation.setThinnedCollectionID(orphanHandle.id());
    event.emplace(ThinningProducerBase<Collection>::thinnedOutToken_, std::move(thinnedAssociation));
  }
}  // namespace edm
#endif
