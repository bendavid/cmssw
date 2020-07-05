#ifndef FWCore_Framework_ThinningProducerBase_h
#define FWCore_Framework_ThinningProducerBase_h

/** \class edm::ThinningProducerBase
\author W. David Dagenhart, created 11 June 2014
*/

#include "FWCore/Framework/interface/stream/EDProducer.h"

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

  template <typename Collection>
  class ThinningProducerBase : public stream::EDProducer<> {
  public:
    explicit ThinningProducerBase(ParameterSet const& pset);
    ~ThinningProducerBase() override;

    void registerThinnedAssociations(ProductRegistry const& productRegistry,
                                     ThinnedAssociationsHelper& thinnedAssociationsHelper) override;

  protected:
    edm::EDGetTokenT<Collection> inputToken_;
    edm::InputTag inputTag_;
    edm::EDPutTokenT<Collection> outputToken_;
    edm::EDPutTokenT<ThinnedAssociation> thinnedOutToken_;
  };

  template <typename Collection>
  ThinningProducerBase<Collection>::ThinningProducerBase(ParameterSet const& pset) {
    inputTag_ = pset.getParameter<InputTag>("inputTag");
    inputToken_ = consumes<Collection>(inputTag_);

    outputToken_ = produces<Collection>();
    thinnedOutToken_ = produces<ThinnedAssociation>();
  }

  template <typename Collection>
  ThinningProducerBase<Collection>::~ThinningProducerBase() {}

  template <typename Collection>
  void ThinningProducerBase<Collection>::registerThinnedAssociations(
      ProductRegistry const& productRegistry, ThinnedAssociationsHelper& thinnedAssociationsHelper) {
    BranchID associationID;
    BranchID thinnedCollectionID;

    // If the InputTag does not specify the process name, it is
    // possible that there will be more than one match found below.
    // For a particular event only one match is correct and the
    // others will be false. It even possible for some events one
    // match is correct and for others another is correct. This is
    // a side effect of the lookup mechanisms when the process name
    // is not specified.
    // When using the registry this generates one would have to
    // check the ProductIDs in ThinnedAssociation product to get
    // the correct association. This ambiguity will probably be
    // rare and possibly never occur in practice.
    std::vector<BranchID> parentCollectionIDs;

    ProductRegistry::ProductList const& productList = productRegistry.productList();
    for (auto const& product : productList) {
      BranchDescription const& desc = product.second;
      if (desc.unwrappedType().typeInfo() == typeid(Collection)) {
        if (desc.produced() && desc.moduleLabel() == moduleDescription().moduleLabel() &&
            desc.productInstanceName().empty()) {
          thinnedCollectionID = desc.branchID();
        }
        if (desc.moduleLabel() == inputTag_.label() && desc.productInstanceName() == inputTag_.instance()) {
          if (inputTag_.willSkipCurrentProcess()) {
            if (!desc.produced()) {
              parentCollectionIDs.push_back(desc.branchID());
            }
          } else if (inputTag_.process().empty() || inputTag_.process() == desc.processName()) {
            if (desc.produced()) {
              parentCollectionIDs.push_back(desc.originalBranchID());
            } else {
              parentCollectionIDs.push_back(desc.branchID());
            }
          }
        }
      }
      if (desc.produced() && desc.unwrappedType().typeInfo() == typeid(ThinnedAssociation) &&
          desc.moduleLabel() == moduleDescription().moduleLabel() && desc.productInstanceName().empty()) {
        associationID = desc.branchID();
      }
    }
    if (parentCollectionIDs.empty()) {
      // This could happen if the input collection was dropped. Go ahead and add
      // an entry and let the exception be thrown only if the module is run (when
      // it cannot find the product).
      thinnedAssociationsHelper.addAssociation(BranchID(), associationID, thinnedCollectionID);
    } else {
      for (auto const& parentCollectionID : parentCollectionIDs) {
        thinnedAssociationsHelper.addAssociation(parentCollectionID, associationID, thinnedCollectionID);
      }
    }
  }
}  // namespace edm
#endif
