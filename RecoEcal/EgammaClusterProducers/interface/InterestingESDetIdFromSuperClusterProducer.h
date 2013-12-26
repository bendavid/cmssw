#ifndef _INTERESTINGESDETIDFROMSUPERCLUSTERPRODUCER_H
#define _INTERESTINGESDETIDFROMSUPERCLUSTERPRODUCER_H

// -*- C++ -*-
//
// Package:    InterestingESDetIdFromSuperClusterProducer
// Class:      InterestingESDetIdFromSuperClusterProducer
// 
/**\class InterestingESDetIdFromSuperClusterProducer 
Adapted from InterestingDetIdCollectionProducer by J.Bendavid
 
Make a collection of ES detids to be kept typically in a AOD rechit collection

The following classes of "interesting id" are considered

    1.All rechits included in all ES subclusters within the SuperCluster
*/



// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

class CaloTopology;
class EcalSeverityLevelAlgo;

class InterestingESDetIdFromSuperClusterProducer : public edm::EDProducer {
   public:
      //! ctor
      explicit InterestingESDetIdFromSuperClusterProducer(const edm::ParameterSet&);
      virtual void beginRun (edm::Run const&, const edm::EventSetup&) override final;
      //! producer
      virtual void produce(edm::Event &, const edm::EventSetup&);

   private:
      // ----------member data ---------------------------
      edm::EDGetTokenT<EcalRecHitCollection>         recHitsToken_;
      edm::EDGetTokenT<reco::SuperClusterCollection> superClustersToken_;
      std::string interestingDetIdCollection_;

};

#endif
