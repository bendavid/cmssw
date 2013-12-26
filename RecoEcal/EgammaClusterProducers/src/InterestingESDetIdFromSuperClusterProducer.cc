#include "RecoEcal/EgammaClusterProducers/interface/InterestingESDetIdFromSuperClusterProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "DataFormats/DetId/interface/DetIdCollection.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"


InterestingESDetIdFromSuperClusterProducer::InterestingESDetIdFromSuperClusterProducer(const edm::ParameterSet& iConfig) 
{

  recHitsToken_ = 
	  consumes<EcalRecHitCollection>(iConfig.getParameter< edm::InputTag > ("recHitsLabel"));
  superClustersToken_ = 
	  consumes<reco::SuperClusterCollection>(iConfig.getParameter< edm::InputTag >("superClustersLabel"));

  interestingDetIdCollection_ = iConfig.getParameter<std::string>("interestingDetIdCollection");

   //register your products
  produces< DetIdCollection > (interestingDetIdCollection_) ;

}


void InterestingESDetIdFromSuperClusterProducer::beginRun (edm::Run const& run, const edm::EventSetup & iSetup)  
{

}

// ------------ method called to produce the data  ------------
void
InterestingESDetIdFromSuperClusterProducer::produce (edm::Event& iEvent, 
                                const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   // take SuperClusters
  Handle<reco::SuperClusterCollection> pClusters;
  iEvent.getByToken(superClustersToken_, pClusters);
  
  // take EcalRecHits
  Handle<EcalRecHitCollection> recHitsHandle;
  iEvent.getByToken(recHitsToken_,recHitsHandle);

  //Create empty output collections
  std::vector<DetId> indexToStore;
  indexToStore.reserve(1000);

  reco::SuperClusterCollection::const_iterator sclusIt;

  std::vector<DetId> xtalsToStore;
  xtalsToStore.reserve(50);
  
  //loop over superclusters
  for (sclusIt=pClusters->begin(); sclusIt!=pClusters->end(); sclusIt++) {
    //loop over subclusters
    for (reco::CaloCluster_iterator clusIt = sclusIt->preshowerClustersBegin(); clusIt!=sclusIt->preshowerClustersEnd(); ++clusIt) {   
      std::vector<std::pair<DetId,float> > clusterDetIds = (*clusIt)->hitsAndFractions();
      std::vector<std::pair<DetId,float> >::iterator posCurrent;

      xtalsToStore.clear();
      const std::vector<std::pair<DetId,float > > &xtalsInClus=(*clusIt)->hitsAndFractions();
      
      for (unsigned int ii=0;ii<xtalsInClus.size();ii++)
        {
            xtalsToStore.push_back(xtalsInClus[ii].first);
        }
      
      indexToStore.insert(indexToStore.end(),xtalsToStore.begin(),xtalsToStore.end());
    }
  }

  //unify the vector
  std::sort(indexToStore.begin(),indexToStore.end());
  std::unique(indexToStore.begin(),indexToStore.end());
  
  std::auto_ptr< DetIdCollection > detIdCollection (new DetIdCollection(indexToStore) ) ;

 
  iEvent.put( detIdCollection, interestingDetIdCollection_ );

}
