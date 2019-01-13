// -*- C++ -*-
//
// Package:    AnalysisMTD/BackProp
// Class:      BackProp
//
/**\class BackProp BackProp.cc AnalysisMTD/BackProp/plugins/BackProp.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue, 13 Nov 2018 09:00:35 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TH1D.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class BackProp : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit BackProp(const edm::ParameterSet&);
      ~BackProp();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::ValueMap<float> > simt0Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > simsigmat0Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > trkt0Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > trksigmat0Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > trkpidt0Token_;
      edm::EDGetTokenT<edm::ValueMap<float> > trkpidsigmat0Token_;
      edm::EDGetTokenT<float> gent0Token_;
      edm::EDGetTokenT<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genxyz0Token_;
      edm::EDGetTokenT<reco::VertexCollection> vtxsToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxspidToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxsfastToken_;
      edm::EDGetTokenT<edm::PSimHitContainer> btlSimHitToken_;
      edm::EDGetTokenT<edm::PSimHitContainer> etlSimHitToken_;
      
/*      TH1D *hdttrackbtl_;
      TH1D *hdttracketl_;
      TH1D *hdthitbtl_;
      TH1D *hdthitetl_;
      TH1D *hdtdtbtl_;
      TH1D *hdtdtetl_;
      TH1D *hdxbtl_;
      TH1D *hdxetl_;
      TH1D *hdybtl_;
      TH1D *hdyetl_;
      TH1D *hresxbtl_;
      TH1D *hresxetl_;
      TH1D *hresybtl_;
      TH1D *hresyetl_;      
      TH1D *hsimresxbtl_;
      TH1D *hsimresxetl_;
      TH1D *hsimresybtl_;
      TH1D *hsimresyetl_;    */  

    TH1D *hdtvtxgen_;
    TH1D *hdtvtxpidgen_;
    TH1D *hdtvtxfastgen_;
    TH1D *hdttrkvtx_;
    TH1D *hdttrkpidvtx_;
    TH1D *hdttrkgen_;
    TH1D *hdttrkpidgen_;
    TH1D *hdttrksim_;
    TH1D *hdttrkpidsim_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BackProp::BackProp(const edm::ParameterSet& iConfig)
 :
//   tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
  tracksToken_(consumes<TrackCollection>(edm::InputTag("generalTracks"))),
  simt0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("trackTimeValueMapProducer","generalTracksPerfectResolutionModel"))),
  simsigmat0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("trackTimeValueMapProducer","generalTracksPerfectResolutionModelResolution"))),
  trkt0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD","generalTrackt0"))),
  trksigmat0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("trackExtenderWithMTD","generalTracksigmat0"))),
  trkpidt0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("TOFPIDProducer","t0"))),
  trkpidsigmat0Token_(consumes<edm::ValueMap<float> >(edm::InputTag("TOFPIDProducer","sigmat0"))),
  gent0Token_(consumes<float>(edm::InputTag("genParticles","t0"))),
  genxyz0Token_(consumes<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >(edm::InputTag("genParticles","xyz0"))),
  vtxsToken_(consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices4DWithBS"))),
  vtxspidToken_(consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices4DwithPIDWithBS"))),
  vtxsfastToken_(consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices4DfastsimWithBS"))),
  btlSimHitToken_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","FastTimerHitsBarrel"))),
  etlSimHitToken_(consumes<edm::PSimHitContainer>(edm::InputTag("g4SimHits","FastTimerHitsEndcap")))
{
   //now do what ever initialization is needed

  edm::Service<TFileService> fs;
  
  hdtvtxgen_ = fs->make<TH1D>("hdtvtxgen","",50,-0.2,0.2);
  hdtvtxpidgen_ = fs->make<TH1D>("hdtvtxpidgen","",50,-0.2,0.2);
  hdtvtxfastgen_ = fs->make<TH1D>("hdtvtxfastgen","",50,-0.2,0.2);
  
  hdttrkvtx_ = fs->make<TH1D>("hdttrkvtx","",50,-0.2,0.2);
  hdttrkpidvtx_ = fs->make<TH1D>("hdttrkpidvtx","",50,-0.2,0.2);
  
  hdttrkgen_ = fs->make<TH1D>("hdttrkgen","",50,-0.2,0.2);
  hdttrkpidgen_ = fs->make<TH1D>("hdttrkpidgen","",50,-0.2,0.2);
  
  hdttrksim_ = fs->make<TH1D>("hdttrksim","",50,-0.2,0.2);
  hdttrkpidsim_ = fs->make<TH1D>("hdttrkpidsim","",50,-0.2,0.2);
  
//   hdttrackbtl_ = fs->make<TH1D>("hdttrackbtl","",50,-0.2,0.2);
//   hdttracketl_ = fs->make<TH1D>("hdttracketl","",50,-0.2,0.2);
// 
//   hdthitbtl_ = fs->make<TH1D>("hdthitbtl","",50,-0.2,0.2);
//   hdthitetl_ = fs->make<TH1D>("hdthitetl","",50,-0.2,0.2);
// 
//   hdtdtbtl_ = fs->make<TH1D>("hdtdtbtl","",50,-0.2,0.2);
//   hdtdtetl_ = fs->make<TH1D>("hdtdtetl","",50,-0.2,0.2);
// 
//   hdxbtl_ = fs->make<TH1D>("hdxbtl","",50,-2.,2.);
//   hdxetl_ = fs->make<TH1D>("hdxetl","",50,-5.,5.);
// 
//   hdybtl_ = fs->make<TH1D>("hdybtl","",50,-2.,2.);
//   hdyetl_ = fs->make<TH1D>("hdyetl","",50,-5.,5.);
//   
//   hresxbtl_ = fs->make<TH1D>("hresxbtl","",50,-2.,2.);
//   hresxetl_ = fs->make<TH1D>("hresxetl","",50,-2.,2.);
//   
//   hresybtl_ = fs->make<TH1D>("hresybtl","",50,-2.,2.);
//   hresyetl_ = fs->make<TH1D>("hresyetl","",50,-2.,2.);
//   
//   hsimresxbtl_ = fs->make<TH1D>("hsimresxbtl","",50,-2.,2.);
//   hsimresxetl_ = fs->make<TH1D>("hsimresxetl","",50,-2.,2.);
//   
//   hsimresybtl_ = fs->make<TH1D>("hsimresybtl","",50,-2.,2.);
//   hsimresyetl_ = fs->make<TH1D>("hsimresyetl","",50,-2.,2.);  
  
}


BackProp::~BackProp()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BackProp::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::ESHandle<MTDGeometry> geom;
  iSetup.get<MTDDigiGeometryRecord>().get(geom);   
  
//   edm::ESHandle<MTDDetLayerGeometry> geo;
//   iSetup.get<MTDRecoGeometryRecord>().get(geo);

  Handle<TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);
  
  Handle<edm::ValueMap<float> > hsimt0;
  iEvent.getByToken(simt0Token_, hsimt0);

  Handle<edm::ValueMap<float> > hsimsigmat0;
  iEvent.getByToken(simsigmat0Token_, hsimsigmat0);
  
  Handle<edm::ValueMap<float> > htrkt0;
  iEvent.getByToken(trkt0Token_, htrkt0);

  Handle<edm::ValueMap<float> > htrksigmat0;
  iEvent.getByToken(trksigmat0Token_, htrksigmat0);  
  
  Handle<edm::ValueMap<float> > htrkpidt0;
  iEvent.getByToken(trkpidt0Token_, htrkpidt0);

  Handle<edm::ValueMap<float> > htrkpidsigmat0;
  iEvent.getByToken(trkpidsigmat0Token_, htrkpidsigmat0);  
  
  Handle<float> gent0;
  iEvent.getByToken(gent0Token_, gent0);
  double tgen = *gent0.product();
  
  Handle<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > genxyz0;
  iEvent.getByToken(genxyz0Token_, genxyz0);
  
  Handle<reco::VertexCollection> vtxs;
  iEvent.getByToken(vtxsToken_, vtxs);

  Handle<reco::VertexCollection> vtxspid;
  iEvent.getByToken(vtxspidToken_, vtxspid);
  
  Handle<reco::VertexCollection> vtxsfast;
  iEvent.getByToken(vtxsfastToken_, vtxsfast);
  
//   Handle<edm::PSimHitContainer> btlSimHits;
//   iEvent.getByToken(btlSimHitToken_,btlSimHits);
// 
//   Handle<edm::PSimHitContainer> etlSimHits;
//   iEvent.getByToken(etlSimHitToken_,etlSimHits);    
  
//   if (vtxs->size()==0) {
//     return;
//   }
  
//   printf("filling vertex plots, %i, %i, %i\n",int(vtxs->size()), int(vtxspid->size()),int(vtxsfast->size()));
  
  hdtvtxgen_->Fill(vtxs->front().t()-tgen);
  hdtvtxpidgen_->Fill(vtxspid->front().t()-tgen);
  hdtvtxfastgen_->Fill(vtxsfast->front().t()-tgen);
  
//   printf("done filling vertex plots\n");
  
  const reco::Vertex &vtx = vtxspid->front();
  
//   double dzvtxgen = std::abs(vtx.position().z()-genxyz0->z());
//   if (dzvtxgen>0.1) {
//     return;
//   }
  
//   if (vtx.tError()>0.) {
//     double dtvtxgen = vtx.t() - tgen;
//     hdtvtxgen_->Fill(dtvtxgen);
//   }
    
    
  for (unsigned int itrk=0; itrk<tracks->size(); ++itrk) {
    const reco::Track &track = (*tracks)[itrk];

    const reco::TrackRef trackref(tracks,itrk);
    
    double simt0 = (*hsimt0)[trackref];
    double simsigmat0 = (*hsimsigmat0)[trackref];
    
    double trkt0 = (*htrkt0)[trackref];
    double trksigmat0 = (*htrksigmat0)[trackref];
    
    double trkpidt0 = (*htrkpidt0)[trackref];
    double trkpidsigmat0 = (*htrkpidsigmat0)[trackref];

    if (vtx.tError()>0. && vtx.tError()<0.05 && std::abs(track.dz(vtx.position()))<0.1) {
      if (trksigmat0>=0.) {
        hdttrkvtx_->Fill(trkt0 - vtx.t());
      }
      if (trkpidsigmat0>=0.) {
        hdttrkpidvtx_->Fill(trkpidt0 - vtx.t());
      }
      
    }
    else {
      continue;
    }
    
    if (simsigmat0>=0.) {
      if (trksigmat0>=0.) {
        double dttrkgen = trkt0 - tgen;
        hdttrkgen_->Fill(dttrkgen);
        
        double dttrksim = trkt0 - simt0;
        hdttrksim_->Fill(dttrksim);
      }
      if (trkpidsigmat0>=0.) {
        double dttrkpidgen = trkpidt0 - tgen;
        hdttrkpidgen_->Fill(dttrkpidgen);
        
        double dttrkpidsim = trkpidt0 - simt0;
        hdttrkpidsim_->Fill(dttrkpidsim);
      }
    }
    
  }
  
//   for (const reco::Track &track : *tracks.product()) {
//     double ttrack = track.t0();
//     double dttrack = ttrack - tgen;
//     double simsigmat0 = (*trksigmat0)
    

    
//   }
    
//     for(TrackCollection::const_iterator itTrack = tracks->begin();
//         itTrack != tracks->end();
//         ++itTrack) {
      // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
//     }

// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
// 
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
}


// ------------ method called once each job just before starting event loop  ------------
void
BackProp::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
BackProp::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BackProp::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BackProp);
