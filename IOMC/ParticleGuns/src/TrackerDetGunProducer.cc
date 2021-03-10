#include <ostream>

#include "IOMC/ParticleGuns/interface/TrackerDetGunProducer.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "CLHEP/Random/RandFlat.h"

using namespace edm;
using namespace std;

TrackerDetGunProducer::TrackerDetGunProducer(const ParameterSet& pset) : 
   BaseFlatGunProducer(pset)
{


   ParameterSet defpset ;
   ParameterSet pgun_params = 
      pset.getParameter<ParameterSet>("PGunParameters") ;
  
  fProbParticle_ = pgun_params.getParameter<std::vector<double> >("ProbParts");
  if (fProbParticle_.size() != fPartIDs.size()) throw cms::Exception("Configuration") << "Not all probabilities given for all particle types " << fProbParticle_.size() << ":" << fPartIDs.size() << " need them to match\n";
  for (unsigned int k=1; k<fProbParticle_.size(); ++k)
    fProbParticle_[k] += fProbParticle_[k-1];
  for (unsigned int k=0; k<fProbParticle_.size(); ++k)
    fProbParticle_[k] /= fProbParticle_[fProbParticle_.size()-1];
  
   fMinPt = pgun_params.getParameter<double>("MinPt");
   fMaxPt = pgun_params.getParameter<double>("MaxPt");
   fMinX = pgun_params.getParameter<double>("MinX");
   fMaxX = pgun_params.getParameter<double>("MaxX");
   fMinY = pgun_params.getParameter<double>("MinY");
   fMaxY = pgun_params.getParameter<double>("MaxY");
  
   detId_ = DetId(pgun_params.getParameter<unsigned int>("detId"));
  
  produces<HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
}

TrackerDetGunProducer::~TrackerDetGunProducer()
{
   // no need to cleanup GenEvent memory - done in HepMCProduct
}

void TrackerDetGunProducer::produce(Event &e, const EventSetup& es) 
{
   edm::Service<edm::RandomNumberGenerator> rng;
   CLHEP::HepRandomEngine* engine = &rng->getEngine(e.streamID());

  edm::ESHandle<GlobalTrackingGeometry> globalGeometry;
  es.get<GlobalTrackingGeometryRecord>().get(globalGeometry);
  
  const GeomDet* det = globalGeometry->idToDet(detId_);
  
  const double xlocal = fMaxX > fMinX ? CLHEP::RandFlat::shoot(engine, fMinX, fMaxX) : fMinX;
  const double ylocal = fMaxY > fMinY ? CLHEP::RandFlat::shoot(engine, fMinY, fMaxY) : fMinY;
  
  const GlobalPoint pos = det->surface().toGlobal(LocalPoint(xlocal, ylocal));
//   const GlobalPoint pos = det->surface().position();
  const GlobalVector posv(pos.x(), pos.y(), pos.z());
   
   // event loop (well, another step in it...)
  // no need to clean up GenEvent memory - done in HepMCProduct
  // here re-create fEvt (memory)
  //
  fEvt = new HepMC::GenEvent() ;
   
  // now actualy, cook up the event from PDGTable and gun parameters
  //

  // 1st, primary vertex
  //
  HepMC::GenVertex* Vtx = new HepMC::GenVertex( HepMC::FourVector(0.,0.,0.));
  
  // loop over particles
  //
  int    barcode(0), PartID(fPartIDs[0]);
  double r1      = CLHEP::RandFlat::shoot(engine, 0., 1.);
  for (unsigned int ip=0; ip<fPartIDs.size(); ip++) {
    if (r1 <= fProbParticle_[ip]) {
      PartID = fPartIDs[ip];
      break;
    }
  }
#ifdef DebugLog
  if (fVerbosity > 0) std::cout << "Random " << r1 << " PartID " << PartID
				<< std::endl;
#endif
  double pt     = fMinPt;
  if (fMaxPt > fMinPt) {
    pt = CLHEP::RandFlat::shoot(engine, fMinPt, fMaxPt);
  }
  double theta  = posv.theta();
  double phi    = posv.phi();
  const HepPDT::ParticleData* 
    PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
  double mass   = PData->mass().value() ;
  double mom    = pt/sin(theta) ;
  double px     = pt*cos(phi) ;
  double py     = pt*sin(phi) ;
  double pz     = mom*cos(theta) ;
  double energy = sqrt(mom*mom + mass*mass);

  HepMC::FourVector p(px,py,pz,energy) ;
  HepMC::GenParticle* Part = new HepMC::GenParticle(p,PartID,1);
  barcode++ ;
  Part->suggest_barcode(barcode);
  Vtx->add_particle_out(Part);

  if (fAddAntiParticle) {
    HepMC::FourVector ap(-px,-py,-pz,energy) ;
    int APartID = (PartID == 22 || PartID == 23) ? PartID : -PartID;
    HepMC::GenParticle* APart = new HepMC::GenParticle(ap,APartID,1);
    barcode++ ;
    APart->suggest_barcode(barcode);
    Vtx->add_particle_out(APart);
  }

  fEvt->add_vertex(Vtx) ;
  fEvt->set_event_number(e.id().event()) ;
  fEvt->set_signal_process_id(20) ;  
   
#ifdef DebugLog 
  if (fVerbosity > 0) fEvt->print() ;  
#endif

  std::unique_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
  BProduct->addHepMCData( fEvt );
  e.put(std::move(BProduct), "unsmeared");

  std::unique_ptr<GenEventInfoProduct> genEventInfo(new GenEventInfoProduct(fEvt));
  e.put(std::move(genEventInfo));
#ifdef DebugLog    
   if (fVerbosity > 0) 
     std::cout << " FlatRandomMultiParticlePtGunProducer : Event Generation Done " << std::endl;
#endif
}
//#include "FWCore/Framework/interface/MakerMacros.h"
//DEFINE_FWK_MODULE(TrackerDetGunProducer);
