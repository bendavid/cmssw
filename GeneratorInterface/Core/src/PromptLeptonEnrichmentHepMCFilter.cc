
#include "GeneratorInterface/Core/interface/PromptLeptonEnrichmentHepMCFilter.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <iostream>

using namespace edm;
using namespace std;


PromptLeptonEnrichmentHepMCFilter::PromptLeptonEnrichmentHepMCFilter(const edm::ParameterSet& iConfig) :
        particleIDs(iConfig.getParameter<std::vector<int> >("ParticleIDs")),
        chargeconju(iConfig.getParameter<bool>("ChargeConjugation")),
        maybeprompt(iConfig.getParameter<bool>("TreatDubiousAsPrompt")),
        minnparticles(iConfig.getParameter<int>("MinNumberOfParticles")),
        maxnparticles(iConfig.getParameter<int>("MaxNumberOfParticles")),
        minptcut(iConfig.getParameter<double>("MinPt")),
        maxptcut(iConfig.getParameter<double>("MaxPt")),
        minetacut(iConfig.getParameter<double>("MinEta")),
        maxetacut(iConfig.getParameter<double>("MaxEta")),
        masscut(iConfig.getParameter<double>("MinPairMass"))
{

}


PromptLeptonEnrichmentHepMCFilter::~PromptLeptonEnrichmentHepMCFilter()
{

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool PromptLeptonEnrichmentHepMCFilter::filter(const HepMC::GenEvent* evt)
{
  //printf("PromptLeptonEnrichmentHepMCFilter::filter!\n");
  int count = 0; 
  std::vector<math::XYZTLorentzVectorD> p4s; 
  for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
                  p != evt->particles_end(); ++p ) {

          int partid = (*p)->pdg_id();
          if (chargeconju) partid = abs(partid);
          if (std::find(particleIDs.begin(), particleIDs.end(), partid) == particleIDs.end()) continue;

          const auto p4 = math::XYZTLorentzVectorD((*p)->momentum());
          if(p4.pt()  >  minptcut  && p4.pt() <  maxptcut  &&
             p4.eta() >  minetacut && p4.eta()  <  maxetacut ) {
        
             //printf("particle pdgId %d status %d \n", (*p)->pdg_id(), (*p)->status());
             if (isPrompt(*p)) {
                 bool dup = false;
                 for (const auto &otherp4 : p4s) {
                    if ((p4+otherp4).mass() < masscut) {
                        dup = true;
                        break;
                    }
                 }
                 if (!dup) {
                    p4s.push_back(p4);
                    ++count;
                 }
             }
          }

  }
  //printf("In total found %d prompt leptons passing selection: %d\n", count, (count >= minnparticles && (maxnparticles == -1 || count <= maxnparticles)));
  return (count >= minnparticles && (maxnparticles == -1 || count <= maxnparticles));
}

bool PromptLeptonEnrichmentHepMCFilter::isPrompt(const HepMC::GenParticle *p) const
{
    
    if ( !p->production_vertex() ) return maybeprompt;

    for ( HepMC::GenVertex::particle_iterator 
            des=p->production_vertex()->particles_begin(HepMC::parents);
            des != p->production_vertex()->particles_end(HepMC::parents);
            ++des ) {
        //printf("particle pdgId %d status %d from mother pdgId %d status %d\n", p->pdg_id(), p->status(), (*des)->pdg_id(), (*des)->status());
        int id = abs((*des)->pdg_id()), stat = (*des)->status();
        if ((11 <= id && id <= 15) && isPrompt(*des)) return true; // Lepton
        if (id >= 1000001 || (22 <= id && id <= 39)) return true; // BSM or EWK
        if (id >= 100) return false;                        // Hadronic
        if ((id <= 5 || id == 21) && (stat == 2)) return false; // Shower
        // boh?
    }
    return maybeprompt;
}
