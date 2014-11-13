#ifndef PromptLeptonEnrichmentHepMCFilter_h
#define PromptLeptonEnrichmentHepMCFilter_h
// -*- C++ -*-
//
// Package:    PromptLeptonEnrichmentHepMCFilter
// Class:      PromptLeptonEnrichmentHepMCFilter
//
/**\class PromptLeptonEnrichmentHepMCFilter PromptLeptonEnrichmentHepMCFilter.cc

Description: Filter events requiring N prompt leptons (i.e. not from hadronic showers)

*/
//
// Original Author:  Giovanni Pertrucciani (CERN)
//         Created:  11 Nov 2014
//
//


// system include files
#include <memory>

// user include files
#include "GeneratorInterface/Core/interface/BaseHepMCFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <HepMC/GenParticle.h>


//
// class decleration
//

class PromptLeptonEnrichmentHepMCFilter : public BaseHepMCFilter {
  public:
     PromptLeptonEnrichmentHepMCFilter(const edm::ParameterSet&);
    ~PromptLeptonEnrichmentHepMCFilter();

    virtual bool filter(const HepMC::GenEvent* evt);
    
  private:
    // ----------memeber function----------------------
    bool isPrompt(const HepMC::GenParticle *p) const ;

    // ----------member data ---------------------------

    std::vector<int> particleIDs;
    bool chargeconju;
    bool maybeprompt;
    int minnparticles;
    int maxnparticles;
    double minptcut;
    double maxptcut;
    double minetacut;
    double maxetacut;
    double masscut;
    };
#endif
