/** \class EcalMaxSampleUncalibRecHitProducer
 *   produce ECAL uncalibrated rechits from dataframes 
 *
 *  \author G. Franzoni, E. Di Marco
 *
 */
#include "RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMaxSample.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

EcalUncalibRecHitWorkerMaxSample::EcalUncalibRecHitWorkerMaxSample(const edm::ParameterSet& ps, edm::ConsumesCollector& c) :
  EcalUncalibRecHitWorkerBaseClass( ps ,c)
{
}


void
EcalUncalibRecHitWorkerMaxSample::set(const edm::EventSetup& es)
{
  es.get<EcalGainRatiosRcd>().get(gains);
  es.get<EcalPedestalsRcd>().get(peds);
  
  weightseb.resize(10);
  weightseb[0] = 0.;
  weightseb[1] = 0.;
  weightseb[2] = 0.;
  weightseb[3] = 0.;
  weightseb[4] = 0.18463;
  weightseb[5] = 0.424762;
  weightseb[6] = 0.349301;
  weightseb[7] = 0.176883;
  weightseb[8] = 0.0129613;
  weightseb[9] = 0.;
  
  weightsee.resize(10);
  weightsee[0] = 0.;
  weightsee[1] = 0.;
  weightsee[2] = 0.;
  weightsee[3] = 0.;
  weightsee[4] = 0.204788;
  weightsee[5] = 0.415302;
  weightsee[6] = 0.339506;
  weightsee[7] = 0.172017;
  weightsee[8] = 0.0145885;
  weightsee[9] = 0.;  
  
  {
    std::vector<double> &pulse = pulseeb;
    double alpha = 1.138;
    double beta = 1.655;
    double alfabeta = alpha*beta;
    
    pulse.resize(10);
    for (int isample=0; isample<10; ++isample) {
      double dt = isample-5;
      if(dt > -alfabeta)  {
        double dtsbeta=dt/beta;
        double variable=1.+dt/alfabeta ;
        double puiss=pow(variable,alpha);
        double val = puiss*exp(-dtsbeta);
        pulse[isample] = val;
      }
      else {
        pulse[isample] = 0.;
      }
    }
  }
  
  {
    std::vector<double> &pulse = pulseee;
    double alpha = 1.890;
    double beta = 1.400;
    double alfabeta = alpha*beta;
    
    pulse.resize(10);
    for (int isample=0; isample<10; ++isample) {
      double dt = isample-5;
      if(dt > -alfabeta)  {
        double dtsbeta=dt/beta;
        double variable=1.+dt/alfabeta ;
        double puiss=pow(variable,alpha);
        double val = puiss*exp(-dtsbeta);
        pulse[isample] = val;
      }
      else {
        pulse[isample] = 0.;
      }
    }
  }  
  
}

bool
EcalUncalibRecHitWorkerMaxSample::run( const edm::Event & evt, 
                const EcalDigiCollection::const_iterator & itdg, 
                EcalUncalibratedRecHitCollection & result )
{
        DetId detid(itdg->id());

        if ( detid.subdetId() == EcalBarrel ) {
                const EcalPedestals::Item * aped = &peds->barrel(EBDetId(detid).hashedIndex());
                const EcalMGPAGainRatio * aGain  = &gains->barrel(EBDetId(detid).hashedIndex());
                result.push_back( ebAlgo_.makeRecHit(*itdg, aped, aGain, weightseb, pulseeb) );
        } else {
                const EcalPedestals::Item * aped = &peds->endcap(EEDetId(detid).hashedIndex());
                const EcalMGPAGainRatio * aGain  = &gains->endcap(EEDetId(detid).hashedIndex());          
                result.push_back( eeAlgo_.makeRecHit(*itdg, aped, aGain, weightsee, pulseee) );
        }
            
        return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalUncalibRecHitWorkerFactory, EcalUncalibRecHitWorkerMaxSample, "EcalUncalibRecHitWorkerMaxSample" );
