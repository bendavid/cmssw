#ifndef RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMaxSampleAlgo_HH
#define RecoLocalCalo_EcalRecAlgos_EcalUncalibRecHitMaxSampleAlgo_HH

/** \class EcalUncalibRecHitMaxSampleAlgo
  *  Amplitude reconstucted by the difference MAX_adc - min_adc
  *  jitter is sample number of MAX_adc, pedestal is min_adc
  *
  *  \author G. Franzoni, E. Di Marco
  */

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitRecAbsAlgo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

template<class C> class EcalUncalibRecHitMaxSampleAlgo : public EcalUncalibRecHitRecAbsAlgo<C>
{
  
 public:
  
  virtual ~EcalUncalibRecHitMaxSampleAlgo<C>() { };
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const std::vector<double> &weights, const std::vector<double> &pulse);
  virtual EcalUncalibratedRecHit makeRecHit(const C& dataFrame, const double* pedestals,
                                            const double* gainRatios,
                                            const EcalWeightSet::EcalWeightMatrix** weights, 
                                            const EcalWeightSet::EcalChi2WeightMatrix** chi2Matrix) { return EcalUncalibratedRecHit(); }

 private:

};

/// compute rechits
template<class C> EcalUncalibratedRecHit  
EcalUncalibRecHitMaxSampleAlgo<C>::makeRecHit(const C& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const std::vector<double> &weights, const std::vector<double> &pulse) {

  double maxamplitude = -std::numeric_limits<double>::max();
  double maxpedestal  = 4095;
  double maxjitter    = -1;
  double maxchi2      = -1;
  double maxpederr = 0.;
  //bool isSaturated = 0;
  uint32_t maxflags = 0;
  int maxgain = 1.;
  std::vector<float> amplitudes(C::MAXSAMPLES);
  for(int16_t iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
    
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double pederr = 0.;
    double gainratio = 1.;
    int gain = 1;
    
    uint32_t flags = 0;
    
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      pederr = aped->rms_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gain = 1;
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      pederr = aped->rms_x12;
      gainratio = 1.;
      gain = 12;
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      pederr = aped->rms_x6;
      gainratio = aGain->gain12Over6();
      gain = 6;
    }

    amplitude = double(((double)(sample.adc()) - pedestal) * gainratio);
    
    
    if (gainId == 0) {
      flags = EcalUncalibratedRecHit::kSaturated;
      amplitude = double((4095. - pedestal) * gainratio);
    }
    
    //printf("isamp = %i, adc = %i, pedestal = %5f, gainratio = %5f, amplitude = %5f\n",int(iSample),sample.adc(),pedestal,gainratio,amplitude);
    
    amplitudes[iSample] = amplitude;
    
    if (amplitude>maxamplitude) {
    //if (iSample==5) {
      maxamplitude = amplitude;
      maxpedestal = pedestal;
      maxjitter = (iSample-5);
      maxflags = flags;
      maxpederr = pederr*gainratio;
      maxgain = gain;
    }
    
    

  }// loop on samples
  
  //loop on samples again to compute multisample pulse for each bx
//   std::vector<float> pulseAmplitudes(C::MAXSAMPLES);
//   maxamplitude = -std::numeric_limits<double>::max();
//   for(int16_t iSample = 0; iSample < C::MAXSAMPLES; iSample++) {
//     //loop over available samples
//     double amplitude = 0.;
//     double sumweight = 0.;
//     for(int16_t jSample = std::max(0,iSample-5); jSample < std::min(C::MAXSAMPLES,iSample+5); jSample++) {
//       double weight = weights[jSample-iSample+5];
//       double pulseval = pulse[jSample-iSample+5];
//       amplitude += weight*amplitudes[jSample];
//       sumweight += weight*pulseval;
//     }
//     amplitude/=sumweight;
//     pulseAmplitudes[iSample] = amplitude;
//     if (amplitude>maxamplitude) {
//       maxamplitude = amplitude;
//       maxjitter = (iSample-5);
//     }
//   }

  

  
//   if (maxjitter!=0) {
//     printf("out of time hit with jitter = %i, amplitude = %5f\n",int(maxjitter),maxamplitude);
//   }
      
  EcalUncalibratedRecHit rh( dataFrame.id(), maxamplitude , maxpedestal, maxjitter, maxchi2, maxflags );
  //rh.setAmplitudes(pulseAmplitudes);
  rh.setAmplitudes(amplitudes);
  rh.setPedrms(maxpederr);
  rh.setGain(maxgain);
  return rh;
}

#endif
