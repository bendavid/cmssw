#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "TFile.h"

EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() :
  npulseEB(0),
  npulseEE(0),
  sumpulseEB(10),
  sumpulseEE(10),
  sumx0x1EB(10),
  sumx0x1EE(10),
  sumx0EB(10),
  sumx0EE(10),
  noisecovwsumEB(10),
  noisecovwsumEE(10) {
    
  }

EcalUncalibRecHitMultiFitAlgo::~EcalUncalibRecHitMultiFitAlgo() {
  
  TFile *fpulse = new TFile("fpulse.root","RECREATE");
  
  const unsigned int nsample = 10;
  
  {
    bool barrel = true;

    const int &npulse = barrel ? npulseEB : npulseEE;
    const TVectorD &sumpulse = barrel ? sumpulseEB : sumpulseEE;
    const TVectorD &sumx0 = barrel ? sumx0EB : sumx0EE;
    const TMatrixDSym &sumx0x1 = barrel ? sumx0x1EB : sumx0x1EE;
    TMatrixDSym &noisecovwsum = barrel ? noisecovwsumEB : noisecovwsumEE;
    
    sumpulse.Write("pulseEB");
    
    TVectorD npulsev(1);
    npulsev[0] = double(npulse);
    npulsev.Write("npulseEB");
    sumx0.Write("sumx0EB");
    sumx0x1.Write("sumx0x1EB");
    noisecovwsum.Write("noisecovwsumEB");    
  }

  {
    bool barrel = false;

    const int &npulse = barrel ? npulseEB : npulseEE;
    const TVectorD &sumpulse = barrel ? sumpulseEB : sumpulseEE;
    const TVectorD &sumx0 = barrel ? sumx0EB : sumx0EE;
    const TMatrixDSym &sumx0x1 = barrel ? sumx0x1EB : sumx0x1EE;
    TMatrixDSym &noisecovwsum = barrel ? noisecovwsumEB : noisecovwsumEE;
    
    sumpulse.Write("pulseEE");
    
    TVectorD npulsev(1);
    npulsev[0] = double(npulse);
    npulsev.Write("npulseEE");
    sumx0.Write("sumx0EE");
    sumx0x1.Write("sumx0x1EE");
    noisecovwsum.Write("noisecovwsumEE");    
  }  
  
  
  if (0)
  {
    bool barrel = true;

    const int &npulse = barrel ? npulseEB : npulseEE;
    const TVectorD &sumpulse = barrel ? sumpulseEB : sumpulseEE;
    const TVectorD &sumx0 = barrel ? sumx0EB : sumx0EE;
    const TMatrixDSym &sumx0x1 = barrel ? sumx0x1EB : sumx0x1EE;
    TMatrixDSym &noisecovwsum = barrel ? noisecovwsumEB : noisecovwsumEE;
    
    TMatrixDSym sumx0sumx1(nsample);
    
    
    for (unsigned int isample=0; isample<nsample; ++isample) {
      for (unsigned int jsample=0; jsample<nsample; ++jsample) {
        sumx0sumx1(isample,jsample) = sumx0(isample)*sumx0(jsample);
      }
    }
      
    double scale = 1./double(npulse);
    TMatrixDSym pulsecov = scale*sumx0x1 - scale*scale*sumx0sumx1;
    pulsecov.Write("pulsecovEB");
    
    noisecovwsum *= scale*scale;
    noisecovwsum.Write("noisecovwsumEB");
    
    TMatrixDSym pulsecovcor = pulsecov - noisecovwsum;
    pulsecovcor.Write("pulsecovcorEB");
    
    sumpulse.Write("pulseEB");
    
    printf("npulseEB = %i\n",npulse);

  } 
  
  if (0)
  {
    bool barrel = false;

    const int &npulse = barrel ? npulseEB : npulseEE;
    const TVectorD &sumpulse = barrel ? sumpulseEB : sumpulseEE;
    const TVectorD &sumx0 = barrel ? sumx0EB : sumx0EE;
    const TMatrixDSym &sumx0x1 = barrel ? sumx0x1EB : sumx0x1EE;
    TMatrixDSym &noisecovwsum = barrel ? noisecovwsumEB : noisecovwsumEE;
    
    TMatrixDSym sumx0sumx1(nsample);
    
    
    for (unsigned int isample=0; isample<nsample; ++isample) {
      for (unsigned int jsample=0; jsample<nsample; ++jsample) {
        sumx0sumx1(isample,jsample) = sumx0(isample)*sumx0(jsample);
      }
    }
      
    double scale = 1./double(npulse);
    TMatrixDSym pulsecov = scale*sumx0x1 - scale*scale*sumx0sumx1;
    pulsecov.Write("pulsecovEE");
    
    noisecovwsum *= scale*scale;
    noisecovwsum.Write("noisecovwsumEE");
    
    TMatrixDSym pulsecovcor = pulsecov - noisecovwsum;
    pulsecovcor.Write("pulsecovcorEE");
    
    sumpulse.Write("pulseEE");
    
    printf("npulseEE = %i\n",npulse);

  }  
  
  
  fpulse->Close();
  
  
}

  
/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const TMatrixDSym &noisecor, const TVectorD &fullpulse, const TMatrixDSym &fullpulsecov, std::set<int> activeBX) {

  uint32_t flags = 0;
  
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
  double maxamplitude = -std::numeric_limits<double>::max();
  
  double pedval = 0.;
  double pedrms = 0.;
  
  std::vector<double> amplitudes(nsample);
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {
    
    const EcalMGPASample &sample = dataFrame.sample(iSample);
    
    double amplitude = 0.;
    int gainId = sample.gainId();
    
    double pedestal = 0.;
    double pederr = 0.;
    double gainratio = 1.;
        
    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      pederr = aped->rms_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      pederr = aped->rms_x12;
      gainratio = 1.;
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      pederr = aped->rms_x6;
      gainratio = aGain->gain12Over6();
    }

    amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    
    if (gainId == 0) {
      //saturation
      amplitude = (4095. - pedestal) * gainratio;
    }
        
    amplitudes[iSample] = amplitude;
    
    if (amplitude>maxamplitude) {
    //if (iSample==5) {
      maxamplitude = amplitude;
      pedval = pedestal;
      pedrms = pederr*gainratio;
    }    
        
  }
    
  if (1) {
    if (maxamplitude/pedrms > 1000.) {
    //if (maxamplitude/pedrms > 3.) {
    //if (true) {
      DetId detid(dataFrame.id());
      bool barrel = detid.subdetId()==EcalBarrel;
    
      int &npulse = barrel ? npulseEB : npulseEE;
      TVectorD &sumpulse = barrel ? sumpulseEB : sumpulseEE;
      TVectorD &sumx0 = barrel ? sumx0EB : sumx0EE;
      TMatrixDSym &sumx0x1 = barrel ? sumx0x1EB : sumx0x1EE;
      TMatrixDSym &noisecovwsum = barrel ? noisecovwsumEB : noisecovwsumEE;
      
      double scale = 1./maxamplitude;
      for (unsigned int isample=0; isample<nsample; ++isample) {
        printf("isample = %i, amplitude = %5f\n",isample,amplitudes[isample]);
        sumpulse[isample] += amplitudes[isample];
        sumx0[isample] += scale*amplitudes[isample];
        noisecovwsum += scale*scale*pedrms*pedrms*noisecor;
        for (unsigned int jsample=0; jsample<nsample; ++jsample) {
          sumx0x1(isample,jsample) += scale*scale*amplitudes[isample]*amplitudes[jsample];
        }
      }
      ++npulse;
    }
    return EcalUncalibratedRecHit();
  }
  

  bool status = _pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
  double chisq = _pulsefunc.ChiSq();
  
  if (!status) {
    edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
  }

  unsigned int ipulseintime = std::distance(activeBX.begin(),activeBX.find(0));
  double amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
  double amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
  
  double jitter = 0.;
  
  //printf("amplitude = %5f +- %5f, chisq = %5f\n",amplitude,amperr,chisq);
  
  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);
  for (std::set<int>::const_iterator bxit = activeBX.begin(); bxit!=activeBX.end(); ++bxit) {
    int ipulse = std::distance(activeBX.begin(),bxit);
    if(*bxit==0) {
      rh.setOutOfTimeAmplitude(*bxit + 5,0.);
    } else {
      rh.setOutOfTimeAmplitude(*bxit + 5, status ? _pulsefunc.X()[ipulse] : 0.);
    }
  }

  return rh;
}

