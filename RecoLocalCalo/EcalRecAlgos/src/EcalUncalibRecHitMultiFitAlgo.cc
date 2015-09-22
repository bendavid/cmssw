#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() : 
  _computeErrors(true),
  _doPrefit(false),
  _prefitMaxChiSq(1.0) { 
    
  _singlebx.resize(1);
  _singlebx << 0;
  
  _pulsefuncSingle.disableErrorCalculation();
  _pulsefuncSingle.setMaxIters(1);
  _pulsefuncSingle.setMaxIterWarnings(false);
    
}

/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrix &noisecor, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX) {

  uint32_t flags = 0;
  
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;
  
  double maxamplitude = -std::numeric_limits<double>::max();
  
  double pedval = 0.;
  double pedrms = 0.;
  
  SampleVector amplitudes;
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
  
  double amplitude, amperr, chisq;
  bool status = false;
  
  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (_doPrefit) {
    status = _pulsefuncSingle.DoFit(amplitudes,noisecor,pedrms,_singlebx,fullpulse,fullpulsecov);
    amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
    amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
    chisq = _pulsefuncSingle.ChiSq();
    
    if (chisq < _prefitMaxChiSq) {
      usePrefit = true;
    }
  }
  
  if (!usePrefit) {
    
//     printf("detid = %i\n",dataFrame.id().rawId());
  
    if(!_computeErrors) _pulsefunc.disableErrorCalculation();
    status = _pulsefunc.DoFit(amplitudes,noisecor,pedrms,activeBX,fullpulse,fullpulsecov);
    chisq = _pulsefunc.ChiSq();
    
    if (!status) {
      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      if (_pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }
    
    amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
    amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;
    
    
    DetId id = dataFrame.id();
    bool isbarrel = id.subdetId() == EcalBarrel;    
    
    
//     if (false && isbarrel && amplitude > 40.) {
//     if (id.rawId()==838865670) { 
    if (id.rawId() == 838882443) {
    
      gROOT->SetBatch();
      gStyle->SetOptStat(0);

      
  //     std::cout << _pulsefunc.pulsemat() << std::endl;
  //     std::cout << fullpulse << std::endl;
//       std::cout << fullpulsecov << std::endl;
      
      TH1D *hsample = new TH1D("hsample","",10,-0.5,9.5);
      for (unsigned int isample=0; isample<nsample; ++isample) {
        int ibin = hsample->FindBin(double(isample));
        hsample->SetBinContent(ibin,amplitudes[isample]);
        hsample->SetBinError(ibin,pedrms);
      }
      
      hsample->GetXaxis()->SetTitle("sample");
      hsample->GetYaxis()->SetTitle("amplitude (ADC counts)");
      hsample->SetLineColor(kBlack);
      hsample->SetMarkerColor(kBlack);
      
      TCanvas *cpulse = new TCanvas;
      hsample->Draw("E");
      
      
      const unsigned int npulse = _pulsefunc.BXs().rows();
      std::vector<TH1D*> hpulses(npulse,0);
      std::vector<double> sumerrsq(nsample,0.);
//       TH1D *hpulseit = 0;
      TH1D *hpulsesum = new TH1D("hpulsesum","",10,-0.5,9.5);
      for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
        int bx = _pulsefunc.BXs().coeff(ipulse);
        int idx = bx + 5;
        int offset = 7-3-bx;
        
        double amp = status ? _pulsefunc.X().coeff(ipulse) : 0.;
              
        TH1D *hpulse = hpulses[ipulse];
        
/*        if (bx==0) {
          hpulseit = hpulse;
        }   */   
        
        hpulse = new TH1D(TString::Format("hpulse_%i",bx),"",10,-0.5,9.5);
        for (unsigned int isample=0; isample<nsample; ++isample) {
          int ibin = hpulse->FindBin(double(isample));
          hpulsesum->Fill(double(isample),amp*fullpulse(offset+isample));
          sumerrsq[isample] += amp*amp*fullpulsecov(offset+isample,offset+isample);
          hpulse->SetBinContent(ibin,amp*fullpulse(offset+isample));          
          hpulse->SetBinError(ibin,amp*sqrt(fullpulsecov(offset+isample,offset+isample)));
        }
        
        if (bx==0) {
          hpulse->SetLineColor(kRed);
        }      
        else {
          hpulse->SetLineColor(idx+2);
        }
        
        hpulse->Draw("HISTSAME");
        

      }
      
      
      for (unsigned int isample=0; isample<nsample; ++isample) {
        int ibin = hpulsesum->FindBin(double(isample));
//         printf("isample = %i, error = %5e\n",sqrt(sumerrsq[isample]));
        hpulsesum->SetBinError(ibin,sqrt(sumerrsq[isample]));
      }      
      hpulsesum->SetLineColor(kBlue);
      hpulsesum->Draw("HISTSAME");
      
      cpulse->SaveAs(TString::Format("hpulse_%i_%i_%i_%5f.pdf",isbarrel,npulse,id.rawId(),amplitude));
      delete cpulse;
      
      TCanvas *cpulseit = new TCanvas;
      hsample->Draw("E");
//       hpulsesum->SetFillColor(kBlue);
//       hpulsesum->SetFillStyle(3344);
      hpulsesum->Draw("E2HISTSAME");
      cpulseit->SaveAs(TString::Format("hpulsesum_%i_%i_%i_%5f.pdf",isbarrel,npulse,id.rawId(),amplitude));
      delete cpulseit;
      
      delete hsample;
      for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
        delete hpulses[ipulse];
      }
    
    }
    
    
  
  }
  
  double jitter = 0.;
  
  //printf("status = %i\n",int(status));
  //printf("amplitude = %5f +- %5f, chisq = %5f\n",amplitude,amperr,chisq);
  
  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);
  
  if (!usePrefit) {
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
    }
  }
  
  return rh;
}

