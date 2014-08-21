#include "RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerMulti.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "CondFormats/DataRecord/interface/EcalTimeBiasCorrectionsRcd.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompBK.h"
#include "TDecompChol.h"
#include "vdt/vdtMath.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TH1D.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

EcalRecHitWorkerMulti::EcalRecHitWorkerMulti(const edm::ParameterSet&ps, edm::ConsumesCollector& c) :
  EcalRecHitWorkerBaseClass(ps,c), invsamplecorEBg12(10), invsamplecorEEg12(10),
  invsamplecorEBg6(10), invsamplecorEEg6(10),
  invsamplecorEBg1(10), invsamplecorEEg1(10), sampleunitcov(10)
{
    rechitMaker_ = new EcalRecHitSimpleAlgo();
        v_chstatus_ = ps.getParameter<std::vector<int> >("ChannelStatusToBeExcluded");
        v_DB_reco_flags_ = ps.getParameter<std::vector<int> >("flagsMapDBReco");
        killDeadChannels_ = ps.getParameter<bool>("killDeadChannels");
        laserCorrection_ = ps.getParameter<bool>("laserCorrection");
	EBLaserMIN_ = ps.getParameter<double>("EBLaserMIN");
	EELaserMIN_ = ps.getParameter<double>("EELaserMIN");
	EBLaserMAX_ = ps.getParameter<double>("EBLaserMAX");
	EELaserMAX_ = ps.getParameter<double>("EELaserMAX");

	
// 	// Traslate string representation of flagsMapDBReco into enum values 
// 	const edm::ParameterSet & p=ps.getParameter< edm::ParameterSet >("flagsMapDBReco");
// 	std::vector<std::string> recoflagbitsStrings = p.getParameterNames();
// 	v_DB_reco_flags_.resize(32); 
// 
// 	for (unsigned int i=0;i!=recoflagbitsStrings.size();++i){
// 	  EcalRecHit::Flags recoflagbit = (EcalRecHit::Flags)
// 	    StringToEnumValue<EcalRecHit::Flags>(recoflagbitsStrings[i]);
// 	  std::vector<std::string> dbstatus_s =  
// 	    p.getParameter<std::vector<std::string> >(recoflagbitsStrings[i]);
// 	  std::vector<uint32_t> dbstatuses;
// 	  for (unsigned int j=0; j!= dbstatus_s.size(); ++j){
// 	    EcalChannelStatusCode::Code  dbstatus  = (EcalChannelStatusCode::Code)
// 	      StringToEnumValue<EcalChannelStatusCode::Code>(dbstatus_s[j]);
// 	    dbstatuses.push_back(dbstatus);
// 	  }
// 
// 	  v_DB_reco_flags_[recoflagbit]=dbstatuses;
// 	}  
// 
//     flagmask_=0;
//     flagmask_|=    0x1<<EcalRecHit::kNeighboursRecovered;
//     flagmask_|=    0x1<<EcalRecHit::kTowerRecovered;
//     flagmask_|=    0x1<<EcalRecHit::kDead;
//     flagmask_|=    0x1<<EcalRecHit::kKilled;
//     flagmask_|=    0x1<<EcalRecHit::kTPSaturated;
//     flagmask_|=    0x1<<EcalRecHit::kL1SpikeFlag; 
    
    const int nsample = 10;
    
    double samplecorvEBg12[nsample] = {1.00000, 0.71073, 0.55721, 0.46089, 0.40449, 0.35931, 0.33924, 0.32439, 0.31581, 0.30481};
    double samplecorvEEg12[nsample] = {1.00000, 0.71373, 0.44825, 0.30152, 0.21609, 0.14786, 0.11772, 0.10165, 0.09465, 0.08098};
    
    double samplecorvEBg6[nsample] = {1.00000, 0.70946, 0.58021, 0.49846, 0.45006, 0.41366, 0.39699, 0.38478, 0.37847, 0.37055};
    double samplecorvEEg6[nsample] = {1.00000, 0.71217, 0.47464, 0.34056, 0.26282, 0.20287, 0.17734, 0.16256, 0.15618, 0.14443};
    
    double samplecorvEBg1[nsample] = {1.00000, 0.73354, 0.64442, 0.58851, 0.55425, 0.53082, 0.51916, 0.51097, 0.50732, 0.50409};
    double samplecorvEEg1[nsample] = {1.00000, 0.72698, 0.62048, 0.55691, 0.51848, 0.49147, 0.47813, 0.47007, 0.46621, 0.46265};
    
    //fill correlation matrices
    for (int i=0; i<nsample; ++i) {
      for (int j=0; j<nsample; ++j) {
        int vidx = std::abs(j-i);
        invsamplecorEBg12(i,j) = samplecorvEBg12[vidx];
        invsamplecorEEg12(i,j) = samplecorvEEg12[vidx];
        invsamplecorEBg6(i,j) = samplecorvEBg6[vidx];
        invsamplecorEEg6(i,j) = samplecorvEEg6[vidx];
        invsamplecorEBg1(i,j) = samplecorvEBg1[vidx];
        invsamplecorEEg1(i,j) = samplecorvEEg1[vidx];        
      }
    }
//     invsamplecorEBg12.Invert();
//     invsamplecorEEg12.Invert();
//     invsamplecorEBg6.Invert();
//     invsamplecorEEg6.Invert();    
//     invsamplecorEBg1.Invert();
//     invsamplecorEEg1.Invert();    
    
    for (int i=0; i<nsample; ++i) {
      sampleunitcov(i,i) = 1.;
    }
    
    
    hpulseEB = new TH1D("hpulseEB","",10,-0.5,9.5);
    hpulseEE = new TH1D("hpulseEE","",10,-0.5,9.5);
    
    hpulseprofEB = new TProfile("hpulseprofEB","",10,-0.5,9.5,"s");
    hpulseprofEE = new TProfile("hpulseprofEE","",10,-0.5,9.5,"s");
    
}


void EcalRecHitWorkerMulti::set(const edm::EventSetup& es)
{
        es.get<EcalIntercalibConstantsRcd>().get(ical);
        es.get<EcalTimeCalibConstantsRcd>().get(itime);
        es.get<EcalTimeOffsetConstantRcd>().get(offtime);
        // for the time correction methods
        es.get<EcalTimeBiasCorrectionsRcd>().get(timeCorrBias);        
        es.get<EcalADCToGeVConstantRcd>().get(agc);
        es.get<EcalChannelStatusRcd>().get(chStatus);
        if ( laserCorrection_ ) es.get<EcalLaserDbRecord>().get(laser);
        es.get<CaloTopologyRecord>().get(theCaloTopology);
        es.get<CaloGeometryRecord>().get(caloGeometry);   
}


void
EcalRecHitWorkerMulti::run( const edm::Event & evt,
                const EcalUncalibratedRecHitCollection& uncalibRHs,
                EcalRecHitCollection & result)
{
    std::vector<const EcalUncalibratedRecHit*> outuncalib;
    std::map<DetId, std::vector<double> > calibpulses;
    std::map<DetId, double> pedrmss;
  
    std::vector<std::string> pulsenames;
    pulsenames.push_back("pulse0");
    pulsenames.push_back("pulse1");
    pulsenames.push_back("pulse2");
    pulsenames.push_back("pulse3");
    pulsenames.push_back("pulse4");
    pulsenames.push_back("pulse5");
    pulsenames.push_back("pulse6");
    pulsenames.push_back("pulse7");
    pulsenames.push_back("pulse8");
    pulsenames.push_back("pulse9");
    pulsenames.push_back("pedestal");
    
    DetId globalmax;
    double globalmaxval = 0.;
    
    for (const EcalUncalibratedRecHit &uncalibRH : uncalibRHs) {

      
      DetId detid=uncalibRH.id();
      bool barrel = detid.subdetId()==EcalBarrel;

      
        EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
        EcalChannelStatusCode chStatusCode = 1;
        if ( chit != chStatus->end() ) {
                chStatusCode = *chit;
        } else {
                edm::LogError("EcalRecHitError") << "No channel status found for xtal " 
                        << detid.rawId() 
                        << "! something wrong with EcalChannelStatus in your DB? ";
        }
        if ( v_chstatus_.size() > 0) {
                uint16_t code = chStatusCode.getStatusCode() & 0x001F;
                std::vector<int>::const_iterator res = std::find( v_chstatus_.begin(), v_chstatus_.end(), code );
                if ( res != v_chstatus_.end() ) {
                        continue;
                }
        }

        // find the proper flag for the recHit
        // from a configurable vector
        // (see cfg file for the association)
        uint32_t recoFlag = 0;
        uint16_t statusCode = chStatusCode.getStatusCode() & 0x001F;
        if ( statusCode < v_DB_reco_flags_.size() ) {
                // not very nice...
                recoFlag = v_DB_reco_flags_[ statusCode ];  
        } else {
                edm::LogError("EcalRecHitError") << "Flag " << statusCode 
                        << " in DB exceed the allowed range of " << v_DB_reco_flags_.size();
        }

          //float offsetTime = 0; // the global time phase
          const EcalIntercalibConstantMap& icalMap = ical->getMap();  
//       if ( detid.subdetId() == EcalEndcap ) {
//           rechitMaker_->setADCToGeVConstant( float(agc->getEEValue()) );
//                   //offsetTime = offtime->getEEValue();
//       } else {
//           rechitMaker_->setADCToGeVConstant( float(agc->getEBValue()) );
//                   //offsetTime = offtime->getEBValue();
//       }

      double agcval = barrel ? agc->getEBValue() : agc->getEEValue();
          
      // first intercalibration constants
      EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detid);
      EcalIntercalibConstant icalconst = 1;
      if( icalit!=icalMap.end() ) {
          icalconst = (*icalit);
      } else {
          edm::LogError("EcalRecHitError") << "No intercalib const found for xtal "
                                          << detid.rawId()
                                          << "! something wrong with EcalIntercalibConstants in your DB? ";
      }

      // get laser coefficient
      float lasercalib = 1.;
      if ( laserCorrection_ )	lasercalib = laser->getLaserCorrection( detid, evt.time());
          

//       double offsetTime = barrel ? offtime->getEBValue() : offtime->getEEValue();
// 
//       
//       // get time calibration coefficient
//       const EcalTimeCalibConstantMap & itimeMap = itime->getMap();  
//       EcalTimeCalibConstantMap::const_iterator itime = itimeMap.find(detid);
//       EcalTimeCalibConstant itimeconst = 0;
//       if( itime!=itimeMap.end() ) {
//           itimeconst = (*itime);
//       } else {
//           edm::LogError("EcalRecHitError") << "No time calib const found for xtal "
//                           << detid.rawId()
//                           << "! something wrong with EcalTimeCalibConstants in your DB? ";
//           }
//            
//       printf("barrel = %i, offsetTime = %5f, timeconst = %5f\n",int(barrel),offsetTime,float(itimeconst));
          
      // make the rechit and put in the output collection, unless recovery has to take care of it
      if (recoFlag<=EcalRecHit::kLeadingEdgeRecovered || !killDeadChannels_) {
        outuncalib.push_back(&uncalibRH);
        std::vector<double> &calibpulse =  calibpulses[detid];
        
        for (float amplitude : uncalibRH.amplitudes()) {
          calibpulse.push_back(amplitude*icalconst*lasercalib*agcval);
        }
        
        double pedrms = uncalibRH.pedrms()*icalconst*lasercalib*agcval;
        pedrmss[detid] = pedrms;
        
        if (calibpulse[5]>globalmaxval) {
          globalmax = detid;
          globalmaxval = calibpulse[5];
        }
        

      }
    }
    
    for (const EcalUncalibratedRecHit *uncalibRHp : outuncalib) {
      const EcalUncalibratedRecHit &uncalibRH = *uncalibRHp;
      DetId detid=uncalibRH.id();
      
      bool barrel = detid.subdetId()==EcalBarrel;
     
      const std::vector<double> &calibpulse = calibpulses[detid];
      double pedrms = pedrmss[detid];
      
      
      bool fillpulse = detid==globalmax;
      
      double sumpulse = 0.;
      for (double pulseval : calibpulse) {
        sumpulse += pulseval;
      }
      //fillpulse = false;
      fillpulse = barrel ? calibpulse[5]>50. : calibpulse[5]>100.;
      //fillpulse = calibpulse[5]>5. && calibpulse[5]<15.;
      
      fillpulse = false;
      
      if (fillpulse) {
        //double pulseweight = sumpulse/pedrms/pedrms;
        
        const double alpha = barrel ? 1.3789 : 1.6836;
        const double beta = barrel ? 1.5006 : 1.4173; 
        
        std::vector<double> alphav(1,alpha);
        std::vector<double> betav(1,beta);
        std::vector<double> tmaxv(1,5.);
        
        TMatrixDSym inparamcov(3);
        
        double pulseweight = 1.;
        TH1D *hpulse = barrel ? hpulseEB : hpulseEE;
        TProfile *hpulseprof = barrel ? hpulseprofEB : hpulseprofEE;
        PulseChiSqGlobal &pulseglobal = barrel ? pulseglobalEB : pulseglobalEE;
        for (int isample = 0; isample<10; ++isample) {
          hpulse->Fill(isample,pulseweight*calibpulse[isample]);
          hpulseprof->Fill(isample,calibpulse[isample]/sumpulse);
        }
        pulseglobal.samples().emplace_back(calibpulse.size(),calibpulse.data());
        pulseglobal.invcovs().emplace_back((1./(pedrms*pedrms))*invsamplecor(barrel,uncalibRH.gain()));
        
        const TMatrixDSym &theinvsamplecor = invsamplecor(barrel,uncalibRH.gain());
        TMatrixDSym invsamplecov = (1./(pedrms*pedrms))*theinvsamplecor;        
        
        PulseChiSqTime pulsefunc(calibpulse,invsamplecov,alphav,betav,tmaxv,inparamcov);

        ROOT::Minuit2::Minuit2Minimizer minim;
        minim.SetPrintLevel(2);
        minim.SetStrategy(2);
        //if (doprint) minim.SetPrintLevel(9);
        minim.SetFunction(pulsefunc);
        
        minim.SetLowerLimitedVariable(0,"amp",calibpulse[5],0.001,0.);
        minim.SetLimitedVariable(1,"alpha",alpha,0.001,1.,3.);
        minim.SetLimitedVariable(2,"beta",beta,0.001,1.,3.);
        minim.SetLimitedVariable(3,"tmax",5.,0.001,4.8,5.2);
        
        minim.Minimize();
        
        
      }      
      
    }
    
    //return;
//     const double eseedEB = 0.23;
//     //const double etseedEB = 0.;
//     
//     const double eseedEE = 0.6;
//     //const double etseedEE = 0.15;
//     
//     const double eseedlowEB = 0.23;
//     const double eseedlowEE = 0.6;
    
    std::vector<std::pair<DetId, unsigned int> > localmaxima;
    std::map<DetId,std::set<int> > activeBXs;
    
    std::set<int> activeAll;
    activeAll.insert(-5);
    activeAll.insert(-4);
    activeAll.insert(-3);
    activeAll.insert(-2);
    activeAll.insert(-1);
    activeAll.insert(0);
    activeAll.insert(1);
    activeAll.insert(2);
    activeAll.insert(3);
    activeAll.insert(4);
    
    //activeAll.insert(0);
    
    //make rechits
    for (const EcalUncalibratedRecHit *uncalibRHp : outuncalib) {
      //continue;
      const EcalUncalibratedRecHit &uncalibRH = *uncalibRHp;
      DetId detid=uncalibRH.id();
      
      bool barrel = detid.subdetId()==EcalBarrel;
      
      const std::vector<double> &calibpulse = calibpulses[detid];
      //std::set<int> &activeBX = activeBXs[detid];
      //const std::set<int> &activeBX = activeAll;
      std::set<int> activeBX;
      std::set<int> &activeBXfinal = activeBXs[detid];
      
      double pedrms = pedrmss[detid];
      
      const double alpha = barrel ? 1.36745 : 1.59918;
      const double beta = barrel ? 1.51341 : 1.46568;      
      const double tmax = barrel ? 5.04159 : 5.09323;      
  
//       const double alpha = barrel ? 1.3789 : 1.6836;
//       const double beta = barrel ? 1.5006 : 1.4173;      
//       const double tmax = barrel ? 5.04583 : 5.11208;
//      double tmaxerr = barrel ? 0.0140 : 0.0109;
//      tmaxerr = 0.05;

      const unsigned int nsample = calibpulse.size();
      
      TVectorD sampvec(nsample,calibpulse.data());
      
      double shapeerr = 1e-2;
      
      TMatrixDSym invsamplecov = pedrms*pedrms*invsamplecor(barrel,uncalibRH.gain());
      for (unsigned int isample = 0; isample<nsample; ++isample) {
        invsamplecov(isample,isample) += pow(shapeerr*calibpulse[isample],2);
      }
      std::vector<double> caliberr(nsample);
      for (unsigned int isample=0; isample<nsample; ++isample) {
        caliberr[isample] = sqrt(invsamplecov(isample,isample));
      }
      invsamplecov.Invert();      
      
      double chisq = invsamplecov.Similarity(sampvec);
      
      
      std::vector<double> calibpulsenow = calibpulse;
      unsigned int npulse = 0;
      //double chisq = 0.;
      std::vector<double> fitvals;
      std::vector<double> fiterrs;
      TMatrixD pulsemat;
      while (true) {
        double maxsample = -std::numeric_limits<double>::max();
        std::set<int> maxsamples;
        for (unsigned int isample = 0; isample<nsample; ++isample) {
          int bx = int(isample)-5;
          bool localmax = !activeBX.count(bx) && !activeBX.count(bx-1) && !activeBX.count(bx+1);
          if (localmax && calibpulsenow[isample]>maxsample) {
            maxsample = calibpulsenow[isample];
          }
        }
        if (maxsample<0.) {
          break;
        }
        for (unsigned int isample = 0; isample<nsample; ++isample) {
          int bx = int(isample)-5;
          bool localmax = !activeBX.count(bx) && !activeBX.count(bx-1) && !activeBX.count(bx+1);
          if (!localmax) continue;
          if (calibpulsenow[isample]==maxsample || (isample>0 && calibpulsenow[isample-1]==maxsample) || (isample<(nsample-1) && calibpulsenow[isample+1]==maxsample) ) {
            maxsamples.insert(isample);
          }
        }
        
        if (!activeBX.count(0)) maxsamples.insert(5);
        
        unsigned int npulsenow = activeBX.size() + 1;        
        
        double sigmax = 0.;
        int bxmax = 0;
        double chisqmax = 0.;
        std::vector<double> fitvalsmax(npulsenow);
        std::vector<double> fiterrsmax(npulsenow);
        TMatrixD pulsematmax(nsample,npulsenow);        
        double chisqintime = 0.;
        
        for (int isamplemax : maxsamples) {
          //printf("start fit attempt\n");
          std::set<int> activeBXnow(activeBX);
          int newbx = int(isamplemax)-5;
          activeBXnow.insert(newbx);
          
          TMatrixD pulsemat(nsample,npulsenow);
          for (std::set<int>::const_iterator pulsebx = activeBXnow.begin(); pulsebx!=activeBXnow.end(); ++pulsebx) {
            int ipulse = std::distance(activeBXnow.begin(),pulsebx);
            for (unsigned int isample=0; isample<nsample; ++isample) {
              double t = double(isample);
              double val = pulseShape(t,alpha,beta,tmax+double(*pulsebx));
              pulsemat[isample][ipulse] = val;
            }
          }
          
          
          PulseChiSq pulsefunc(calibpulse,pulsemat,invsamplecov);
          
          ROOT::Minuit2::Minuit2Minimizer minim;
          minim.SetStrategy(2);
          minim.SetFunction(pulsefunc);
          
          for (std::set<int>::const_iterator pulsebx = activeBXnow.begin(); pulsebx!=activeBXnow.end(); ++pulsebx) {
            int ipulse = std::distance(activeBXnow.begin(),pulsebx);
            minim.SetLowerLimitedVariable(ipulse,pulsenames[ipulse],calibpulsenow[*pulsebx+5],0.001,0.);
          }
                    
          bool status = minim.Minimize();
          if (!status) {
            continue;
          }
          
          const double chisqnow = minim.MinValue();
          if (newbx==0) {
            chisqintime = chisqnow;
          }
          double sig = sqrt(std::max(0.,chisq-chisqnow));
          //printf("chisq0 = %5f, chisqnow = %5f, sig = %5f\n",chisq0,chisqnow,sig);
          if (sig>sigmax) {
            sigmax = sig;
            bxmax = newbx;
            pulsematmax = pulsemat;
            chisqmax = chisqnow;
            for (unsigned int ipulse=0; ipulse<npulsenow; ++ipulse) {
              fitvalsmax[ipulse] = minim.X()[ipulse];
              fiterrsmax[ipulse] = minim.Errors()[ipulse];
            }
          }
        }
        
        //printf("sigmax = %5f\n",sigmax);
        //bool significant = sigmax>5 || (activeBX.size()==0 && chisq>20.);
        //bool significant = sigmax > 3. || (activeBX.size()==0 && sigmax>2.);
        
        if (!activeBX.count(0)) {
          sigmax = sqrt(std::max(0.,chisqintime-chisqmax));
        }
        
        bool significant = sigmax > 5.;
        
        
        if (!significant) break;
        
        activeBX.insert(bxmax);
        npulse = activeBX.size();
        chisq = chisqmax;
        fitvals = fitvalsmax;
        fiterrs = fiterrsmax;    
        pulsemat.ResizeTo(nsample,npulse);
        pulsemat = pulsematmax;
        for (unsigned int isample=0; isample<nsample; ++isample) {
          calibpulsenow[isample] = calibpulse[isample];
          for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
            int ipulse = std::distance(activeBX.begin(),pulsebx);
            calibpulsenow[isample] -= fitvals[ipulse]*pulsemat(isample,ipulse);
          }
        }
        
      }
      
      double maxpulseval = -std::numeric_limits<double>::max();
      for (unsigned int isample = 0; isample<nsample; ++isample) {
        if (calibpulse[isample]>maxpulseval) {
          maxpulseval = calibpulse[isample];
        }
      }
      
      bool doprint = false;
      //doprint = maxpulseval>0.6 && (activeBX.size()>1 || !activeBX.count(0));
      //doprint = maxpulseval>5. && (activeBX.size()>1);
      //bool doprint = fitvals[5]>5.;
//       for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
//         double sig = fitvals[ipulse]/fiterrs[ipulse];
//         if (sig>5.) doprint = true;
//       }

      
      //TProfile *hpulseprof = barrel ? hpulseprofEB : hpulseprofEE;
      double maxdev = 0.;
//       for (unsigned int isample=0; isample<nsample; ++isample) {
//         double sampval = calibpulse[isample]/sumpulse;
//         double pulseval = hpulseprof->GetBinContent(hpulseprof->FindFixBin(isample));
//         double pulserms = hpulseprof->GetBinError(hpulseprof->FindFixBin(isample));
//         double pedrmsnorm = pedrms/sumpulse;
//         double fullerr = sqrt(pulserms*pulserms+pedrmsnorm*pedrmsnorm);
//         double dev = std::abs(sampval-hpulseprof->GetBinContent(hpulseprof->FindFixBin(isample)))/fullerr;
//         if (doprint) printf("isample = %i, sampval = %5f, fitpulseval = %5f, pulseval = %5f, pulserms = %5f, pedrmsnorm = %5f, fullerr = %5f, dev = %5f\n",isample,sampval,fitpulsevals[isample],pulseval,pulserms,pedrmsnorm,fullerr, dev);
//         if (dev>maxdev) {
//           maxdev = dev;
//         }
// //         if (dev > 5.) {
// //           doprint = true;
// //         }
//       }
      

      
      //doprint &= fitvals[5]>5.;
      
      //bool doprint = fitvals[5]>5.;
      //bool doprint = detid==globalmax && fitvals[5]>5.;
      
      
      
      
      double ootsig = 0;
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        int ipulse = std::distance(activeBX.begin(),pulsebx);   
        double pulsesig = fitvals[ipulse]/fiterrs[ipulse];
        if ((*pulsebx)!=0 && pulsesig>ootsig) ootsig=pulsesig;
        //if ( (*pulsebx)!=0 && pulsesig>5.) doprint = true; 
      }

      //double eseed = barrel ? 0.23 : 0.6;
      
      double maxfitval = 0.;
      int maxfitbx = 0;
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        int ipulse = std::distance(activeBX.begin(),pulsebx);   
        if (fitvals[ipulse]>maxfitval) {
          maxfitval = fitvals[ipulse];
          maxfitbx = *pulsebx;
        }
      }
      
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        if ((*pulsebx)==0) continue;
        int ipulse = std::distance(activeBX.begin(),pulsebx);   
        double pulsesig = fitvals[ipulse]/fiterrs[ipulse];
        //bool maxpulse = ( (ipulse==0 || fitvals[ipulse]>fitvals[ipulse-1]) && ( ipulse==(int(npulse)-1) || fitvals[ipulse]>fitvals[ipulse+1]) );
        if (pulsesig>10. || ipulse==0 || (*pulsebx)==maxfitbx) {
          //doprint = fitvals[5]>5.;
          //doprint = true;
          activeBXfinal.insert(*pulsebx);
          
          if (pulsesig>5.) {
            const int xtalradius = 5;
            //3x3 neighbours for local maxima search
            CaloNavigator<DetId> cursor = CaloNavigator<DetId>( detid, barrel ? theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel) : theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap) );
            for (int dx = -int(xtalradius)/2; dx<=int(xtalradius)/2; ++dx) {
              for (int dy = -int(xtalradius)/2; dy<=int(xtalradius)/2; ++dy) {
                  
                if (dx==0 && dy==0) continue;
                
                cursor.home();
                cursor.offsetBy( dx, dy );
                const DetId &neighbour = *cursor;

                const auto &pulsepair = calibpulses.find(neighbour);
                if (pulsepair==calibpulses.end()) continue;            
                              
                activeBXs[neighbour].insert(*pulsebx);
              }
            }
          }
          
        }
      }      
      
      if (doprint) {
        printf("prefit barrel = %i, pederr = %5f, maxdev = %5f, ", int(barrel), uncalibRH.pedrms(), maxdev);
        if (barrel) {
          EBDetId ebid(detid);
          printf("ieta = %i, iphi=%i\n",ebid.ieta(),ebid.iphi());
        }
        else {
          EEDetId eeid(detid);
          printf("ix = %i, iy=%i\n",eeid.ix(),eeid.iy());          
        }
        for (unsigned int isample = 0; isample<calibpulse.size(); ++isample) {
          printf(" %5f",calibpulse[isample]);
        }
        printf("\n");    
        printf("active BXs:");
        for (int bx : activeBX) {
          printf(" %i",bx);
        }
        printf("\n");
        
        printf("chisq = %5f\n", chisq);
        for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
          int ipulse = std::distance(activeBX.begin(),pulsebx);
          printf("ipulse = %i, bx = %i, energy = %5f +- %5f, sig = %5f\n",ipulse,*pulsebx,fitvals[ipulse],fiterrs[ipulse],fitvals[ipulse]/fiterrs[ipulse]);
          //printf("ipulse = %i, bx = %i, energy = %5f +- %5f, sig = %5f, tmax = %5f +- %5f\n",ipulse,*pulsebx,fitvals[ipulse],fiterrs[ipulse],fitvals[ipulse]/fiterrs[ipulse],minim.X()[1+2*ipulse],minim.Errors()[1+2*ipulse]);
        }   
        
/*        for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
          int ipulse = std::distance(activeBX.begin(),pulsebx);
          if (fitvals[ipulse]/fiterrs[ipulse]<0.1) continue;
          for (std::set<int>::const_iterator jpulsebx = activeBX.begin(); jpulsebx!=activeBX.end(); ++jpulsebx) {
            int jpulse = std::distance(activeBX.begin(),jpulsebx);
            if (jpulse==ipulse) {
              minim.SetFixedVariable(jpulse,pulsenames[jpulse],0.);
            }
            else {
              minim.SetLowerLimitedVariable(jpulse,pulsenames[jpulse],fitvals[jpulse],0.001,0.);
            }
          }
          minim.Minimize();
          double chisqalt = minim.MinValue();
          double sigalt = sqrt(chisqalt-chisq);

          printf("ipulse = %i, sigalt = %5f\n",ipulse,sigalt);
        }*/        
        
        bool dodraw = true;
        dodraw = false;
        if (dodraw) {
          TCanvas *c = new TCanvas;
          TH1D *hsamples = new TH1D("hsamples","",10,-0.5,9.5);
          TH1D *hfullpulse = new TH1D("hfullpulse","",10,-0.5,9.5);
          std::vector<TH1D*> hists;
          for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
            TH1D *hpulse = new TH1D(TString::Format("hpulse_%i",ipulse),"",10,-0.5,9.5);
            hpulse->SetLineColor(ipulse+1);
            for (unsigned int isample=0; isample<nsample; ++isample) {
              //double pulseval = pulseShape(double(isample),alpha, beta,tmaxv[ipulse]);
              double pulseval = pulsemat(isample,ipulse);
//               hpulse->Fill(isample,fitvals[ipulse]*pulsemat(isample,ipulse));
//               hfullpulse->Fill(isample,fitvals[ipulse]*pulsemat(isample,ipulse));
              hpulse->Fill(isample,fitvals[ipulse]*pulseval);
              hfullpulse->Fill(isample,fitvals[ipulse]*pulseval);              
            }
            hists.push_back(hpulse);
          }
          for (unsigned int isample=0; isample<nsample; ++isample) {
            hsamples->Fill(isample,calibpulse[isample]);
            //hsamples->SetBinError(hsamples->FindFixBin(isample),pedrms);
            hsamples->SetBinError(hsamples->FindFixBin(isample),caliberr[isample]);
          }
          
          hfullpulse->Draw("HIST");
          for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
            hists[ipulse]->Draw("HISTSAME");
          }
          hsamples->Draw("ESAME");
          c->SaveAs(TString::Format("pulse_%i_%5f.eps",int(barrel),ootsig));
          delete c;
          delete hsamples;
          delete hfullpulse;
          for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
            delete hists[ipulse];
          }
          //return;
        }        

      }       

    }
    
    //return;
    
    std::vector<int> colors(10);
    colors[0] = 46;
    colors[1] = 4;
    colors[2] = 6;
    colors[3] = 9;
    colors[4] = 8;
    colors[5] = 2;
    colors[6] = kOrange+10;
    colors[7] = 42;
    colors[8] = kViolet+10;
    colors[9] = kTeal+3;
    
    for (const EcalUncalibratedRecHit *uncalibRHp : outuncalib) {
      const EcalUncalibratedRecHit &uncalibRH = *uncalibRHp;
      DetId detid=uncalibRH.id();
      
      bool barrel = detid.subdetId()==EcalBarrel;
      
      const std::vector<double> &calibpulse = calibpulses[detid];
      
      std::set<int> &activeBX = activeBXs[detid];
      activeBX.insert(0);
       
      bool doprint = activeBX.size()>1;         
      //doprint = false;   
      
      const double alpha = barrel ? 1.36745 : 1.59918;
      const double beta = barrel ? 1.51341 : 1.46568;      
      const double tmax = barrel ? 5.04159 : 5.09323;      
      
      const unsigned int npulse = activeBX.size() ;
      const unsigned int nsample = calibpulse.size();
      
      TMatrixD pulsemat(nsample,npulse);
            
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        int ipulse = std::distance(activeBX.begin(),pulsebx);
        for (unsigned int isample=0; isample<nsample; ++isample) {
          double t = double(isample);
          double val = pulseShape(t,alpha,beta,tmax + double(*pulsebx));
          pulsemat[isample][ipulse] = val;
        }
      }
      
      double pedrms = pedrmss[detid];      

      bool status = false;
      
      double shapeerr = 1e-2;
      
      TMatrixDSym invsamplecov = pedrms*pedrms*invsamplecor(barrel,uncalibRH.gain());
      for (unsigned int isample = 0; isample<nsample; ++isample) {
        invsamplecov(isample,isample) += pow(shapeerr*calibpulse[isample],2);
      }
      std::vector<double> caliberr(nsample);
      for (unsigned int isample=0; isample<nsample; ++isample) {
        caliberr[isample] = sqrt(invsamplecov(isample,isample));
      }      
      invsamplecov.Invert();
      PulseChiSq pulsefunc(calibpulse,pulsemat,invsamplecov);
      
      ROOT::Minuit2::Minuit2Minimizer minim;
      //minim.SetStrategy(2);
      //if (doprint) minim.SetPrintLevel(9);
      minim.SetFunction(pulsefunc);
      
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        int ipulse = std::distance(activeBX.begin(),pulsebx);
        minim.SetLowerLimitedVariable(ipulse,pulsenames[ipulse],calibpulse[*pulsebx+5],0.001,0.);
      }
      
      std::vector<bool> floating(activeBX.size(),true);
      int nfloating = activeBX.size();
      while (!status) {        
        status = minim.Minimize();
        
        if (!status && nfloating==1) break;
        
        double minval = std::numeric_limits<double>::max();
        int ipulsemin = 0;
        if (!status) {
          for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
            int ipulse = std::distance(activeBX.begin(),pulsebx);
            if (!floating[ipulse]) continue;
            double pulseval = calibpulse[*pulsebx+5];
            if (pulseval<minval) {
              minval = pulseval;
              ipulsemin = ipulse;
            }
            minim.SetVariableValue(ipulse,pulseval);
          }
          minim.SetFixedVariable(ipulsemin,pulsenames[ipulsemin],0.);
          floating[ipulsemin] = false;
          --nfloating;
        }
      }
      
      double maxpulseval = -std::numeric_limits<double>::max();
      for (unsigned int isample = 0; isample<nsample; ++isample) {
        if (calibpulse[isample]>maxpulseval) {
          maxpulseval = calibpulse[isample];
        }
      }
      
      double maxpulsevaloot = -std::numeric_limits<double>::max();
      int maxbxoot = 0;
      for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
        if ( (*pulsebx)==0 ) continue;
        int ipulse = std::distance(activeBX.begin(),pulsebx);      
        if (minim.X()[ipulse]>maxpulsevaloot) {
          maxpulsevaloot = minim.X()[ipulse];
          maxbxoot = *pulsebx;
        }        
      }
      
      if (!status) {
        printf("final fit failed: barrel = %i, maxpulseval = %5f\n",int(barrel),maxpulseval);
        continue;
      }      
         
      
      const double chisq = minim.MinValue();
         
    
         
      //make rechit from in time pulse
      int ipulseintime = std::distance(activeBX.begin(),activeBX.find(0));
      double energy = minim.X()[ipulseintime];
      //double energy = 0.;
      double time = 0.;      
            
      
      doprint &= (barrel && energy>1.) || (!barrel && energy>5.);
      
      gStyle->SetErrorX(0.);
      gStyle->SetOptStat(0.);
      
      if (doprint) {
        printf("barrel = %i, pederr = %5f, ", int(barrel), uncalibRH.pedrms());
        if (barrel) {
          EBDetId ebid(detid);
          printf("ieta = %i, iphi=%i\n",ebid.ieta(),ebid.iphi());
        }
        else {
          EEDetId eeid(detid);
          printf("ix = %i, iy=%i\n",eeid.ix(),eeid.iy());          
        }
        for (unsigned int isample = 0; isample<calibpulse.size(); ++isample) {
          printf(" %5f",calibpulse[isample]);
        }
        printf("\n");    
        printf("active BXs:");
        for (int bx : activeBX) {
          printf(" %i",bx);
        }
        printf("\n");        
        
        printf("chisq = %5f\n", chisq);
        for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
          int ipulse = std::distance(activeBX.begin(),pulsebx);
          printf("ipulse = %i, bx = %i, energy = %5f +- %5f, sig = %5f\n",ipulse,*pulsebx,minim.X()[ipulse],minim.Errors()[ipulse],minim.X()[ipulse]/minim.Errors()[ipulse]);
        }
        
        bool dodraw = true;
        //dodraw = false;
        if (dodraw) {
          TCanvas *c = new TCanvas;
          TH1D *hsamples = new TH1D("hsamples","",10,-0.5,9.5);
          TH1D *hfullpulse = new TH1D("hfullpulse","",10,-0.5,9.5);
          std::vector<TH1D*> hists;
          for (std::set<int>::const_iterator pulsebx = activeBX.begin(); pulsebx!=activeBX.end(); ++pulsebx) {
            int ipulse = std::distance(activeBX.begin(),pulsebx);
            TH1D *hpulse = new TH1D(TString::Format("hpulse_%i",ipulse),"",10,-0.5,9.5);
            hpulse->SetLineColor(colors[*pulsebx+5]);
            hpulse->SetMarkerColor(colors[*pulsebx+5]);
            hpulse->SetLineWidth(2);
            for (unsigned int isample=0; isample<nsample; ++isample) {
              //double pulseval = pulseShape(double(isample),alpha, beta,tmaxv[ipulse]);
              double pulseval = pulsemat(isample,ipulse);
//               hpulse->Fill(isample,fitvals[ipulse]*pulsemat(isample,ipulse));
//               hfullpulse->Fill(isample,fitvals[ipulse]*pulsemat(isample,ipulse));
              hpulse->Fill(isample,minim.X()[ipulse]*pulseval);
              hfullpulse->Fill(isample,minim.X()[ipulse]*pulseval);              
            }
            hists.push_back(hpulse);
          }
          double minsample = std::numeric_limits<double>::max();
          double minsampleerr = 0.;
          for (unsigned int isample=0; isample<nsample; ++isample) {
            if (calibpulse[isample]<minsample) {
              minsample = calibpulse[isample];
              minsampleerr = caliberr[isample];
            }
            hsamples->Fill(isample,calibpulse[isample]);
            //hsamples->SetBinError(hsamples->FindFixBin(isample),pedrms);
            hsamples->SetBinError(hsamples->FindFixBin(isample),caliberr[isample]);
          }
          double plotmin = std::min(0.,minsample-2.*minsampleerr);
          hfullpulse->SetLineWidth(2);
          hfullpulse->SetMinimum(plotmin);
          hfullpulse->GetXaxis()->SetTitle("Sample");
          hfullpulse->GetYaxis()->SetTitle("Calibrated/Pedestal-Subtracted Energy (GeV)");
          hfullpulse->Draw("HIST");
          for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
            hists[ipulse]->Draw("HISTSAME");
          }
          hsamples->SetMarkerStyle(8);
          hsamples->Draw("ESAME");
          c->SaveAs(TString::Format("pulse_%i_%5f_%i_%5f.eps",int(barrel),energy,maxbxoot,maxpulsevaloot));
          delete c;
          delete hsamples;
          delete hfullpulse;
          for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
            delete hists[ipulse];
          }
          //return;
        }        
        
      }               
      
      EcalRecHit rh( uncalibRH.id(), energy, time );
      //rh.setChi2( uncalibRH.chi2() );
      //rh.setOutOfTimeEnergy( uncalibRH.outOfTimeEnergy() * adcToGeVConstant_ * intercalibConstant );
      /* rh.setOutOfTimeChi2( uncalibRH.outOfTimeChi2() ); */
      //rh.setTimeError(uncalibRH.jitterErrorBits());

      // Now fill flags
      

      float lasercalib = 1.;
      if ( laserCorrection_ )   lasercalib = laser->getLaserCorrection( detid, evt.time());      
      
//           
      if (detid.subdetId() == EcalBarrel && (lasercalib < EBLaserMIN_ || lasercalib > EBLaserMAX_)) 
          rh.setFlag(EcalRecHit::kPoorCalib);
      if (detid.subdetId() == EcalEndcap && (lasercalib < EELaserMIN_ || lasercalib > EELaserMAX_)) 
          rh.setFlag(EcalRecHit::kPoorCalib);      
      
      bool good=true;
      
      if ( uncalibRH.checkFlag(EcalUncalibratedRecHit::kLeadingEdgeRecovered) ){
              rh.setFlag(EcalRecHit::kLeadingEdgeRecovered);
              good=false; 
      } 
      if ( uncalibRH.checkFlag(EcalUncalibratedRecHit::kSaturated) ) { 
              // leading edge recovery failed - still keep the information 
              // about the saturation and do not flag as dead 
              rh.setFlag(EcalRecHit::kSaturated); 
              good=false;
      } 
      if( uncalibRH.isSaturated() ) {
              rh.setFlag(EcalRecHit::kSaturated);
              good=false;
      } 
      if ( uncalibRH.checkFlag(EcalUncalibratedRecHit::kOutOfTime) ) {
              rh.setFlag(EcalRecHit::kOutOfTime) ;
              good=false;   
      } 
      if ( uncalibRH.checkFlag(EcalUncalibratedRecHit::kPoorReco) ) {
              rh.setFlag(EcalRecHit::kPoorReco);
              good=false;
      }
      if ( uncalibRH.checkFlag( EcalUncalibratedRecHit::kHasSwitchToGain6 ) ) {
              rh.setFlag(EcalRecHit::kHasSwitchToGain6);
      }
      if (uncalibRH.checkFlag( EcalUncalibratedRecHit::kHasSwitchToGain1 ) ) {
              rh.setFlag(EcalRecHit::kHasSwitchToGain1);
      }
      
      if (good) rh.setFlag(EcalRecHit::kGood);
      
      result.push_back(rh);
      
    }
    
    
    
    //return true;
}

// // Take our association map of dbstatuses-> recHit flagbits and return the apporpriate flagbit word
// uint32_t EcalRecHitWorkerMulti::setFlagBits(const std::vector<std::vector<uint32_t> >& map, 
// 					     const uint32_t& status  ){
//   
//   for (unsigned int i = 0; i!=map.size(); ++i){
//     if (std::find(map[i].begin(), map[i].end(),status)!= map[i].end()) 
//       return 0x1 << i;
//   }
// 
//   return 0;
// }


EcalRecHitWorkerMulti::~EcalRecHitWorkerMulti(){

  TFile *fpulse = new TFile("fpulse.root","RECREATE");
  hpulseEB->SetDirectory(fpulse);
  hpulseEE->SetDirectory(fpulse);
  hpulseprofEB->SetDirectory(fpulse);
  hpulseprofEE->SetDirectory(fpulse);
  
  hpulseEB->Write();
  hpulseEE->Write();
  hpulseprofEB->Write();
  hpulseprofEE->Write();
  fpulse->Close();
  
  delete rechitMaker_;

  return;  
  
  {
    ROOT::Minuit2::Minuit2Minimizer minim;
    minim.SetPrintLevel(9);
    minim.SetFunction(pulseglobalEB);

    PulseChiSqGlobal &pulseglobal = pulseglobalEB;
    
    unsigned int npulse = pulseglobal.samples().size();
    
    minim.SetLimitedVariable(0,"alphaEB",1.5,0.01,1.,3.);
    minim.SetLimitedVariable(1,"betaEB",1.5,0.01,1.,3.);
    for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
      double pulseval = pulseglobal.samples()[ipulse][5];
      minim.SetLimitedVariable(2+2*ipulse,TString::Format("amp_%i",ipulse).Data(),pulseval,0.01,0.,10.*pulseval);
      minim.SetLimitedVariable(3+2*ipulse,TString::Format("tmax_%i",ipulse).Data(),5.,0.01,4.8,5.2);
    }
    
    bool status = minim.Minimize(); 
    
    double sumtmax = 0.;
    double sumtmaxsq = 0.;
    double sumw = 0;
    for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
      double tmax = minim.X()[3+2*ipulse];
      sumtmax += tmax;
      sumtmaxsq += tmax*tmax;
      sumw += 1.;
    }
    double meantmax = sumtmax/sumw;
    double meantmaxerr = sqrt(sumtmaxsq/sumw - meantmax*meantmax);
    
    printf("status = %i, alpha = %5f +- %5f, beta = %5f +- %5f, tmax = %5f +- %5f\n",int(status),minim.X()[0],minim.Errors()[0], minim.X()[1],minim.Errors()[1],meantmax,meantmaxerr);
    
  }
  
  {
    ROOT::Minuit2::Minuit2Minimizer minim;
    minim.SetPrintLevel(9);
    minim.SetFunction(pulseglobalEE);

    PulseChiSqGlobal &pulseglobal = pulseglobalEE;
    
    unsigned int npulse = pulseglobal.samples().size();
    
    minim.SetLimitedVariable(0,"alphaEE",1.5,0.01,1.,3.);
    minim.SetLimitedVariable(1,"betaEE",1.5,0.01,1.,3.);
    for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
      double pulseval = pulseglobal.samples()[ipulse][5];
      minim.SetLimitedVariable(2+2*ipulse,TString::Format("amp_%i",ipulse).Data(),pulseval,0.01,0.,10.*pulseval);
      minim.SetLimitedVariable(3+2*ipulse,TString::Format("tmax_%i",ipulse).Data(),5.,0.01,4.8,5.2);
    }
    
    bool status = minim.Minimize(); 
    
    double sumtmax = 0.;
    double sumtmaxsq = 0.;
    double sumw = 0;
    for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
      double tmax = minim.X()[3+2*ipulse];
      sumtmax += tmax;
      sumtmaxsq += tmax*tmax;
      sumw += 1.;
    }
    double meantmax = sumtmax/sumw;
    double meantmaxerr = sqrt(sumtmaxsq/sumw - meantmax*meantmax);
    
    printf("status = %i, alpha = %5f +- %5f, beta = %5f +- %5f, tmax = %5f +- %5f\n",int(status),minim.X()[0],minim.Errors()[0], minim.X()[1],minim.Errors()[1],meantmax,meantmaxerr);
    
  }  
  
  

}


double EcalRecHitWorkerMulti::pulseShape(double t, double alpha, double beta, double tmax) const {  
  const double alfabeta = alpha*beta;
  const double tr = t-tmax;  
  
  if(tr > -alfabeta)  {
    double dtsbeta=tr/beta;
    double variable=1.+tr/alfabeta ;
    double puiss=pow(variable,alpha);
    double val = puiss*vdt::fast_exp(-dtsbeta);
    //printf("val = %5f\n",val);
    return val;
  }
  else {
    return 0.;
  }
  
  
}

const TMatrixDSym &EcalRecHitWorkerMulti::invsamplecor(bool barrel, int gain) const {
  if (barrel) {
    if (gain==6) {
      return invsamplecorEBg6;
    }
    else if (gain==1) {
      return invsamplecorEBg1;
    }
    else {
      return invsamplecorEBg12;
    }    
  }
  else {
    if (gain==6) {
      return invsamplecorEEg6;
    }
    else if (gain==1) {
      return invsamplecorEEg1;
    }
    else {
      return invsamplecorEEg12;
    }        
  }
  
  return invsamplecorEBg12;
  
}

EcalRecHitWorkerMulti::PulseChiSq::PulseChiSq(const std::vector<double> &samples, const TMatrixD &pulsemat, const TMatrixDSym &invsamplecov) :
  _sampvec(samples.size(),samples.data()),
  _pulsemat(pulsemat),
  _invsamplecov(invsamplecov)
{
  
}

double EcalRecHitWorkerMulti::PulseChiSq::DoEval(const double *invals) const {
  _workvec.ResizeTo(_pulsemat.GetNcols());
  _workvec.SetElements(invals);
  
  _workvec *= _pulsemat;
  _workvec -= _sampvec;
  
//   double chisq = 0.;
//   int nsample = _sampvec.GetNrows();
//   for (int i=0; i<nsample; ++i) {
//     chisq += _workvec(i)*_workvec(i)*_invsamplecov(i,i);
//     for (int j=0; j<nsample; ++j) {
//       chisq += _workvec(i)*_workvec(j)*_invsamplecov(i,j);
//     }    
// //     for (int j=i+1; j<nsample; ++j) {
// //       chisq += 2.0*_workvec(i)*_workvec(j)*_invsamplecov(i,j);
// //     }
//   }
//   return chisq;
  
  return _invsamplecov.Similarity(_workvec);
  
  //return _workvec.Norm2Sqr();
  
}

EcalRecHitWorkerMulti::PulseChiSqTime::PulseChiSqTime(const std::vector<double> &samples, const TMatrixDSym &invsamplecov, const std::vector<double> &alphanom, const std::vector<double> &betanom, const std::vector<double> &tmaxnom, const TMatrixDSym &invparamcov) :
  _sampvec(samples.size(),samples.data()),
  _invsamplecov(invsamplecov),
  _workvec(samples.size()),
  _paramworkvec(3),
  _alphanom(alphanom),
  _betanom(betanom),
  _tmaxnom(tmaxnom),
  _invparamcov(invparamcov)
{
  
}

double EcalRecHitWorkerMulti::PulseChiSqTime::DoEval(const double *invals) const {

  const unsigned int npulse = _tmaxnom.size();
  const unsigned int nsample = _sampvec.GetNrows();
  
  
  
  _workvec.Zero();
  
  double chisq = 0.;
  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    const double amp = invals[0 + 4*ipulse];
    const double alpha = invals[1 + 4*ipulse];
    const double beta = invals[2 + 4*ipulse];
    const double tmax = invals[3 + 4*ipulse];
    const double alphabeta = alpha*beta;
    
    _paramworkvec[0] = alpha - _alphanom[ipulse];
    _paramworkvec[1] = beta - _betanom[ipulse];
    _paramworkvec[2] = tmax - _tmaxnom[ipulse];
    
    chisq += _invparamcov.Similarity(_paramworkvec);
    
    //chisq += pow(tmax-_tmaxnom[ipulse],2)*_invtmaxerrsq;
    if (amp<=0.) continue;
    for (unsigned int isample = 0; isample<nsample; ++isample) {
      const double tr = double(isample)-tmax;
      _workvec[isample] += (tr > -alphabeta) ? amp*pow(1.+tr/alphabeta,alpha)*vdt::fast_exp(-tr/beta) : 0.;
    }
  }
  _workvec -= _sampvec;
  chisq += _invsamplecov.Similarity(_workvec);

  
  return chisq;
  
}

double EcalRecHitWorkerMulti::PulseChiSqGlobal::DoEval(const double *invals) const {

  const unsigned int npulse = samples_.size();
  const unsigned int nsample = workvec_.GetNrows();
  
  const double alpha = invals[0];
  const double beta = invals[1];
  const double alphabeta = alpha*beta;
  
  double chisq = 0.;
  for (unsigned int ipulse=0; ipulse<npulse; ++ipulse) {
    const double amp = invals[2 + 2*ipulse];
    const double tmax = invals[3 + 2*ipulse];
    const TMatrixDSym &invsamplecov = invcovs_[ipulse];
    const TVectorD &sample = samples_[ipulse];
    for (unsigned int isample = 0; isample<nsample; ++isample) {
      const double tr = double(isample)-tmax;
      workvec_[isample] = (tr > -alphabeta) ? amp*pow(1.+tr/alphabeta,alpha)*exp(-tr/beta) : 0.;
    }
    workvec_ -= sample;
    chisq += invsamplecov.Similarity(workvec_);
  }
  
  return chisq;
  
}

#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalRecHitWorkerFactory, EcalRecHitWorkerMulti, "EcalRecHitWorkerMulti" );
