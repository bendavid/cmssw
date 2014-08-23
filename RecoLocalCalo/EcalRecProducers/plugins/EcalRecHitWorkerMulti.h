#ifndef RecoLocalCalo_EcalRecProducers_EcalRecHitWorkerMulti_hh
#define RecoLocalCalo_EcalRecProducers_EcalRecHitWorkerMulti_hh

/** \class EcalRecHitSimpleAlgo
  *  Simple algoritm to make rechits from uncalibrated rechits
  *
  *  \author Shahram Rahatlou, University of Rome & INFN, March 2006
  */

#include "RecoLocalCalo/EcalRecProducers/interface/EcalRecHitWorkerBaseClass.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalRecHitSimpleAlgo.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTimeOffsetConstant.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "CondFormats/EcalObjects/interface/EcalTimeBiasCorrections.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "Math/IFunction.h"

#include "TH1D.h"
#include "TProfile.h"

class EcalRecHitWorkerMulti : public EcalRecHitWorkerBaseClass {
        public:
                EcalRecHitWorkerMulti(const edm::ParameterSet&, edm::ConsumesCollector& c);
				EcalRecHitWorkerMulti(const edm::ParameterSet&);
                virtual ~EcalRecHitWorkerMulti();                       
        
                void set(const edm::EventSetup& es);
                void run(const edm::Event& evt, const EcalUncalibratedRecHitCollection& uncalibRHs, EcalRecHitCollection & result);



        protected:

		double EBLaserMIN_;
		double EELaserMIN_;
		double EBLaserMAX_;
		double EELaserMAX_;


        edm::ESHandle<EcalIntercalibConstants> ical;
        edm::ESHandle<EcalTimeCalibConstants> itime;
        edm::ESHandle<EcalTimeOffsetConstant> offtime;
        edm::ESHandle<EcalTimeBiasCorrections> timeCorrBias;        
        edm::ESHandle<EcalADCToGeVConstant> agc;
        edm::ESHandle<EcalChannelStatus> chStatus;
        std::vector<int> v_chstatus_;
        edm::ESHandle<EcalLaserDbService> laser;
        edm::ESHandle<CaloTopology> theCaloTopology;
        edm::ESHandle<CaloGeometry> caloGeometry;


		// Associate reco flagbit ( outer vector) to many db status flags (inner vector)
        std::vector<int> v_DB_reco_flags_;

// 		uint32_t setFlagBits(const std::vector<std::vector<uint32_t> >& map, 
// 				     const uint32_t& status  );

        double pulseShape(double t, double alpha, double beta, double tmax) const;
        const TMatrixDSym &invsamplecor(bool barrel, int gain) const;
                
//         uint32_t flagmask_; // do not propagate channels with these flags on

        bool killDeadChannels_;
        bool laserCorrection_;
        bool blindtagging_;
        
        EcalRecHitSimpleAlgo * rechitMaker_;
        
        TMatrixDSym invsamplecorEBg12;
        TMatrixDSym invsamplecorEEg12;
        TMatrixDSym invsamplecorEBg6;
        TMatrixDSym invsamplecorEEg6;
        TMatrixDSym invsamplecorEBg1;
        TMatrixDSym invsamplecorEEg1;        
        
        TMatrixDSym sampleunitcov;
        
        TH1D *hpulseEB;
        TH1D *hpulseEE;
        
        TProfile *hpulseprofEB;
        TProfile *hpulseprofEE;
        

        
        class PulseChiSq : public ROOT::Math::IBaseFunctionMultiDim {
          public:
            PulseChiSq(const std::vector<double> &samples, const TMatrixD &pulsemat, const TMatrixDSym &invsamplecov);
            unsigned int NDim() const { return _pulsemat.GetNcols(); }
            IBaseFunctionMultiDim *Clone() const { return new PulseChiSq(*this); }
            
          protected:
            TVectorD _sampvec;
            TMatrixD _pulsemat;
            TMatrixDSym _invsamplecov;
            mutable TVectorD _workvec;
            
          private:
            double DoEval(const double *invals) const;
        };
        
        class PulseChiSqTime : public ROOT::Math::IBaseFunctionMultiDim {
          public:
            PulseChiSqTime(const std::vector<double> &samples, const TMatrixDSym &invsamplecov, const std::vector<double> &alphanom, const std::vector<double> &betanom, const std::vector<double> &tmaxnom, const TMatrixDSym &invparamcov);
            unsigned int NDim() const { return 2*_tmaxnom.size(); }
            IBaseFunctionMultiDim *Clone() const { return new PulseChiSqTime(*this); }
            
          protected:
            TVectorD _sampvec;
            TMatrixDSym _invsamplecov;
            mutable TVectorD _workvec;
            mutable TVectorD _paramworkvec;
            std::vector<double> _alphanom;
            std::vector<double> _betanom;
            std::vector<double> _tmaxnom;
            TMatrixDSym _invparamcov;
            
          private:
            double DoEval(const double *invals) const;
        };        
        
        class PulseChiSqGlobal : public ROOT::Math::IBaseFunctionMultiDim {
          public:
            PulseChiSqGlobal() : workvec_(10) {}
            unsigned int NDim() const { return 2 + 2*samples_.size(); }
            IBaseFunctionMultiDim *Clone() const { return new PulseChiSqGlobal(*this); }
            
            std::vector<TVectorD> &samples() { return samples_; }
            std::vector<TMatrixDSym> &invcovs() { return invcovs_; }
            
          protected:
            std::vector<TVectorD> samples_;
            std::vector<TMatrixDSym> invcovs_;
            mutable TVectorD workvec_;
            
          private:
            double DoEval(const double *invals) const;
        };       
        
        PulseChiSqGlobal pulseglobalEB;
        PulseChiSqGlobal pulseglobalEE;        

};

#endif
