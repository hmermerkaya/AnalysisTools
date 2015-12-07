#ifndef Ztotautau_pi_spin_h
#define Ztotautau_pi_spin_h

#include "Selection.h"
#include <vector>
#include "TString.h"

class Ztotautau_pi_spin : public Selection {

 public:
  Ztotautau_pi_spin(TString Name_, TString id_);
  virtual ~Ztotautau_pi_spin();

  virtual void  Configure();

  enum cuts {TriggerOk=0, 
	     hasTag,
	     TagPtmin,
             TagPtmax,
	     TagIso,
	     NJets,
	     MaxTracksinJet,
	     MinTracksinJet,
	     charge,
	     deltaPhi,
	     MT,
	     MET,
	     JetPt,
	     etaq,
	     JetTrackPtMax,
	     ZMassmax,
             ZMassmin,
	     HT,
	     NCuts};

  enum Channel{muontag,electontag,rhotag,threepiontag,NChannels};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH2D> TagEtaPT;

  std::vector<TH1D> TauCandFound;
  std::vector<TH2D> TauCandEtaPhi;

  std::vector<TH1D> TauCandPhi;
  std::vector<TH1D> TauCandPhiRes;
  std::vector<TH1D> TauCandEta;
  std::vector<TH1D> TauCandEtaRes;
  std::vector<TH1D> TauCandE;
  std::vector<TH1D> TauCandERes;

  std::vector<TH1D> TauCandP;
  std::vector<TH1D> TauCandPRes;
  std::vector<TH1D> TauCandPT;
  std::vector<TH1D> TauCandPTRes;

  std::vector<TH1D> TauCandMass;
  std::vector<TH1D> TauCandPiMass;
  std::vector<TH1D> TauCandNuMass;

  std::vector<TH1D> TauSolutionResult;
  std::vector<TH1D> EstimatedTauE;
  std::vector<TH1D> EstimatedTauPhi;
  std::vector<TH1D> EstimatedTauEta;

  std::vector<TH1D> EstimatedTauERes;
  std::vector<TH1D> EstimatedTauPhiRes;
  std::vector<TH1D> EstimatedTauEtaRes;

  std::vector<TH1D> EstimatedTauDirPhiRes;
  std::vector<TH1D> EstimatedTauDirEtaRes;

  std::vector<TH1D> KFTau_Fit_chiprob;
  std::vector<TH1D> KFTau_Fit_a1mass;
  std::vector<TH1D> KFTau_Fit_chi2;
  std::vector<TH1D> KFTau_Fit_ndf;
  std::vector<TH1D> KFTau_Fit_ambiguity;
  std::vector<TH1D> KFTau_Fit_csum;
  std::vector<TH1D> KFTau_Fit_iterations;
  std::vector<TH1D> KFTau_Fit_TauEnergyFraction;
  std::vector<TH1D> KFTau_Fit_PV_PV_significance;
  std::vector<TH1D> KFTau_Fit_SV_PV_significance;

  std::vector<TH1D> MaxTauTrackPt;
  std::vector<TH1D> MinTauTrackPt;
  std::vector<TH1D> TauTrackdphitheta;

  std::vector<TH1D>  PFTauCandPhi; 
  std::vector<TH1D>  PFMTauCandPhi;
  std::vector<TH1D>  KFTauCandPhi; 
  std::vector<TH1D>  KFSTauCandPhi;
  std::vector<TH1D>  KFTauVCandPhi;
  std::vector<TH1D>  KFSTauVCandPhi;


  int channel;
  double jeteta,muoneta,TauTrackPtThreshold;

};
#endif
