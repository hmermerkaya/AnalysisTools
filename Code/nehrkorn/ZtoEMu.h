#ifndef ZtoEMu_h
#define ZtoEMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

class ZtoEMu : public Selection {

 public:
  ZtoEMu(TString Name_, TString id_);
  virtual ~ZtoEMu();

  virtual void Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
		 mll,
		 drEMu,
		 diMuonVeto,
		 triLeptonVeto,
		 charge,
		 jetVeto,
		 MtMu,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> RelIsoE;
  std::vector<TH1D> RelIsoMu;
  std::vector<TH1D> EPt;
  std::vector<TH1D> EEt;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> mtMu;
  std::vector<TH1D> mtE;
  std::vector<TH1D> etaMu;
  std::vector<TH1D> etaE;
  std::vector<TH1D> jetsum;
  std::vector<TH1D> NJets;
  std::vector<TH1D> NJetsLoose;
  std::vector<TH1D> NJetsMedium;
  std::vector<TH1D> NJetsTight;
  std::vector<TH1D> PUJetId;
  std::vector<TH1D> chargesum;
  std::vector<TH1D> drmue;
  std::vector<TH1D> deltaphi;
  std::vector<TH1D> ptbal;
  std::vector<TH1D> chargesumsigned;
  std::vector<TH1D> FirstJetPt;
  std::vector<TH1D> SecondJetPt;
  
  std::vector<TH1D> invmass_zmass;
  std::vector<TH1D> invmass_ptbalance;
  std::vector<TH1D> invmass_mtmu;
  std::vector<TH1D> invmass_jetveto;
  std::vector<TH1D> invmass_vetos;
  std::vector<TH1D> invmass_only_object_id;
  
  std::vector<TH1D> invmass_zmass_m;
  std::vector<TH1D> invmass_ptbalance_m;
  std::vector<TH1D> invmass_mtmu_m;
  std::vector<TH1D> invmass_jetveto_m;
  std::vector<TH1D> invmass_vetos_m;
  std::vector<TH1D> invmass_only_object_id_m;

  std::vector<TH1D> invmass_dremu_only;
  std::vector<TH1D> invmass_dimuon_only;
  std::vector<TH1D> invmass_trilepton_only;
  std::vector<TH1D> invmass_charge_only;
  std::vector<TH1D> invmass_jetveto_only;
  std::vector<TH1D> invmass_mtmu_only;
  std::vector<TH1D> invmass_ptbal_only;

  std::vector<TH1D> nm0_met;
  std::vector<TH1D> nm0_jetsum;
  std::vector<TH1D> nm0_onejet;
  std::vector<TH1D> nm0_mtmu;
  std::vector<TH1D> nm0_ptbalance;
  
  std::vector<TH1D> NPV;
  std::vector<TH1D> NPV3d;
  std::vector<TH1D> NPVfine;
  std::vector<TH1D> evtweight;
  
  std::vector<TH1D> met;
  std::vector<TH1D> met_uncorr;
  std::vector<TH1D> onejet;
  std::vector<TH1D> mte_mtmu;
  std::vector<TH1D> NbJets;
  std::vector<TH1D> NbJetsVtxL;
  std::vector<TH1D> NbJetsVtxM;
  std::vector<TH1D> NbJetsVtxT;

  // binning tests
  std::vector<TH1D> etaE_offBins;
  std::vector<TH1D> etaE_manyBins;
  std::vector<TH1D> etaMu_offBins;
  std::vector<TH1D> etaMu_manyBins;

  // cross checks
  std::vector<TH1D> mtmu_metgr30;
  std::vector<TH1D> mtmu_metsm30;
  std::vector<TH1D> jet1E;
  std::vector<TH1D> jet2E;
  std::vector<TH1D> jet1Mu;
  std::vector<TH1D> jet2Mu;

  // comparison of generators

  std::vector<TH1D> zpt;
  std::vector<TH1D> zeta;
  std::vector<TH1D> zmass;
  std::vector<TH1D> leadingjet_pt;
  std::vector<TH1D> subleadingjet_pt;
  std::vector<TH1D> leadingjet_eta;
  std::vector<TH1D> subleadingjet_eta;
  std::vector<TH1D> jetsumcustom;

  std::vector<TH1D> ptbal_chargepass;
  std::vector<TH1D> ptbal_chargefail;

  std::vector<TH1D> Dxy_trig;
  std::vector<TH1D> Dz_trig;
  std::vector<TH1D> Dxy_nontrig;
  std::vector<TH1D> Dz_nontrig;
  std::vector<TH1D> Dxy_trignoip;
  std::vector<TH1D> Dz_trignoip;
  std::vector<TH2D> eta_mu_e;
  std::vector<TH2D> pt_vs_eta_mu;
  std::vector<TH2D> pt_vs_eta_e;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax,mtmu,ptbalance,mmin;
  int n_mu,n_e;
  bool doHiggsObjects;
  bool doWWObjects;
  bool useMadgraphZ;
  
  double csvl,csvm,csvt;

  double calculatePzeta(int muiterator, int eiterator);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  bool matchTrigger(unsigned int i, double dr, std::string trigger, std::string object);
  int matchTruth(TLorentzVector tvector);
  bool matchTruth(TLorentzVector tvector, int pid, double dr);
  int findBin(TGraphAsymmErrors* graph, double xval);
  
  bool isTightMuon(unsigned int idx);
  bool isTightMuon(unsigned int idx, unsigned int vtx);
  bool isHiggsMuon(unsigned int idx, unsigned int vtx);
  bool isLooseMuon(unsigned int idx);
  bool isFakeMuon(unsigned int idx);
  bool isFakeMuon(unsigned int idx, unsigned int vtx);
  double Muon_RelIso(unsigned int idx);
  
  bool isTrigPreselElectron(unsigned int idx);
  bool isTrigNoIPPreselElectron(unsigned int idx);
  bool isMVATrigElectron(unsigned int idx);
  bool isMVATrigNoIPElectron(unsigned int idx);
  bool isMVANonTrigElectron(unsigned int idx, unsigned int vtx);
  bool isHiggsElectron(unsigned int idx, unsigned int vtx);
  bool isWWElectron(unsigned int idx, unsigned int vtx);
  bool isTightElectron(unsigned int idx);
  bool isTightElectron(unsigned int idx, unsigned int vtx);
  bool isLooseElectron(unsigned int idx);
  bool isFakeElectron(unsigned int idx);
  bool isFakeElectron(unsigned int idx, unsigned int vtx);
  double Electron_RelIso(unsigned int idx);
  double Electron_Aeff_R04(double Eta);
  double Electron_Aeff_R03(double Eta);
  
  double MuonIDeff(unsigned int idx);
  double MuonIDerrUp(unsigned int idx);
  double MuonIDerrDown(unsigned int idx);
  double MuonHiggsIDeff(unsigned int idx);
  double MuonTriggerEff(unsigned int idx);
  double MuonTriggerErr(unsigned int idx);
  double ElectronIDeff(unsigned int idx, std::string id);
  double ElectronIDerr(unsigned int idx, std::string id);
  double ElectronTrigIDeff(unsigned int idx);
  double ElectronTrigIDerr(unsigned int idx);
  double ElectronNonTrigIDeff(unsigned int idx);
  double ElectronNonTrigIDerr(unsigned int idx);
  double ElectronHiggsIDeff(unsigned int idx);
  double ElectronTriggerEff(unsigned int idx);
  double ElectronTriggerErr(unsigned int idx);
  double ElectronEmbeddedEff(unsigned int idx);
  
  double TriggerEff(unsigned int muid, unsigned int eid, TString path);
  double SingleEle(unsigned int idx);
  double DoubleEleLeading(unsigned int idx);
  double DoubleEleTrailing(unsigned int idx);
  double SingleMu(unsigned int idx);
  double DoubleMuLeading(unsigned int idx);
  double DoubleMuTrailing(unsigned int idx);

  double TrackingEff(double eta);

  double ElectronMassScale(unsigned int idx);
  double ZPtReweight(double zpt);
  double PowhegReweight(double zpt);
  double rundependentJetPtCorrection(double jeteta, int runnumber);

  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  double FakerateWW(unsigned int idx, std::string type);
  double FakerateWWerror(unsigned int idx, std::string type);
  
  TFile* FRFile;
  TFile* EmbEffFile;
  TFile* ZptCorrFile;
  TH1D* ZptCorrection;
  TH2D* ElectronFakeRate;
  TH2D* MuonFakeRate;
  TH2D* EmbEff;
  double fakeRate;
  double fakeRateMu;
  double fakeRateE;
  
  TFile* MuIdEffFile;
  TFile* MuIsoEffFile;
  TFile* ETrigIdEffFile;
  TFile* ENonTrigIdEffFile;
  TFile* TriggerEfficiencies;
  TFile* FakeRates;
  TFile* ENonTrigIdRecoEffFile;

  TH2D* ElectronTrigEff;
  TH2D* ElectronNonTrigEff;
  TH2D* ElectronNonTrigRecoEff;
  TGraphAsymmErrors* MuIdEff09;
  TGraphAsymmErrors* MuIdEff12;
  TGraphAsymmErrors* MuIdEff21;
  TGraphAsymmErrors* MuIdEff24;
  TGraphAsymmErrors* MuIsoEff09;
  TGraphAsymmErrors* MuIsoEff12;
  TGraphAsymmErrors* MuIsoEff21;
  TGraphAsymmErrors* MuIsoEff24;

  /*TGraphAsymmErrors* SingleEle15;
  TGraphAsymmErrors* SingleEle25;
  TGraphAsymmErrors* DoubleEleLead15;
  TGraphAsymmErrors* DoubleEleLead25;
  TGraphAsymmErrors* DoubleEleTrail15;
  TGraphAsymmErrors* DoubleEleTrail25;
  TGraphAsymmErrors* SingleMu08;
  TGraphAsymmErrors* SingleMu12;
  TGraphAsymmErrors* SingleMu21;
  TGraphAsymmErrors* SingleMu25;
  TGraphAsymmErrors* DoubleMuLead12;
  TGraphAsymmErrors* DoubleMuLead21;
  TGraphAsymmErrors* DoubleMuLead25;
  TGraphAsymmErrors* DoubleMuTrail12;
  TGraphAsymmErrors* DoubleMuTrail21;
  TGraphAsymmErrors* DoubleMuTrail25;*/

  TH1D* SingleEle15;
  TH1D* SingleEle25;
  TH1D* DoubleEleLead15;
  TH1D* DoubleEleLead25;
  TH1D* DoubleEleTrail15;
  TH1D* DoubleEleTrail25;
  TH1D* SingleMu08;
  TH1D* SingleMu12;
  TH1D* SingleMu21;
  TH1D* SingleMu25;
  TH1D* DoubleMuLead12;
  TH1D* DoubleMuLead21;
  TH1D* DoubleMuLead25;
  TH1D* DoubleMuTrail12;
  TH1D* DoubleMuTrail21;
  TH1D* DoubleMuTrail25;

  TGraphAsymmErrors* EleFake1;
  TGraphAsymmErrors* EleFake15;
  TGraphAsymmErrors* EleFake2;
  TGraphAsymmErrors* EleFake25;
  TGraphAsymmErrors* MuFake1;
  TGraphAsymmErrors* MuFake15;
  TGraphAsymmErrors* MuFake2;
  TGraphAsymmErrors* MuFake25;

  TF1* gause;
  TF1* gausmu;
  TH1D* eres;
  TH1D* mures;
  double eleres;
  double muonres;

};
#endif

