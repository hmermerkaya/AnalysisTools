#ifndef Tau_SpinExample_h
#define Tau_SpinExample_h
      
#include "Selection.h"
#include <vector>
#include "TString.h"

class Tau_SpinExample : public Selection {

 public:
  Tau_SpinExample(TString Name_, TString id_);
  virtual ~Tau_SpinExample();

  virtual void  Configure();
  virtual void  Finish();
  enum cuts {isZtautautopimu=0,NCuts};


//std::cout<<"hello";
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  std::vector<TH1D> LongitudinalPolarization;
  std::vector<TH1D> LongitudinalPolarization_Spin;
  ////// Gen-Infos /////// 
 /////// Lab_Frame ///////
  // X1 = Pi_Pt / Tau_Pt
  std::vector<TH1D> pi_PtRatio;  
  std::vector<TH1D> pi_PtRatio_hplus;
  std::vector<TH1D> pi_PtRatio_hminus;
  std::vector<TH1D> pi_PtRatio_Spin;
  std::vector<TH1D> pi_PtRatio_hplus_Spin;
  std::vector<TH1D> pi_PtRatio_hminus_Spin;
  // X2 = Pi_Pt / Z_Mt
  std::vector<TH1D> pi_PtoverZmt;
  std::vector<TH1D> pi_PtoverZmt_hplus;
  std::vector<TH1D> pi_PtoverZmt_hminus;
  std::vector<TH1D> pi_PtoverZmt_Spin;
  std::vector<TH1D> pi_PtoverZmt_hplus_Spin;
  std::vector<TH1D> pi_PtoverZmt_hminus_Spin;
  
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  std::vector<TH1D> pi_EoverEtau;
  std::vector<TH1D> pi_EoverEtau_hplus;
  std::vector<TH1D> pi_EoverEtau_hminus;
  std::vector<TH1D> pi_EoverEtau_Spin;
  std::vector<TH1D> pi_EoverEtau_hplus_Spin;
  std::vector<TH1D> pi_EoverEtau_hminus_Spin;
  // X4 = Pi_E / Z_M
  std::vector<TH1D> pi_EoverZm;
  std::vector<TH1D> pi_EoverZm_hplus;
  std::vector<TH1D> pi_EoverZm_hminus;
  std::vector<TH1D> pi_EoverZm_Spin;
  std::vector<TH1D> pi_EoverZm_hplus_Spin;
  std::vector<TH1D> pi_EoverZm_hminus_Spin;  

  ////// Rec-Infos /////// 
 /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  std::vector<TH1D> pirec_PtRatio;  
  std::vector<TH1D> pirec_PtRatio_hplus;
  std::vector<TH1D> pirec_PtRatio_hminus;
  std::vector<TH1D> pirec_PtRatio_Spin;
  std::vector<TH1D> pirec_PtRatio_hplus_Spin;
  std::vector<TH1D> pirec_PtRatio_hminus_Spin;
  // X2 = Pirec_Pt / Zrec_Mt
  std::vector<TH1D> pirec_PtoverZmt;
  std::vector<TH1D> pirec_PtoverZmt_hplus;
  std::vector<TH1D> pirec_PtoverZmt_hminus;
  std::vector<TH1D> pirec_PtoverZmt_Spin;
  std::vector<TH1D> pirec_PtoverZmt_hplus_Spin;
  std::vector<TH1D> pirec_PtoverZmt_hminus_Spin;
  
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  std::vector<TH1D> pirec_EoverEtau;
  std::vector<TH1D> pirec_EoverEtau_hplus;
  std::vector<TH1D> pirec_EoverEtau_hminus;
  std::vector<TH1D> pirec_EoverEtau_Spin;
  std::vector<TH1D> pirec_EoverEtau_hplus_Spin;
  std::vector<TH1D> pirec_EoverEtau_hminus_Spin;
  // X4 = Pirec_E / Zrec_M
  std::vector<TH1D> pirec_EoverZm;
  std::vector<TH1D> pirec_EoverZm_hplus;
  std::vector<TH1D> pirec_EoverZm_hminus;
  std::vector<TH1D> pirec_EoverZm_Spin;
  std::vector<TH1D> pirec_EoverZm_hplus_Spin;
  std::vector<TH1D> pirec_EoverZm_hminus_Spin;  

  //// Rec Taus /////
 std::vector<TH1D> M_TauPi_rec;
  std::vector<TH1D> Pt_TauPi_rec;
  std::vector<TH1D> Eta_TauPi_rec;
  std::vector<TH1D> Phi_TauPi_rec;
  std::vector<TH1D> Theta_TauPi_rec;
  
  std::vector<TH1D> M_TauMu_rec;   
  std::vector<TH1D> Pt_TauMu_rec;
  std::vector<TH1D> Eta_TauMu_rec;
  std::vector<TH1D> Phi_TauMu_rec;
  std::vector<TH1D> Theta_TauMu_rec;

  /////////////////////   dR  ////////////////////////////////
  //Delta_R(Mu_rec , Mu_gen) , Delta_R(Pi_rec , Pi_gen)
  std::vector<TH1D> dR_recMu_genMu;
  std::vector<TH1D> mindR_recMu_genMu;
  std::vector<TH1D> dR_recPi_genPi;
  std::vector<TH1D> mindR_recPi_genPi;

  //////////// some other histos
  std::vector<TH1D> mu_ExoverEtau;
  std::vector<TH1D> mu_ExoverEtau_hplus;
  std::vector<TH1D> mu_ExoverEtau_hminus;
 
  std::vector<TH1D> pi_ExoverEtau;
  std::vector<TH1D> pi_ExoverEtau_hplus;
  std::vector<TH1D> pi_ExoverEtau_hminus;


  bool verbose; 

  double Zstoa(double zs);

  int zsbins;
  float zsmin,zsmax;


};
#endif
   
