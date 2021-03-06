#ifndef MMC_PolarizedTau_h
#define MMC_PolarizedTau_h
         
#include "Selection.h"
#include <vector>
#include "TString.h"

class MMC_PolarizedTau : public Selection {

 public:
  MMC_PolarizedTau(TString Name_, TString id_);
  virtual ~MMC_PolarizedTau();

  virtual void  Configure();
  virtual void  Finish();

//  TF1 DeltaR_had(double pt);
//  TF1 DeltaR_lep(double pt);
double DeltaR_had(double pt,double deltaR);
double DeltaR_lep(double pt,double deltaR);
  /* double f_const1 (double *x,double *par);         */
  /* double f_const2 (double *x,double *par);   */
  /* double f_const2_had (double *x,double *par); */

  
  enum cuts {isZtautautopimu=0,NCuts};
  
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  
 private:
  // Selection Variables
  std::vector<TH1D> LongitudinalPolarization;
  std::vector<TH1D> LongitudinalPolarization_Spin;

  std::vector<TH1D> spin_WT;

  //// Gen Taus /////
  std::vector<TH1D> M_TauPi_gen;
  std::vector<TH1D> Pt_TauPi_gen;
  std::vector<TH1D> Eta_TauPi_gen;
  std::vector<TH1D> Phi_TauPi_gen;
  std::vector<TH1D> Theta_TauPi_gen;
  
  std::vector<TH1D> M_TauMu_gen;   
  std::vector<TH1D> Pt_TauMu_gen;
  std::vector<TH1D> Eta_TauMu_gen;
  std::vector<TH1D> Phi_TauMu_gen;
  std::vector<TH1D> Theta_TauMu_gen;

  std::vector<TH1D> Mz_gen;
  std::vector<TH1D> Mz_gen1;

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
 
  /////////////////////  dR  ////////////////////////////////
  //Delta_R(Mu_rec , Mu_gen) , Delta_R(Pi_rec , Pi_gen)
  std::vector<TH1D> dR_recMu_genMu;
  std::vector<TH1D> mindR_recMu_genMu;
  std::vector<TH1D> dR_recPi_genPi;
  std::vector<TH1D> mindR_recPi_genPi;
int mass_dum_count_gen;
  int mass_dum_count;
int sol_gen_count;
int sol_count;  
int n_sel_ev;
// MET infos
  std::vector<TH1D> MetX;
  std::vector<TH1D> MetY;
  std::vector<TH1D> MetX_gen;
  std::vector<TH1D> MetY_gen;
  
  // Di-Tau mass
  std::vector<TH1D> mtautau11_gen;
  std::vector<TH1D> mtautau12_gen;
  std::vector<TH1D> mtautau21_gen;
  std::vector<TH1D> mtautau22_gen;
  std::vector<TH1D> mtautau3_gen; 
  std::vector<TH1D> mtautau4_gen; 
  std::vector<TH1D> mtautau5_gen; 
  std::vector<TH1D> Lw_gen; 
  std::vector<TH1D> mtautau_gen; 

  // Di-Tau mass
  std::vector<TH1D> mtautau11;
  std::vector<TH1D> mtautau12;
  std::vector<TH1D> mtautau21;
  std::vector<TH1D> mtautau22;
  std::vector<TH1D> mtautau3; 
  std::vector<TH1D> mtautau4; 
  std::vector<TH1D> mtautau5; 
  std::vector<TH1D> Lw; 
  std::vector<TH1D> mtautau; 
  // Di Tau mass
  std::vector<TH1D> mtautau_all_1_gen;
  std::vector<TH1D> mtautau_all_2_gen;
  std::vector<TH1D> mtautau_bm_1_gen;
  std::vector<TH1D> mtautau_bm_2_gen;
  std::vector<TH1D> mtautau_bw_1_gen;
  std::vector<TH1D> mtautau_bw_2_gen;

  std::vector<TH1D> mtautau_all_1;
  std::vector<TH1D> mtautau_all_2;
  std::vector<TH1D> mtautau_bm_1;
  std::vector<TH1D> mtautau_bm_2;
  std::vector<TH1D> mtautau_bw_1;
  std::vector<TH1D> mtautau_bw_2;

  std::vector<TH1D> TauE_all_1;
  std::vector<TH1D> Taupt_all_1;
std::vector<TH1D> H_SPIN_OBS_X1;
std::vector<TH1D> H_SPIN_OBS_X2;
std::vector<TH1D> H_SPIN_OBS_X3;
std::vector<TH1D> H_SPIN_OBS_X4;
std::vector<TH1D> H_SPIN_OBS_X1_1;
std::vector<TH1D> H_SPIN_OBS_X1_gen;
std::vector<TH1D> H_SPIN_OBS_X2_gen;
std::vector<TH1D> H_SPIN_OBS_X3_gen;


std::vector<TH1D> Diff_Spin_obs_X1;
std::vector<TH1D> Diff_Spin_obs_X2;
std::vector<TH1D> Diff_Spin_obs_X3;


 std::vector<TH1D> mtautau_all_Res_1;
 std::vector<TH1D>  mtautau_all_Res_2;
 std::vector<TH2D>  mtautau_all_Resp_Matr_1;
 std::vector<TH2D>  mtautau_all_Resp_Matr_2;

 std::vector<TH1D> mtautau_bm_Res_1;
 std::vector<TH1D>  mtautau_bm_Res_2;
 std::vector<TH2D>  mtautau_bm_Resp_Matr_1;
 std::vector<TH2D>  mtautau_bm_Resp_Matr_2;

 std::vector<TH1D> mtautau_bw_Res_1;
 std::vector<TH1D>  mtautau_bw_Res_2;
 std::vector<TH2D>  mtautau_bw_Resp_Matr_1;
 std::vector<TH2D>  mtautau_bw_Resp_Matr_2;

int tt;

  //2D histos
  std::vector<TH2D> Mzgenmmc_gen;
  std::vector<TH2D> Mzrecmmc_gen;
  std::vector<TH2D> Mzrecmmc_genmmc;


   bool mmc_method(const TLorentzVector &p4vis1,const TLorentzVector &p4vis2,const TVector2 MET,const TVector2 MET_sig,
double &mass_z,double &tau_pt, double &Spin_Obs_X1, double &Spin_Obs_X2,double &Spin_Obs_X3,double &Spin_Obs_X4); 
  bool verbose; 
    
};
#endif
   
