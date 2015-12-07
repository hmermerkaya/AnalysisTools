#include "MMC_PolarizedTau.h"  
#include "TLorentzVector.h"   
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>  
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
    
#include <iostream>
#include "TError.h"
#include "TH1D.h"
#include <TH2D.h>
#include <fstream>
#include <string>  
#include <sstream>  
#include <iomanip>
#include <time.h>
#include "TFile.h"
#include "TString.h"
#include "TRandom2.h"
#include <algorithm>  
#include <vector>
#include <map>
#include <cmath>
#include <iterator>  
#include "TMath.h"
#include "TCanvas.h"
#define pi 3.14159

using namespace TMath;
using namespace std;


bool Isnan(double var)
{
  volatile double d = var;

  return d != d;
}  

bool NaN(double l)
{
  volatile double b = l;
  return b != b;
}  


double f_const1 (double *x,double *par);
double f_const2 (double *x,double *par);
double f_const2_had (double *x,double *par);
  
double deltaphi(double x){
  while (x >= 3.14159) x -= 2*3.14159;
  while (x < -3.14159) x += 2*3.14159;
  return x; }

typedef std::pair<double , double > MyPair;

struct CompareByKey {
  bool operator() (const MyPair& a, const MyPair& b) const {
    return a.first > b.first;
  };
};
struct CompareByValue {
  bool operator() (const MyPair& a, const MyPair& b) const {
    return a.second > b.second;
  };
};

float fastSin( float x )
{
    const float B = 4.0f/pi;
    const float C = -4.0f/(pi*pi);

    float y = B * x + C * x * abs(x);

    const float P = 0.225f;

    return P * (y * abs(y) - y) + y;
}

float fastCos(float x)
{
    return fastSin(x + (pi / 2));
}


MMC_PolarizedTau::MMC_PolarizedTau(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}     

MMC_PolarizedTau::~MMC_PolarizedTau() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "MMC_PolarizedTau::~MMC_PolarizedTau Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "MMC_PolarizedTau::~MMC_PolarizedTau()" << std::endl;
}
      
void  MMC_PolarizedTau::Configure() {
  // Setup Cut Values
 
 mass_dum_count_gen=0;
 mass_dum_count=0;
sol_gen_count=0,sol_count=0;
n_sel_ev=0;
 for(int i=0; i<NCuts;i++) {
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isZtautautopimu)          cut.at(isZtautautopimu)=1;
  }

  TString hlabel;
  TString htitle; for(int i=0; i<NCuts; i++) {
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
  
    if(i==isZtautautopimu) {
      title.at(i)="Is $Z\\rightarrow\\tau\\tau\\rightarrow\\mu\\pi$ MC (bool)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Is Z#rightarrow#tau#tau#rightarrow#mu#pi MC (bool)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isZtautautopimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isZtautautopimu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms

  /////// Tau_pi --> Pi + Nu //////  
  spin_WT=HConfig.GetTH1D(Name+"_spin_WT","spin_WT",20,0.0,2.1,"spin_WT","Events"); 

  //////////////////////////////////////
    
  // Setup Extra Histograms
  // Reco Objects
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  mindR_recMu_genMu=HConfig.GetTH1D(Name+"_mindR_recMu_genMu","mindR_recMu_genMu",100,0.,1.,"mindR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  mindR_recPi_genPi=HConfig.GetTH1D(Name+"_mindR_recPi_genPi","mindR_recPi_genPi",100,0.,1.,"mindR_recPi_genPi","Events");

  //Z---> TauPigen                                                                                                                                           
  M_TauPi_gen=HConfig.GetTH1D(Name+"_M_TauPi_gen","M_TauPi_gen",180,0.,4.,"M_TauPi_gen","Events");
  Pt_TauPi_gen=HConfig.GetTH1D(Name+"_Pt_TauPi_gen","Pt_TauPi_gen",30,0.,120.,"Pt_TauPi_gen","Events");
  Eta_TauPi_gen=HConfig.GetTH1D(Name+"_Eta_TauPi_gen","Eta_TauPi_gen",50,-2.5,2.5,"Eta_TauPi_gen","Events");
  Phi_TauPi_gen=HConfig.GetTH1D(Name+"_Phi_TauPi_gen","Phi_TauPi_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_gen","Events");
  //Z---> TauMugen                                                                                                                                           
  M_TauMu_gen=HConfig.GetTH1D(Name+"_M_TauMu_gen","M_TauMu_gen",180,0.,4.,"M_TauMu_gen","Events");
  Pt_TauMu_gen=HConfig.GetTH1D(Name+"_Pt_TauMu_gen","Pt_TauMu_gen",30,0.,120.,"Pt_TauMu_gen","Events");
  Eta_TauMu_gen=HConfig.GetTH1D(Name+"_Eta_TauMu_gen","Eta_TauMu_gen",50,-2.5,2.5,"Eta_TauMu_gen","Events");
  Phi_TauMu_gen=HConfig.GetTH1D(Name+"_Phi_TauMu_gen","Phi_TauMu_gen",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_gen","Events");

  Mz_gen=HConfig.GetTH1D(Name+"_Mz_gen", "Mz_gen", 150.,0.,150.,"Mz_{gen}","Events");
  Mz_gen1=HConfig.GetTH1D(Name+"_Mz_gen1", "Mz_gen1", 150.,0.,150.,"Mz_{gen1}","Events");
  //Z---> TauPirec
  M_TauPi_rec=HConfig.GetTH1D(Name+"_M_TauPi_rec","M_TauPi_rec",180,0.,4.,"M_TauPi_rec","Events");
  Pt_TauPi_rec=HConfig.GetTH1D(Name+"_Pt_TauPi_rec","Pt_TauPi_rec",30,0.,120.,"Pt_TauPi_rec","Events");
  Eta_TauPi_rec=HConfig.GetTH1D(Name+"_Eta_TauPi_rec","Eta_TauPi_rec",50,-2.5,2.5,"Eta_TauPi_rec","Events");
  Phi_TauPi_rec=HConfig.GetTH1D(Name+"_Phi_TauPi_rec","Phi_TauPi_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_rec","Events");
  //Z---> TauMurec
  M_TauMu_rec=HConfig.GetTH1D(Name+"_M_TauMu_rec","M_TauMu_rec",180,0.,4.,"M_TauMu_rec","Events");
  Pt_TauMu_rec=HConfig.GetTH1D(Name+"_Pt_TauMu_rec","Pt_TauMu_rec",30,0.,120.,"Pt_TauMu_rec","Events");
  Eta_TauMu_rec=HConfig.GetTH1D(Name+"_Eta_TauMu_rec","Eta_TauMu_rec",50,-2.5,2.5,"Eta_TauMu_rec","Events");
  Phi_TauMu_rec=HConfig.GetTH1D(Name+"_Phi_TauMu_rec","Phi_TauMu_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_rec","Events");
 
  // MET infos
  MetX=HConfig.GetTH1D(Name+"_MetX","MetX",30,-100.,120.,"E_{X}^{miss}","Events");
  MetY=HConfig.GetTH1D(Name+"_MetY","MetY",30,-100.,120.,"E_{Y}^{miss}","Events");
  MetX_gen=HConfig.GetTH1D(Name+"_MetX_gen","MetX_gen",30,-100.,120.,"E_{X_{gen}}^{miss}","Events");
  MetY_gen=HConfig.GetTH1D(Name+"_MetY_gen","MetY_gen",30,-100.,120.,"E_{Y_{gen}}^{miss}","Events");   

  // Di-Tau mass
  mtautau11_gen=HConfig.GetTH1D(Name+"_ditau11_gen", "ditau11_gen", 150.,0.,150.,"M_tautau_11_gen","Events");
  mtautau12_gen=HConfig.GetTH1D(Name+"_ditau12_gen", "ditau12_gen", 150.,0.,150.,"M_tautau_12_gen","Events");
  mtautau21_gen=HConfig.GetTH1D(Name+"_ditau21_gen", "ditau21_gen", 150.,0.,150.,"M_tautau_21_gen","Events");
  mtautau22_gen=HConfig.GetTH1D(Name+"_ditau22_gen", "ditau22_gen", 150.,0.,150.,"M_tautau_22_gen","Events");
  mtautau3_gen=HConfig.GetTH1D(Name+"_ditau3_gen", "ditau3_gen", 150.,0.,150.,"M_tautau_3_gen","Events");
  mtautau4_gen=HConfig.GetTH1D(Name+"_ditau4_gen", "ditau4_gen", 150.,0.,150.,"M_tautau_4_gen","Events");
  mtautau5_gen=HConfig.GetTH1D(Name+"-ditau5_gen", "ditau5_gen", 150.,0.,150.,"M_tautau_5_gen","Events");

  Lw_gen=HConfig.GetTH1D(Name+"-Lw_gen", "Lw_gen", 150.,-5.,2.,"Likelihood_weight_gen","Events");
  mtautau_gen=HConfig.GetTH1D(Name+"-ditau_gen", "ditau_gen", 150.,0.,150.,"M_tautau_Final_gen","Events");


  mtautau_all_1_gen=HConfig.GetTH1D(Name+"_ditau_all_averagebybin_gen", "ditau_all_averagebybin_gen", 100.,0.,200.,"M_ditau_all_averagebybin_gen","Events"); 
  mtautau_all_2_gen=HConfig.GetTH1D(Name+"_ditau_all_averagebyall_gen", "ditau_all_averagebyall_gen", 100.,0.,200.,"M_ditau_all_averagebyall_gen","Events");
  mtautau_bm_1_gen=HConfig.GetTH1D(Name+"_ditau_bm_averagebybin_gen", "ditau_bm_averagebybin_gen", 100.,0.,200.,"M_ditau_bm_averagebybin_gen","Events");  
  mtautau_bm_2_gen=HConfig.GetTH1D(Name+"_ditau_bm_averagebyall_gen", "ditau_bm_averagebyall_gen", 100.,0.,200.,"M_ditau_bm_averagebyall_gen","Events");  
  mtautau_bw_1_gen=HConfig.GetTH1D(Name+"_ditau_bw_averagebybin_gen", "ditau_bw_averagebybin_gen", 100.,0.,200.,"M_ditau_bw_averagebybin_gen","Events");  
  mtautau_bw_2_gen=HConfig.GetTH1D(Name+"_ditau_bw_averagebyall_gen", "ditau_bw_averagebyall_gen", 100.,0.,200.,"M_ditau_bw_averagebyall_gen","Events");  

  mtautau_all_1=HConfig.GetTH1D(Name+"_ditau_all_averagebybin", "ditau_all_averagebybin", 100.,0.,200.,"M_ditau_all_averagebybin","Events"); 
  mtautau_all_2=HConfig.GetTH1D(Name+"_ditau_all_averagebyall", "ditau_all_averagebyall", 100.,0.,200.,"M_ditau_all_averagebyall","Events");
  mtautau_bm_1=HConfig.GetTH1D(Name+"_ditau_bm_averagebybin", "ditau_bm_averagebybin", 100.,0.,200.,"M_ditau_bm_averagebybin","Events");  
  mtautau_bm_2=HConfig.GetTH1D(Name+"_ditau_bm_averagebyall", "ditau_bm_averagebyall", 100.,0.,200.,"M_ditau_bm_averagebyall","Events");  
  mtautau_bw_1=HConfig.GetTH1D(Name+"_ditau_bw_averagebybin", "ditau_bw_averagebybin", 100.,0.,200.,"M_ditau_bw_averagebybin","Events");  
  mtautau_bw_2=HConfig.GetTH1D(Name+"_ditau_bw_averagebyall", "ditau_bw_averagebyall", 100.,0.,200.,"M_ditau_bw_averagebyall","Events");  

 Taupt_all_1=HConfig.GetTH1D(Name+"_taupt_all_averagebybin", "taupt_all_averagebybin", 100.,0.,200.,"taupt_all_averagebybin","Events");
 TauE_all_1=HConfig.GetTH1D(Name+"_tauE_all_averagebybin", "tauE_all_averagebybin", 200.,0.,400.,"tauE_all_averagebybin","Events");
 H_SPIN_OBS_X1=HConfig.GetTH1D(Name+"_SPIN_OBS_X1", "SPIN_OBS_X1", 60.,0.,1.,"SPIN_OBS_X1","Events");
 H_SPIN_OBS_X2=HConfig.GetTH1D(Name+"_SPIN_OBS_X2", "SPIN_OBS_X2", 60.,0.,1.,"SPIN_OBS_X2","Events");
 H_SPIN_OBS_X3=HConfig.GetTH1D(Name+"_SPIN_OBS_X3", "SPIN_OBS_X3", 60.,0.,1.,"SPIN_OBS_X3","Events");
 H_SPIN_OBS_X4=HConfig.GetTH1D(Name+"_SPIN_OBS_X4", "SPIN_OBS_X4", 60.,0.,1.,"SPIN_OBS_X4","Events");
H_SPIN_OBS_X1_1=HConfig.GetTH1D(Name+"_SPIN_OBS_X1_1", "SPIN_OBS_X1_1", 60.,0.,1.,"SPIN_OBS_X1_1","Events");
H_SPIN_OBS_X1=HConfig.GetTH1D(Name+"_SPIN_OBS_X1", "SPIN_OBS_X1", 60.,0.,1.,"SPIN_OBS_X1","Events");
H_SPIN_OBS_X1_gen=HConfig.GetTH1D(Name+"_SPIN_OBS_X1_gen", "SPIN_OBS_X1_gen", 60.,0.,1.,"SPIN_OBS_X1_gen","Events");
H_SPIN_OBS_X2_gen=HConfig.GetTH1D(Name+"_SPIN_OBS_X2_gen", "SPIN_OBS_X2_gen", 60.,0.,1.,"SPIN_OBS_X2_gen","Events");
H_SPIN_OBS_X3_gen=HConfig.GetTH1D(Name+"_SPIN_OBS_X3_gen", "SPIN_OBS_X3_gen", 60.,0.,1.,"SPIN_OBS_X3_gen","Events");


  mtautau_all_Res_1=HConfig.GetTH1D(Name+"_ditau_all_averagebybin_Res", "ditau_all_averagebybin_Res", 100.,-2,2.,"M_ditau_all_averagebybin_Res","Events");
  mtautau_all_Res_2=HConfig.GetTH1D(Name+"_ditau_all_averagebyall_Res", "ditau_all_averagebyall_Res", 100.,-2,2.,"M_ditau_all_averagebyall_Res","Events");
  mtautau_all_Resp_Matr_1=HConfig.GetTH2D(Name+"_ditau_all_averagebybin_Resp_Matr", "ditau_all_averagebybin_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_all_averagebybin_Resp_Matr","Events");
  mtautau_all_Resp_Matr_2=HConfig.GetTH2D(Name+"_ditau_all_averagebyall_Resp_Matr", "ditau_all_averagebyall_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_all_averagebyall_Resp_Matr","Events");

  mtautau_bm_Res_1=HConfig.GetTH1D(Name+"_ditau_bm_averagebybin_Res", "ditau_bm_averagebybin_Res", 100.,-2,2.,"M_ditau_bm_averagebybin_Res","Events");
  mtautau_bm_Res_2=HConfig.GetTH1D(Name+"_ditau_bm_averagebyall_Res", "ditau_bm_averagebyall_Res", 100.,-2,2.,"M_ditau_bm_averagebyall_Res","Events");

  mtautau_bm_Resp_Matr_1=HConfig.GetTH2D(Name+"_ditau_bm_averagebybin_Resp_Matr", "ditau_bm_averagebybin_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_bm_averagebybin_Resp_Matr","Events");
  mtautau_bm_Resp_Matr_2=HConfig.GetTH2D(Name+"_ditau_bm_averagebyall_Resp_Matr", "ditau_bm_averagebyall_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_bm_averagebyall_Resp_Matr","Events");

  mtautau_bw_Res_1=HConfig.GetTH1D(Name+"_ditau_bw_averagebybin_Res", "ditau_bw_averagebybin_Res", 100.,-2,2.,"M_ditau_bw_averagebybin_Res","Events");
  mtautau_bw_Res_2=HConfig.GetTH1D(Name+"_ditau_bw_averagebyall_Res", "ditau_bw_averagebyall_Res", 100.,-2,2.,"M_ditau_bw_averagebyall_Res","Events");
  mtautau_bw_Resp_Matr_1=HConfig.GetTH2D(Name+"_ditau_bw_averagebybin_Resp_Matr", "ditau_bw_averagebybin_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_bw_averagebybin_Resp_Matr","Events");
  mtautau_bw_Resp_Matr_2=HConfig.GetTH2D(Name+"_ditau_bw_averagebyall_Resp_Matr", "ditau_bw_averagebyall_Resp_Matr", 100.,0,200., 100.,0,200.,"M_ditau_bw_averagebyall_Resp_Matr","Events");


  // Di-Tau mass
  mtautau11=HConfig.GetTH1D(Name+"_ditau11", "ditau11", 150.,0.,150.,"M_tautau_11","Events");
  mtautau12=HConfig.GetTH1D(Name+"_ditau12", "ditau12", 150.,0.,150.,"M_tautau_12","Events");
  mtautau21=HConfig.GetTH1D(Name+"_ditau21", "ditau21", 150.,0.,150.,"M_tautau_21","Events");
  mtautau22=HConfig.GetTH1D(Name+"_ditau22", "ditau22", 150.,0.,150.,"M_tautau_22","Events");
  mtautau3=HConfig.GetTH1D(Name+"_ditau3", "ditau3", 150.,0.,150.,"M_tautau_3","Events");
  mtautau4=HConfig.GetTH1D(Name+"_ditau4", "ditau4", 150.,0.,150.,"M_tautau_4","Events");
  mtautau5=HConfig.GetTH1D(Name+"-ditau5", "ditau5", 150.,0.,150.,"M_tautau_5","Events");
  Diff_Spin_obs_X1=HConfig.GetTH1D(Name+"-diff_Spin_Obs_X1", "diff_Spin_Obs_X1", 60.,-1,1.,"diff_Spin_Obs_X1","Events");
  Diff_Spin_obs_X2=HConfig.GetTH1D(Name+"-diff_Spin_Obs_X2", "diff_Spin_Obs_X2", 60.,-1,1.,"diff_Spin_Obs_X2","Events");
  Diff_Spin_obs_X3=HConfig.GetTH1D(Name+"-diff_Spin_Obs_X3", "diff_Spin_Obs_X3", 60.,-1,1.,"diff_Spin_Obs_X3","Events");
  

Lw=HConfig.GetTH1D(Name+"-Lw", "Lw", 150.,-5.,2.,"Likelihood_weight","Events");
  mtautau=HConfig.GetTH1D(Name+"-ditau", "ditau", 150.,0.,150.,"M_tautau_Final","Events");

  Mzgenmmc_gen=HConfig.GetTH2D(Name+"-Mzgenmmc_gen", "Mzgenmmc_gen", 150.,0.,150.,150.,0.,150.);
  Mzrecmmc_gen=HConfig.GetTH2D(Name+"-Mzrecmmc_gen", "Mzrecmmc_gen", 150.,0.,150.,150.,0.,150.);
  Mzrecmmc_genmmc=HConfig.GetTH2D(Name+"-Mzrecmmc_genmmc", "Mzrecmmc_genmmc", 150.,0.,150.,150.,0.,150.);

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  MMC_PolarizedTau::Store_ExtraDist() {
  ////////////////////////////////////// 

  Extradist1d.push_back(&spin_WT);
  //Gen objects
  Extradist1d.push_back(&M_TauPi_gen);        Extradist1d.push_back(&Pt_TauPi_gen);
  Extradist1d.push_back(&Eta_TauPi_gen);      Extradist1d.push_back(&Phi_TauPi_gen);
  Extradist1d.push_back(&M_TauMu_gen);        Extradist1d.push_back(&Pt_TauMu_gen);
  Extradist1d.push_back(&Eta_TauMu_gen);      Extradist1d.push_back(&Phi_TauMu_gen);        
  Extradist1d.push_back(&Mz_gen);             Extradist1d.push_back(&Mz_gen1);
  //Reco objects
  Extradist1d.push_back(&M_TauPi_rec);        Extradist1d.push_back(&Pt_TauPi_rec);
  Extradist1d.push_back(&Eta_TauPi_rec);      Extradist1d.push_back(&Phi_TauPi_rec);
  Extradist1d.push_back(&M_TauMu_rec);        Extradist1d.push_back(&Pt_TauMu_rec);
  Extradist1d.push_back(&Eta_TauMu_rec);      Extradist1d.push_back(&Phi_TauMu_rec);        
  // Delta_R (,)
  Extradist1d.push_back(&dR_recMu_genMu);     Extradist1d.push_back(&dR_recPi_genPi);
  Extradist1d.push_back(&mindR_recMu_genMu);  Extradist1d.push_back(&mindR_recPi_genPi);
  // MET infos
  Extradist1d.push_back(&MetX);               Extradist1d.push_back(&MetY);
  Extradist1d.push_back(&MetX_gen);           Extradist1d.push_back(&MetY_gen);

  // Di-Tau mass
  Extradist1d.push_back(&mtautau11_gen);      Extradist1d.push_back(&mtautau12_gen);
  Extradist1d.push_back(&mtautau21_gen);      Extradist1d.push_back(&mtautau22_gen);
  Extradist1d.push_back(&mtautau3_gen);       Extradist1d.push_back(&mtautau4_gen); 
  Extradist1d.push_back(&mtautau5_gen);       Extradist1d.push_back(&mtautau_gen); 
  Extradist1d.push_back(&Lw_gen); 

  Extradist1d.push_back(&mtautau_all_1_gen);  Extradist1d.push_back(&mtautau_all_2_gen);
  Extradist1d.push_back(&mtautau_bm_1_gen);   Extradist1d.push_back(&mtautau_bm_2_gen);
  Extradist1d.push_back(&mtautau_bw_1_gen);   Extradist1d.push_back(&mtautau_bw_2_gen);
//Spin Observables
 Extradist1d.push_back(&H_SPIN_OBS_X1);
 Extradist1d.push_back(&H_SPIN_OBS_X2);
 Extradist1d.push_back(&H_SPIN_OBS_X3);
 Extradist1d.push_back(&H_SPIN_OBS_X4);
Extradist1d.push_back(&H_SPIN_OBS_X1_1);
 Extradist1d.push_back(&H_SPIN_OBS_X1_gen);
Extradist1d.push_back(&H_SPIN_OBS_X2_gen);
Extradist1d.push_back(&H_SPIN_OBS_X3_gen);


 Extradist1d.push_back(&Diff_Spin_obs_X1);
 Extradist1d.push_back(&Diff_Spin_obs_X2);
 Extradist1d.push_back(&Diff_Spin_obs_X3);


  // Di-Tau mass
  Extradist1d.push_back(&mtautau11);          Extradist1d.push_back(&mtautau12);
  Extradist1d.push_back(&mtautau21);          Extradist1d.push_back(&mtautau22);
  Extradist1d.push_back(&mtautau3);           Extradist1d.push_back(&mtautau4); 
  Extradist1d.push_back(&mtautau5);           Extradist1d.push_back(&mtautau); 
  Extradist1d.push_back(&Lw); 

  Extradist1d.push_back(&mtautau_all_1);      Extradist1d.push_back(&mtautau_all_2);
  Extradist1d.push_back(&mtautau_bm_1);       Extradist1d.push_back(&mtautau_bm_2);
  Extradist1d.push_back(&mtautau_bw_1);       Extradist1d.push_back(&mtautau_bw_2);

  Extradist1d.push_back(&mtautau_all_Res_1);Extradist1d.push_back(&mtautau_all_Res_2);
  Extradist1d.push_back(&mtautau_bm_Res_1);Extradist1d.push_back(&mtautau_bm_Res_2);
  Extradist1d.push_back(&mtautau_bw_Res_1);Extradist1d.push_back(&mtautau_bw_Res_2);

  Extradist2d.push_back(&mtautau_all_Resp_Matr_1);Extradist2d.push_back(&mtautau_all_Resp_Matr_2);
  Extradist2d.push_back(&mtautau_bm_Resp_Matr_1);Extradist2d.push_back(&mtautau_bm_Resp_Matr_2);
  Extradist2d.push_back(&mtautau_bw_Resp_Matr_1);Extradist2d.push_back(&mtautau_bw_Resp_Matr_2);


  
  Extradist2d.push_back(&Mzgenmmc_gen);       Extradist2d.push_back(&Mzrecmmc_gen);
  Extradist2d.push_back(&Mzrecmmc_genmmc);
   
}

/*

double f_const1 (double *x,double *par){
  if (x[0]<42.5) return TMath::Exp(par[0]+x[0]*par[1]);
  else if( x[0]>42.5 && x[0]<52.6  ) return par[2]+par[3]*x[0];
  else return TMath::Exp(par[4]+x[0]*par[5]) + par[6]+par[7]*x[0];
}

double f_const2 (double *x,double *par){
  if (x[0]<42.5) return TMath::Exp(par[0]+x[0]*par[1]);
  else if( x[0]>42.5 && x[0]<52.6  ) return TMath::Exp(par[2]+par[3]*x[0]);
  else return TMath::Exp(par[4]+x[0]*par[5]) + par[6]+par[7]*x[0];
}

double f_const2_had (double *x,double *par){
  if (x[0]<42.5) return TMath::Exp(par[0]+x[0]*par[1])+par[2]+x[0]*par[3];
  else if( x[0]>42.5 && x[0]<52.6  ) return TMath::Exp(par[4]+par[5]*x[0])+ par[6]+x[0]*par[7];
  else return TMath::Exp(par[8]+x[0]*par[9]) + par[10]+par[11]*x[0];
}

TF1  MMC_PolarizedTau::DeltaR_lep(double pt) {
  TF1 const1 ("const1",f_const1,0,100,8);
  const1.SetParameter(0,7.014037);    const1.SetParameter(1,0.08982361);
  const1.SetParameter(2,260756.1);    const1.SetParameter(3,-4933.671);
  const1.SetParameter(4,21.574);      const1.SetParameter(5,-0.2743597);
  const1.SetParameter(6,874.4552);    const1.SetParameter(7,-7.994477);
  
  TF1 mean("mean","expo(0)+pol1(2)",10.72,85.24);
  mean.SetParameter(0,-1.319604);     mean.SetParameter(1,-0.0698018);
  mean.SetParameter(2,0.05926357);    mean.SetParameter(3,-0.0004089469);
     
  TF1 sigma1("sigma1","expo+pol1(2)",10.72,109);
  sigma1.SetParameter(0,-2.227225);     sigma1.SetParameter(1,-0.04167413);
  sigma1.SetParameter(2,6.679525e-005); sigma1.SetParameter(3,0.0001051946);
    
  TF1 const2("const2",f_const2,0,100,8);
  const2.SetParameter(0,8.641931);   const2.SetParameter(1,0.09052246);
  const2.SetParameter(2,27.88014);   const2.SetParameter(3,-0.3550818);
  const2.SetParameter(4,19.03704);   const2.SetParameter(5,-0.1925955);  
  const2.SetParameter(6,2661.081);   const2.SetParameter(7,-23.63197);
   
  TF1 MPV("MPV","expo+pol1(2)",10.72,96.04);
  MPV.SetParameter(0,-0.8407024);    MPV.SetParameter(1,-0.06564579);
  MPV.SetParameter(2,0.07128014);    MPV.SetParameter(3,-0.0004138105);
   
  TF1 sigma2("sigma2","expo+pol1(2)",11.8,92.8);
  sigma2.SetParameter(0,-2.364371);  sigma2.SetParameter(1,-0.09803685);
  sigma2.SetParameter(2,0.01046975); sigma2.SetParameter(3,-8.072633e-005);
    
  TF1 total ("DeltaR","gaus(0)+landau(3)",0,1);
  double par[6];
  par[0]=const1.Eval(pt);  par[1]=mean.Eval(pt);
  par[2]=sigma1.Eval(pt);  par[3]=const2.Eval(pt);
  par[4]=MPV.Eval(pt);     par[5]=sigma2.Eval(pt);
  total.SetParameters(par);       
  return total; }

TF1  MMC_PolarizedTau::DeltaR_had(double pt) {
  TF1 const2("const2",f_const2_had,0,100,12);
  const2.SetParameter(0,9.048143);    const2.SetParameter(1,0.1029005);
  const2.SetParameter(2,-93188.71);   const2.SetParameter(3,7811.768);
  const2.SetParameter(4,15.14293);    const2.SetParameter(5,-0.02232439);
  const2.SetParameter(6,1898838);     const2.SetParameter(7,-57810.23);
  const2.SetParameter(8,21.52146);    const2.SetParameter(9,-0.212765);
  const2.SetParameter(10,9748.511);   const2.SetParameter(11,-93.85233);
  
  TF1 MPV("MPV","expo+pol1(2)",9.64,98.2);
  MPV.SetParameter(0,-0.7645296);    MPV.SetParameter(1,-0.06920899);
  MPV.SetParameter(2,0.08500593);     MPV.SetParameter(3,-0.0005090078);
  
  TF1 sigma2("sigma2","expo+pol1(2)",1,90.64);
  sigma2.SetParameter(0,-3.007221);    sigma2.SetParameter(1,-0.08821775);
  sigma2.SetParameter(2,0.004511857);  sigma2.SetParameter(3,-2.062199e-05);
  
  TF1 total("DeltaR","gaus(0)+landau(3)",0,0.4);
  double par[6];
  par[0]=0;    par[1]=0;   par[2]=0;
  par[3]=const2.Eval(pt);  par[4]=MPV.Eval(pt);
  par[5]=sigma2.Eval(pt);  total.SetParameters(par);
  return total; }

*/



double MMC_PolarizedTau::DeltaR_had(double pt,double deltaR){
double const2;
if (pt<42.5) const2= TMath::Exp(9.048143+pt*0.1029005)-93188.71+7811.768*pt;

else if( pt>42.5 && pt<52.6  ) const2=TMath::Exp(15.14293-0.02232439*pt)+1898838-57810.23*pt;


else const2= TMath::Exp(21.52146+pt*-0.212765)+9748.511-93.85233*pt;

double MPV= TMath::Exp(-0.7645296-0.06920899*pt)+0.08500593-0.0005090078*pt;

double sigma2=TMath::Exp(-3.007221-0.08821775*pt)+0.004511857-2.062199e-05*pt;

double total=const2*TMath::Landau(deltaR,MPV,sigma2);
 
return  total<0 ? 0:total;
	   
}
	


 
 

double MMC_PolarizedTau::DeltaR_lep(double pt, double deltaR){
   
  double const1,const2;
   if (pt<42.5) {
     const1= TMath::Exp(7.014037+pt*0.08982361); 
     const2=TMath::Exp(8.641931+pt*0.09052246);
     
     
  }
  else if( pt>42.5 && pt<52.6  ) {
    const1= 260756.1-4933.671*pt;
    const2=TMath::Exp(27.88014-0.3550818*pt);
    
    
  }
  else {
    const1= TMath::Exp(21.574-0.2743597*pt) + 874.4552 -7.994477*pt; 
    const2=TMath::Exp(19.03704-0.1925955*pt) +2661.081-23.63197*pt;
  } 
   
 double mean=  TMath::Exp(-1.319604-0.0698018*pt)+0.05926357-0.0004089469*pt;
  double sigma1=TMath::Exp(-2.227225-0.04167413*pt)+6.679525e-005+0.0001051946*pt;
  double MPV=TMath::Exp(-0.8407024-0.06564579*pt)+0.07128014-0.0004138105*pt;
  double sigma2=TMath::Exp(-2.364371-0.09803685*pt)+0.01046975-8.072633e-005*pt;
    
   double  total=const1*TMath::Gaus(deltaR,mean,sigma1)+const2*TMath::Landau(deltaR,MPV,sigma2);
   
	
  return   total<0?0:total;
  
   
   
 
   
}
void  MMC_PolarizedTau::doEvent() {
  unsigned int t(0);
  int id=Ntp->GetMCID();
//  std::cout<<"ID ............................ "<<id<<std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)) {t=2;}// std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,Boson3_idx,Boson4_idx,tau1_idx(0),tau2_idx(0),tau3_idx(0),tau4_idx(0);
  value.at(isZtautautopimu)=0;
  pass.at(isZtautautopimu) = true;
  if(pass.at(isZtautautopimu))value.at(isZtautautopimu)=1;
 tt=t; 
  double wobs=1;
  double w=1;
  //  if(verbose)  std::cout << "MCID=="<< Ntp->GetMCID() << "||Size==" << Npassed.size() << "||t== " << t << "||BosonIdx== " << Boson_idx << "||Tau1== " << tau1_idx << "||Tau2== " << tau1_idx << "||Ntaus== " << Ntp->NMCTaus() << std::endl;
  
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status) {
   // if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    // reject unsupported modes
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson_idx,tau1_idx,tau2_idx)) {
      bool ImpTau=true;
      double jakid1=Ntp->MCTau_JAK(tau1_idx);
      double jakid2=Ntp->MCTau_JAK(tau2_idx);
      if(!(jakid1==TauDecay::JAK_PION || jakid1==TauDecay::JAK_MUON)) ImpTau=false;
      if(!(jakid2==TauDecay::JAK_PION || jakid2==TauDecay::JAK_MUON)) ImpTau=false;
      //if(jakid1!=jakid2) return;
      if(!ImpTau) { std::cout << "Decay modes not implemented in TauSpinner " << std::endl; return;}
    }
    else {
      std::cout << "Not a valid Z0->tau+tau- decay" <<std::endl;
      return;
    }
    //double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double Spin_WT=1;
   // std::cout << "Spin_WT ==" << Spin_WT << std::endl;
    //double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    //double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    double hplus=1;
    double hminus=1;
 //   std::cout << "hplus " << hplus << " hminus " << hminus << std::endl;
//std::cout<<"begin"<<std::endl;
//for (unsigned int i=0;i<Ntp->NPFTaus();i++)  std::cout<<Ntp->PFTau_hpsDecayMode(i)<<" "<<Ntp->PFTau_isHPSByDecayModeFinding(i)<<std::endl;
//std::cout<<"end"<<std::endl;


std::vector<int> tau_idx;
    for(unsigned int i=0;i<Ntp->NPFTaus();i++){
	    if(Ntp->PFTau_isTightIsolationDBSumPtCorr(i) && Ntp->PFTau_isHPSAgainstElectronsTight(i) && Ntp->PFTau_isHPSAgainstMuonTight(i) && Ntp->PFTau_p4(i).Pt()>20){
		    tau_idx.push_back(i);
	    }
    }
// for (unsigned int i=0;i<tau_idx.size();i++)  std::cout<<Ntp->PFTau_hpsDecayMode(tau_idx.at(i))<<std::endl;
    spin_WT.at(t).Fill(Spin_WT,w);
  
    if(verbose)std::cout<< "A " <<std::endl;
    /////////////////////////////////    
    //  Muon && Pion  Channel
    /////////////////////////////////
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson3_idx,TauDecay::JAK_MUON,tau3_idx) && Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson4_idx,TauDecay::JAK_PION,tau4_idx)) { //JAK_MUON && JAK_PION
  //int ddd=Ntp->GetMCID();
  // cout<<"MC type "<<ddd<<std::endl; 
 if(verbose)std::cout<< "pi-mu " <<std::endl;
      
      TLorentzVector GenZ=Ntp->MCSignalParticle_p4(Boson3_idx);
      TLorentzVector GenTauMu(0,0,0,0);       TLorentzVector GenMu(0,0,0,0);  
      TLorentzVector GenNum(0,0,0,0);         TLorentzVector GenNutm(0,0,0,0);
      TLorentzVector Muon(0,0,0,0);           TLorentzVector TauMuRec(0,0,0,0);
      TLorentzVector MuonRec(0,0,0,0);        TLorentzVector Muon_mdR(0,0,0,0);
      TLorentzVector GenNumtm(0,0,0,0);
      
      TLorentzVector GenZ1=Ntp->MCSignalParticle_p4(Boson4_idx);
      TLorentzVector GenTauPi(0,0,0,0);      TLorentzVector GenPi(0,0,0,0);
      TLorentzVector GenNutp(0,0,0,0);       TLorentzVector Pion(0,0,0,0);
      TLorentzVector Zrec(0,0,0,0);          TLorentzVector TauPiRec(0,0,0,0);
      TLorentzVector PionRec(0,0,0,0);       TLorentzVector Pion_mdR(0,0,0,0);  
      TLorentzVector Z_gen(0,0,0,0);         TLorentzVector Z_gen1(0,0,0,0);  
      double mzgen= -100;                    double mzgen1= -100;
      
      double dR_Muons = -1;   double dR_Pions = -1;   double mindRMu = 999; double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index = -1; int mindRPi_Index = -1; //double Zrec_Mt = -1;  double Zgen_Mt = -1;
      
      double m_tau= 1.777;     double mass_Z=91.18;
      TLorentzVector p4miss11_gen(0,0,0,0),p4miss12_gen(0,0,0,0),p4miss21_gen(0,0,0,0),p4miss22_gen(0,0,0,0),p4vis1_gen(0,0,0,0),p4vis2_gen(0,0,0,0);
      TLorentzVector p41_gen(0,0,0,0),p42_gen(0,0,0,0), p43_gen(0,0,0,0),miss(0,0,0,0);
      double METX_gen= -1000;     double METY_gen= -1000;         double pt_miss1_gen= -1000;   double pt_miss2_gen= -1000; 
      double p_vis1_gen= -1000;   double theta_vis1_gen= -1000;   double phi_vis1_gen= -1000;   double m_vis1_gen= -1000; 
      double p_vis2_gen= -1000;   double theta_vis2_gen= -1000;   double phi_vis2_gen= -1000;   double m_vis2_gen= -1000; 
      // double p_miss1_gen= -1000;  double theta_miss1_gen= -1000;  double phi_miss1_gen= -1000;  double m_miss1_gen= 0;
      // double p_miss2_gen= -1000;  double theta_miss2_gen= -1000;  double phi_miss2_gen= -1000;  double m_miss2_gen= -1000; 
      double pz11_gen= -1000;     double pz12_gen= -1000;         double pz21_gen= -1000;       double pz22_gen= -1000;
      double  w11_gen=1.; double w12_gen=1.; double w21_gen=1.; double w22_gen=1.;//weight
      
      TLorentzVector p4miss11(0,0,0,0),p4miss12(0,0,0,0),p4miss21(0,0,0,0),p4miss22(0,0,0,0),p4vis1(0,0,0,0),p4vis2(0,0,0,0);
      TLorentzVector p41(0,0,0,0),p42(0,0,0,0), p43(0,0,0,0);
      double METX= -1000;         double METY= -1000;   double pt_miss1= -1000;       double pt_miss2= -1000;  
      double p_vis1= -1000;       double theta_vis1= -1000;       double phi_vis1= -1000;       double m_vis1= -1000; 
      double p_vis2= -1000;       double theta_vis2= -1000;       double phi_vis2= -1000;       double m_vis2= -1000; 
      //double p_miss1= -1000;      double theta_miss1= -1000;      double p_miss2= -1000;      double theta_miss2= -1000;  
      double phi_miss1= -1000;    double phi_miss2= -1000;        double m_miss1= 0;            double m_miss2= -1000; 
      double pz11= -1000;         double pz12= -1000;             double pz21= -1000;           double pz22= -1000;
     // double  w11 =1.; double w12=1.; double w21=1.; double w22= 1.; 
     // double w11_R,w12_R,w21_R,w22_R; 
      
      

      ////**** Gen- Reco Infos ****////        
      //Gen Muon//
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau3_idx);i++) {
        if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauMu=Ntp->MCTauandProd_p4(tau3_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::mu_plus)) {
          GenMu+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_mu)) {
          GenNum+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(tau3_idx,i))==abs(PDGInfo::nu_tau)) {
          GenNutm+=Ntp->MCTauandProd_p4(tau3_idx,i);
	}
      } 
      //Reco Muon
      for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu);
	dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	if(dR_Muons < mindRMu) { mindRMu = dR_Muons;  mindRMu_Index = i; }
	MuonRec = Ntp->Muon_p4(mindRMu_Index);
      }
      if(mindRMu<0.1 && mindRMu_Index>-1) { 
	mindR_recMu_genMu.at(t).Fill(mindRMu,w);   
	//	std::cout<< "Moun_minDR==" << mindRMu <<std::endl;
	Muon_mdR=MuonRec;
      }
      //Gen Pion//
      for(int i=0; i<Ntp->NMCTauDecayProducts(tau4_idx);i++) { 
        if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::tau_minus)) {
          GenTauPi=Ntp->MCTauandProd_p4(tau4_idx,i);
        }
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::pi_plus) ) {
          GenPi+=Ntp->MCTauandProd_p4(tau4_idx,i);
	}
        else if(abs(Ntp->MCTauandProd_pdgid(tau4_idx,i))==abs(PDGInfo::nu_tau) ) {
	  GenNutp+=Ntp->MCTauandProd_p4(tau4_idx,i);
	}
      }      
      //Reco Pion
      for(unsigned int i=0; i < Ntp->NPFTaus(); i++) {
	dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi);
	dR_recPi_genPi.at(t).Fill(dR_Pions,w);
	if(dR_Pions < mindRPi) { mindRPi = dR_Pions;  mindRPi_Index = i; }
	PionRec = Ntp->PFTau_p4(mindRPi_Index);
      }
      if(mindRPi<0.4 && mindRPi_Index>-1) { 
	mindR_recPi_genPi.at(t).Fill(mindRPi,w);
	Pion_mdR=PionRec;
      }     

      //Gen & Reco MET            
      GenNumtm = GenNum + GenNutm;       // sum of neutrinos from leptonic side.
      miss     = GenNutp + GenNum + GenNutm; // sum of neutrinos from both sides lep & had.
      Z_gen    = GenTauMu + GenTauPi;
      Z_gen1   = GenMu + GenPi + miss;
    std::cout<<"Miss M "<<miss.M()<<std::endl;

      /////////////////////////
      ////////   MET   ////////
      /////////////////////////
      //Gen
      METX_gen = miss.E()*sin(miss.Theta())*cos(miss.Phi());
      METY_gen = miss.E()*sin(miss.Theta())*sin(miss.Phi());
      
      // //Reco
      //METX = Ntp->MET_Uncorr_ex(); //
     // METX = Ntp->MET_CorrMVA_ex();     
      //METY = Ntp->MET_Uncorr_ey(); //
     //METY =  Ntp->MET_CorrMVA_ey();
      METX = Ntp->MET_CorrMVAMuTau_ex();
      METY = Ntp->MET_CorrMVAMuTau_ey();
    double MET=sqrt(METX*METX+METY*METY);
     

     TVector2 MET_(METX,METY);
     


TVector2 MET_sig_(sqrt(Ntp->MET_CorrMVA_significance_xx()),sqrt(Ntp->MET_CorrMVA_significance_yy()));
std::cout<<MET_sig_.X()<<" 1 "<<MET_sig_.Y()<<std::endl;
TVector2 MET_sig_1(sqrt(Ntp->MET_CorrMVAMuTau_significance_xx()),sqrt(Ntp->MET_CorrMVAMuTau_significance_yy()));
std::cout<<MET_sig_1.X()<<" 2 "<<MET_sig_1.Y()<<std::endl;

     if(MET<10) return;  
      ////////////////////////////
      ////**** MMC infos **** ////
      ////////////////////////////
      //Gen
      p_vis1_gen= GenPi.P(); theta_vis1_gen= GenPi.Theta(); phi_vis1_gen= GenPi.Phi(); m_vis1_gen= GenPi.M(); 
      p_vis2_gen= GenMu.P(); theta_vis2_gen= GenMu.Theta(); phi_vis2_gen= GenMu.Phi(); m_vis2_gen= GenMu.M();  
      // Reco
      p_vis1= Pion_mdR.P(); theta_vis1= Pion_mdR.Theta(); phi_vis1= Pion_mdR.Phi(); m_vis1= Pion_mdR.M(); 
      p_vis2= Muon_mdR.P(); theta_vis2= Muon_mdR.Theta(); phi_vis2= Muon_mdR.Phi(); m_vis2= Muon_mdR.M();     
      
     
      ////**** end MMC infos **** ////
      
      // make sure that Spin_WT !=NaN
      if(Spin_WT > 0. &&  Spin_WT < 2.) {
	spin_WT.at(t).Fill(Spin_WT,w);

	//// MET-Plots ////
	/////////////////////////////////
	MetX.at(t).Fill(METX,w);	 MetY.at(t).Fill(METY,w);
	MetX_gen.at(t).Fill(METX_gen,w); MetY_gen.at(t).Fill(METY_gen,w);

	/////////////////////////////////
	////     Z_gen  infos      ////
	/////////////////////////////////                                                                                                                                                   
	mzgen= Z_gen.M();
	mzgen1= Z_gen1.M();

        if (mzgen<50) return;
	if(Z_gen.E()>0) { Mz_gen.at(t).Fill(mzgen,w); }
	if(Z_gen1.E()>0) { Mz_gen1.at(t).Fill(mzgen1,w); }

	/////////////////////////////////
	////     Tau_gen  infos      ////
	/////////////////////////////////                                                                                                                                                   
	if(GenTauMu.E()>0) {
	  M_TauMu_gen.at(t).Fill(GenTauMu.M(),w);	  Pt_TauMu_gen.at(t).Fill(GenTauMu.Pt(),w);
	  Eta_TauMu_gen.at(t).Fill(GenTauMu.Eta(),w);	  Phi_TauMu_gen.at(t).Fill(GenTauMu.Phi(),w); }
	if(GenTauPi.E()>0) {
	  M_TauPi_gen.at(t).Fill(GenTauPi.M(),w);	  Pt_TauPi_gen.at(t).Fill(GenTauPi.Pt(),w);
	  Eta_TauPi_gen.at(t).Fill(GenTauPi.Eta(),w);	  Phi_TauPi_gen.at(t).Fill(GenTauPi.Phi(),w); }

      	/////////////////////////////////
	//// Tau_rec  reconstruction ////
	/////////////////////////////////
	if(mindRMu < 0.1 && mindRPi < 0.4 && mindRMu_Index >-1 && mindRPi_Index >-1 && MET>10 && Pion_mdR.Pt() >10. &&  Muon_mdR.Pt()>10. )
{

n_sel_ev++;

//for (unsigned int i=0;i<tau_idx.size();i++)  std::cout<<Ntp->PFTau_hpsDecayMode(tau_idx.at(i))<<std::endl;
	
//for (unsigned int i=0;i<Ntp->NPFTaus();i++)  std::cout<<Ntp->PFTau_hpsDecayMode(i)<<std::endl;

double SPIN_OBS_X1_gen=GenPi.Pt()/GenTauPi.Pt();
H_SPIN_OBS_X1_gen.at(t).Fill(SPIN_OBS_X1_gen);
double mass_t2=pow(GenTauPi.M(),2)+pow(GenTauMu.M(),2)+2*(GenTauPi.Et()*GenTauMu.Et()-GenTauPi.Px()*GenTauMu.Px()-GenTauPi.Py()*GenTauMu.Py());
double SPIN_OBS_X2_gen=2*GenPi.Pt()/sqrt(mass_t2);

H_SPIN_OBS_X1_gen.at(t).Fill(SPIN_OBS_X1_gen);
H_SPIN_OBS_X2_gen.at(t).Fill(SPIN_OBS_X2_gen);

TLorentzVector GenPi_boosted=GenPi;
GenPi_boosted.Boost(-Z_gen.BoostVector());
double SPIN_OBS_X3_gen=2*GenPi_boosted.E()/Z_gen.M();
H_SPIN_OBS_X3_gen.at(t).Fill(SPIN_OBS_X3_gen);

  TauMuRec = Muon_mdR + GenNum + GenNutm;
	  TauPiRec = Pion_mdR + GenNutp;
	  Zrec = TauMuRec + TauPiRec;
	  if(TauMuRec.E()>0) {
	    M_TauMu_rec.at(t).Fill(TauMuRec.M(),w);
	    Pt_TauMu_rec.at(t).Fill(TauMuRec.Pt(),w);   
	    Eta_TauMu_rec.at(t).Fill(TauMuRec.Eta(),w);
	    Phi_TauMu_rec.at(t).Fill(TauMuRec.Phi(),w);
	  }
	  if(TauPiRec.E()>0) {	
	    M_TauPi_rec.at(t).Fill(TauPiRec.M(),w);
	    Pt_TauPi_rec.at(t).Fill(TauPiRec.Pt(),w);
	    Eta_TauPi_rec.at(t).Fill(TauPiRec.Eta(),w); 
	    Phi_TauPi_rec.at(t).Fill(TauPiRec.Phi(),w);
	  }
	
double Mass_Z=-999,Tau_pt=-999, SPIN_OBS_X1=-999,SPIN_OBS_X2=-999,SPIN_OBS_X3=-999,SPIN_OBS_X4=-999;
if (mmc_method(Pion_mdR,Muon_mdR,MET_,MET_sig_,Mass_Z,Tau_pt,SPIN_OBS_X1,SPIN_OBS_X2,SPIN_OBS_X3,SPIN_OBS_X4)) {

std::cout<<"MassZ "<<Mass_Z<<" Spin_OBs_X1_1 " << Pion_mdR.Pt()/Tau_pt<<" SPIN_OBS_X1 "<<SPIN_OBS_X1<<" SPIN_OBS_X2 "<<SPIN_OBS_X2<<" SPIN_OBS_X3 "<<SPIN_OBS_X3<<" SPIN_OBS_X4 "<<SPIN_OBS_X4<<std::endl;
Diff_Spin_obs_X1.at(t).Fill((SPIN_OBS_X1_gen-SPIN_OBS_X1)/SPIN_OBS_X1_gen);

Diff_Spin_obs_X2.at(t).Fill((SPIN_OBS_X2_gen-SPIN_OBS_X2)/SPIN_OBS_X2_gen);

Diff_Spin_obs_X3.at(t).Fill((SPIN_OBS_X3_gen-SPIN_OBS_X2)/SPIN_OBS_X3_gen);

H_SPIN_OBS_X1.at(t).Fill(SPIN_OBS_X1);
H_SPIN_OBS_X2.at(t).Fill(SPIN_OBS_X2);
H_SPIN_OBS_X3.at(t).Fill(SPIN_OBS_X3);
H_SPIN_OBS_X4.at(t).Fill(SPIN_OBS_X4);
H_SPIN_OBS_X1_1.at(t).Fill(Pion_mdR.Pt()/Tau_pt);



/*
if (Mass_Z>0) mtautau_all_1.at(t).Fill(Mass_Z);
if (Taupt>0) Taupt_all_1.at(t).Fill(Taupt);
if (TauE>0) TauE_all_1.at(t).Fill(TauE);
*/

}	   
	
std::cout<<"N. Sol events  gen and rec "<<n_sel_ev<<" "<<sol_gen_count<<" "<<mass_dum_count<<" "<< (float)mass_dum_count/n_sel_ev<<std::endl;
     
 
      } // if(midRmu........	
	
    } //end if( 0.< Spin_WT <2.)  
} //end JAK_MUON && JAK_PION   ////////////   end Muon && Pion    ////////////    
} //end if(Status)
} //end doEvent()

				
bool  MMC_PolarizedTau::mmc_method(const TLorentzVector &p4vis1,const TLorentzVector &p4vis2,const TVector2 MET,const TVector2 MET_sig,
double &mass_z,double &tau_pt, double &Spin_Obs_X1, double &Spin_Obs_X2,double &Spin_Obs_X3,double &Spin_Obs_X4)
{
     
  
		
TH1D mass_dummy_all("mass_dummy_all","mass_dummy_all",100,0,200);
TH1D hist_dummy_spin_Obs_X1("hist_dummy_spin_Obs_X1","hist_dummy_spin_Obs_X1",60,0,1);
TH1D hist_dummy_spin_Obs_X2("hist_dummy_spin_Obs_X2","hist_dummy_spin_Obs_X2",60,0,1);
TH1D hist_dummy_spin_Obs_X3("hist_dummy_spin_Obs_X3","hist_dummy_spin_Obs_X3",60,0,1);
TH1D hist_dummy_spin_Obs_X4("hist_dummy_spin_Obs_X4","hist_dummy_spin_Obs_X4",60,0,1);
TH1D taupt_dummy_all("taupt_dummy_all","taupt_dummy_all",100,0,200);	



double  phi_miss1, phi_miss2;
double p_vis1,theta_vis1, phi_vis1;
double m_vis1;
double m_tau=1.777;
double p_vis2,theta_vis2,phi_vis2;
double m_vis2;
double m_miss2,m_miss1=0;
TLorentzVector p4miss11,p4miss12,p4miss21,p4miss22,p4Z,p4tau_had,p4tau_lep;
double  w11,w12,w21,w22,w11_R,w12_R,w21_R,w22_R;
double term1,term2,term3;
double pz11,pz12,pz21,pz22;
double m11,m12,m21,m22;
double prob_metx,prob_mety;
double spin_Obs_X1,spin_Obs_X2,spin_Obs_X3,spin_Obs_X4;
bool check=false;
p_vis1=p4vis1.P();theta_vis1=p4vis1.Theta(); phi_vis1=p4vis1.Phi();
p_vis2=p4vis2.P();theta_vis2=p4vis2.Theta(); phi_vis2=p4vis2.Phi();
m_vis1=p4vis1.M();m_vis2=p4vis2.M();
double sig_metx=MET_sig.X();
double sig_mety=MET_sig.Y();

double METX=MET.X();
double METY=MET.Y();	
//double metx=METX;
//double mety=METY;

for (double  ii=phi_vis1-0.4;ii<=phi_vis1+0.4;ii+=0.05)  {
	phi_miss1=ii;					  
	for (double  jj=phi_vis2-0.4;jj<=phi_vis2+0.4;jj+=0.05){
	      phi_miss2=jj;
		for (double metx=METX-3*sig_metx;metx<=METX+3*sig_metx;metx+=2) 
		      for (double mety=METY-3*sig_mety;mety<=METY+3*sig_mety;mety+=2) {
//			        mety=METY;metx=METX;	
                                 if (pow((METX-metx)/sig_metx,2)+pow((METY-mety)/sig_mety,2)>pow(3,2)) continue;     
//                       	      double pt_miss1=(-sin(phi_miss2)*METX+cos(phi_miss2)*METY)/sin(phi_miss1-phi_miss2);
  //                              double pt_miss2=(sin(phi_miss1)*METX-cos(phi_miss1)*METY)/sin(phi_miss1-phi_miss2);

			       double pt_miss1=(-sin(phi_miss2)*metx+cos(phi_miss2)*mety)/sin(phi_miss1-phi_miss2);
				double pt_miss2=(sin(phi_miss1)*metx-cos(phi_miss1)*mety)/sin(phi_miss1-phi_miss2);
				if (pt_miss1<0 || pt_miss2<0) continue;		

				for (double  kk=0;kk<=1.5;kk+= 0.1)  {
					m_miss2=kk;
					
					if (m_miss2>m_tau-m_vis2) continue;
					
					double sqrt_term_1=sqrt((pow(m_vis1,2)+pow(p_vis1,2))*(pow(m_miss1,4)+pow(pow(m_vis1,2)-pow(m_tau,2),2)+4*cos(phi_miss1-phi_vis1)*sin(theta_vis1)*(-pow(m_vis1,2)+pow(m_tau,2))*p_vis1*pt_miss1-4*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(sin(phi_miss1-phi_vis1),2)*pow(p_vis1,2))*pow(pt_miss1,2)-2*pow(m_miss1,2)*(pow(m_vis1,2)+pow(m_tau,2)+2*sin(theta_vis1)*p_vis1*(sin(theta_vis1)*p_vis1+cos(phi_miss1-phi_vis1)*pt_miss1))));  
					
					if (TMath::IsNaN(sqrt_term_1)) continue;
					    
					double sqrt_term_2= sqrt((pow(m_vis2,2)+pow(p_vis2,2))*(pow(m_miss2,4)+pow(pow(m_vis2,2)-pow(m_tau,2),2)+4*cos(phi_miss2-phi_vis2)*sin(theta_vis2)*(-pow(m_vis2,2)+pow(m_tau,2))*p_vis2*pt_miss2-4*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(sin(phi_miss2-phi_vis2),2)*pow(p_vis2,2))*pow(pt_miss2,2)-2*pow(m_miss2,2)*(pow(m_vis2,2)+pow(m_tau,2)+2*sin(theta_vis2)*p_vis2*(sin(theta_vis2)*p_vis2+cos(phi_miss2-phi_vis2)*pt_miss2))));  
					
					if (TMath::IsNaN(sqrt_term_2))continue;
					 term1=(1./(2*(pow(m_vis1,2)+pow(sin(theta_vis1),2)*pow(p_vis1,2))));
					 term2=-cos(theta_vis1)*(pow(m_miss1,2)+pow(m_vis1,2)-pow(m_tau,2))*p_vis1;
					 term3=cos(phi_miss1-phi_vis1)*sin(2*theta_vis1)*pow(p_vis1,2)*pt_miss1;
					 pz11=term1*(term2+term3+sqrt_term_1);
					 pz12=-term1*(-term2-term3+sqrt_term_1);
					
					term1=(1./(2*(pow(m_vis2,2)+pow(sin(theta_vis2),2)*pow(p_vis2,2))));
					term2=-cos(theta_vis2)*(pow(m_miss2,2)+pow(m_vis2,2)-pow(m_tau,2))*p_vis2;
					term3=cos(phi_miss2-phi_vis2)*sin(2*theta_vis2)*pow(p_vis2,2)*pt_miss2;

					pz21=term1*(term2+term3+sqrt_term_2);
					pz22=-term1*(-term2-term3+sqrt_term_2);

//                                         std::cout<<pz21<<" "<<pz22<<" ";


//					std::cout<<pz21<<" "<<pz22<< " "<<std::endl;
					
                                        double tmp_X=pt_miss1*cos(phi_miss1),tmp_Y=pt_miss1*sin(phi_miss1);
                                        p4miss11.SetXYZT(tmp_X,tmp_Y,pz11,sqrt(pow(pz11,2)+pow(pt_miss1,2)+pow(m_miss1,2)));
					p4miss12.SetXYZT(tmp_X,tmp_Y,pz12,sqrt(pow(pz12,2)+pow(pt_miss1,2)+pow(m_miss1,2)));
					
					tmp_X=pt_miss2*cos(phi_miss2);tmp_Y=pt_miss2*sin(phi_miss2);	
					p4miss21.SetXYZT(tmp_X,tmp_Y,pz21,sqrt(pow(pz21,2)+pow(pt_miss2,2)+pow(m_miss2,2)));
					p4miss22.SetXYZT(tmp_X,tmp_Y,pz22,sqrt(pow(pz22,2)+pow(pt_miss2,2)+pow(m_miss2,2)));


                                        

					double taupt11_tmp=(p4vis1+p4miss11).Pt(),deltaR11_tmp=p4vis1.DeltaR(p4miss11);
					double taupt21_tmp=(p4vis2+p4miss21).Pt(),deltaR21_tmp=p4vis2.DeltaR(p4miss21);

					double taupt12_tmp=(p4vis1+p4miss12).Pt(),deltaR12_tmp=p4vis1.DeltaR(p4miss12);
					double taupt22_tmp=(p4vis2+p4miss22).Pt(),deltaR22_tmp=p4vis2.DeltaR(p4miss22);

					w11_R=DeltaR_had(taupt11_tmp,deltaR11_tmp)*DeltaR_lep(taupt21_tmp,deltaR21_tmp);
					w12_R=DeltaR_had(taupt11_tmp,deltaR11_tmp)*DeltaR_lep(taupt22_tmp,deltaR22_tmp);
					w21_R=DeltaR_had(taupt12_tmp,deltaR12_tmp)*DeltaR_lep(taupt21_tmp,deltaR21_tmp);
					w22_R=DeltaR_had(taupt12_tmp,deltaR12_tmp)*DeltaR_lep(taupt22_tmp,deltaR22_tmp);



					prob_metx=TMath::Exp(-pow(METX-metx,2)/(2*sig_metx*sig_metx));
					prob_mety=TMath::Exp(-pow(METY-mety,2)/(2*sig_mety*sig_mety));
					w11=w11_R*prob_metx*prob_mety;
					w12=w12_R*prob_metx*prob_mety;
					w21=w21_R*prob_metx*prob_mety;
					w22=w22_R*prob_metx*prob_mety;

					m11=(p4vis1+p4miss11 +p4vis2+ p4miss21).M();
					m12=(p4vis1+p4miss11 +p4vis2+ p4miss22).M();
					m21=(p4vis1+p4miss12 +p4vis2+ p4miss21).M();
					m22=(p4vis1+p4miss12 +p4vis2+ p4miss22).M();
					//std::cout<<"Mass "<<m11<< " "<<m12<< " "<<m21<<" "<<m22<<" Taupt "<<tau_pt_1<< " "<<tau_pt_2<<" TauE "<<tau_E_1<<" "<<tau_E_2<<std::endl;		 
					mass_dummy_all.Fill(m11,w11);
					mass_dummy_all.Fill(m12,w12);
					mass_dummy_all.Fill(m21,w21);
					mass_dummy_all.Fill(m22,w22);

					double tau_pt_1=(p4vis1+p4miss11).Pt();
					double tau_pt_2=(p4vis1+p4miss12).Pt();					
                                           taupt_dummy_all.Fill(tau_pt_1,w11);
                                           taupt_dummy_all.Fill(tau_pt_1,w12);
                                           taupt_dummy_all.Fill(tau_pt_2,w21);
                                           taupt_dummy_all.Fill(tau_pt_2,w22);
					
	//    Solution 1
					p4Z=(p4vis1+p4miss11 +p4vis2+ p4miss21);
					p4tau_had=p4vis1+p4miss11 ;
					p4tau_lep=p4vis2+p4miss21;
					spin_Obs_X1=p4vis1.Pt()/p4tau_had.Pt(); 

					double mass_t2=pow(p4tau_had.M(),2)+pow(p4tau_lep.M(),2)+2*(p4tau_had.Et()*p4tau_lep.Et()-p4tau_had.Px()*p4tau_lep.Px()-p4tau_had.Py()*p4tau_lep.Py());
		//			std::cout<<"Mass_T2 "<<	mass_t2<<" ";		
					spin_Obs_X2=2*p4vis1.Pt()/sqrt(mass_t2);
					p4tau_had.Boost(-p4Z.BoostVector());
					TLorentzVector p4vis1_Boosted= p4vis1;
					p4vis1_Boosted.Boost(-p4Z.BoostVector());
					spin_Obs_X3=p4vis1_Boosted.E()/p4tau_had.E(); 
					spin_Obs_X4=2*p4vis1_Boosted.E()/p4Z.M(); 					 

					hist_dummy_spin_Obs_X1.Fill(spin_Obs_X1,w11);
					hist_dummy_spin_Obs_X2.Fill(spin_Obs_X2,w11);
					hist_dummy_spin_Obs_X3.Fill(spin_Obs_X3,w11);
					hist_dummy_spin_Obs_X4.Fill(spin_Obs_X4,w11);

					//    Solution 2
					p4Z=(p4vis1+p4miss11 +p4vis2+ p4miss22);
					p4tau_had=p4vis1+p4miss11 ;
					p4tau_lep=p4vis2+ p4miss22;
					spin_Obs_X1=p4vis1.Pt()/p4tau_had.Pt(); 

					mass_t2=pow(p4tau_had.M(),2)+pow(p4tau_lep.M(),2)+2*(p4tau_had.Et()*p4tau_lep.Et()-p4tau_had.Px()*p4tau_lep.Px()-p4tau_had.Py()*p4tau_lep.Py());
		//			std::cout<< mass_t2<<" ";					
					spin_Obs_X2=2*p4vis1.Pt()/sqrt(mass_t2);
					p4tau_had.Boost(-p4Z.BoostVector());
					p4vis1_Boosted= p4vis1;
                                        p4vis1_Boosted.Boost(-p4Z.BoostVector());
					spin_Obs_X3=p4vis1_Boosted.E()/p4tau_had.E(); 
					spin_Obs_X4=2*p4vis1_Boosted.E()/p4Z.M(); 					 

					hist_dummy_spin_Obs_X1.Fill(spin_Obs_X1,w12);
					hist_dummy_spin_Obs_X2.Fill(spin_Obs_X2,w12);
					hist_dummy_spin_Obs_X3.Fill(spin_Obs_X3,w12);
					hist_dummy_spin_Obs_X4.Fill(spin_Obs_X4,w12);
					//   Solution 3
					p4Z=(p4vis1+p4miss12 +p4vis2+ p4miss21);
					p4tau_had=p4vis1+p4miss12 ;
					p4tau_lep=p4vis2+ p4miss21;
					spin_Obs_X1=p4vis1.Pt()/p4tau_had.Pt(); 

					mass_t2=pow(p4tau_had.M(),2)+pow(p4tau_lep.M(),2)+2*(p4tau_had.Et()*p4tau_lep.Et()-p4tau_had.Px()*p4tau_lep.Px()-p4tau_had.Py()*p4tau_lep.Py());
		//			std::cout<< mass_t2<<" ";					
					spin_Obs_X2=2*p4vis1.Pt()/sqrt(mass_t2);
					p4tau_had.Boost(-p4Z.BoostVector());
					p4vis1_Boosted= p4vis1;					
					p4vis1_Boosted.Boost(-p4Z.BoostVector());
					spin_Obs_X3=p4vis1_Boosted.E()/p4tau_had.E(); 
					spin_Obs_X4=2*p4vis1_Boosted.E()/p4Z.M(); 					 

					hist_dummy_spin_Obs_X1.Fill(spin_Obs_X1,w21);
					hist_dummy_spin_Obs_X2.Fill(spin_Obs_X2,w21);
					hist_dummy_spin_Obs_X3.Fill(spin_Obs_X3,w21);
					hist_dummy_spin_Obs_X4.Fill(spin_Obs_X4,w21);

					//   Solution 4
                                        p4Z=(p4vis1+p4miss12 +p4vis2+ p4miss22);
                                        p4tau_had=p4vis1+p4miss12 ;
                                        p4tau_lep=p4vis2+ p4miss22;
                                        spin_Obs_X1=p4vis1.Pt()/p4tau_had.Pt();

					mass_t2=pow(p4tau_had.M(),2)+pow(p4tau_lep.M(),2)+2*(p4tau_had.Et()*p4tau_lep.Et()-p4tau_had.Px()*p4tau_lep.Px()-p4tau_had.Py()*p4tau_lep.Py());
		//			std::cout<< mass_t2<<"\n"; 
					spin_Obs_X2=2*p4vis1.Pt()/sqrt(mass_t2);

					p4tau_had.Boost(-p4Z.BoostVector());
                                        p4vis1_Boosted = p4vis1;
                                        p4vis1_Boosted.Boost(-p4Z.BoostVector());
                                        spin_Obs_X3=p4vis1_Boosted.E()/p4tau_had.E();
                                        spin_Obs_X4=2*p4vis1_Boosted.E()/p4Z.M();

                                        hist_dummy_spin_Obs_X1.Fill(spin_Obs_X1,w22);
                                        hist_dummy_spin_Obs_X2.Fill(spin_Obs_X2,w22);
                                        hist_dummy_spin_Obs_X3.Fill(spin_Obs_X3,w22);
                                        hist_dummy_spin_Obs_X4.Fill(spin_Obs_X4,w22);


       
					
					check=true;

				      }//m_miss2
			            }//metx and mety
			      }//phi_miss2
			} //phi_miss1
      
     
        if (!check) {std::cout<<"No mmc solution"<<std::endl;return check; }
         //All Solution
	

         bool mass_dum_count_sw=false;


        mass_z=-999;
	if (mass_dummy_all.Integral()>0 && mass_dummy_all.GetMaximumBin()<mass_dummy_all.GetNbinsX()-1){
		mass_dum_count++;
		mass_dum_count_sw=true;                      
                double sum_mass=0,sum_weight=0;
                int maxbin_high=(mass_dummy_all.GetMaximumBin()>mass_dummy_all.GetNbinsX()-2)?mass_dummy_all.GetMaximumBin()+1:mass_dummy_all.GetMaximumBin()+3;

                for (int i=mass_dummy_all.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=mass_dummy_all.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=mass_dummy_all.GetBinCenter(i);
                        sum_mass=sum_mass+mss*iweight;
                }
                mass_z=sum_mass/sum_weight;
        }
        
        Spin_Obs_X1=-999;
        if (hist_dummy_spin_Obs_X1.Integral()>0 && hist_dummy_spin_Obs_X1.GetMaximumBin()<hist_dummy_spin_Obs_X1.GetNbinsX()-1){
                double sum_spin_X=0,sum_weight=0;
                int maxbin_high=(hist_dummy_spin_Obs_X1.GetMaximumBin()>hist_dummy_spin_Obs_X1.GetNbinsX()-2)?hist_dummy_spin_Obs_X1.GetMaximumBin()+1:hist_dummy_spin_Obs_X1.GetMaximumBin()+3;

                for (int i=hist_dummy_spin_Obs_X1.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=hist_dummy_spin_Obs_X1.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=hist_dummy_spin_Obs_X1.GetBinCenter(i);
                        sum_spin_X=sum_spin_X+mss*iweight;
                }
                Spin_Obs_X1=sum_spin_X/sum_weight;
   
        }
        
        Spin_Obs_X2=-999;
        if (hist_dummy_spin_Obs_X2.Integral()>0 && hist_dummy_spin_Obs_X2.GetMaximumBin()<hist_dummy_spin_Obs_X2.GetNbinsX()-1){
                double sum_spin_X=0,sum_weight=0;
                int maxbin_high=(hist_dummy_spin_Obs_X2.GetMaximumBin()>hist_dummy_spin_Obs_X2.GetNbinsX()-2)?hist_dummy_spin_Obs_X2.GetMaximumBin()+1:hist_dummy_spin_Obs_X2.GetMaximumBin()+3;

                for (int i=hist_dummy_spin_Obs_X2.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=hist_dummy_spin_Obs_X2.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=hist_dummy_spin_Obs_X2.GetBinCenter(i);
                        sum_spin_X=sum_spin_X+mss*iweight;
                }
               Spin_Obs_X2=sum_spin_X/sum_weight;
		
		
		
        }
       Spin_Obs_X3=-999;
        if (hist_dummy_spin_Obs_X3.Integral()>0 && hist_dummy_spin_Obs_X3.GetMaximumBin()<hist_dummy_spin_Obs_X3.GetNbinsX()-1){
                double sum_spin_X=0,sum_weight=0;
                int maxbin_high=(hist_dummy_spin_Obs_X3.GetMaximumBin()>hist_dummy_spin_Obs_X3.GetNbinsX()-2)?hist_dummy_spin_Obs_X3.GetMaximumBin()+1:hist_dummy_spin_Obs_X3.GetMaximumBin()+3;
                
                for (int i=hist_dummy_spin_Obs_X3.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=hist_dummy_spin_Obs_X3.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=hist_dummy_spin_Obs_X3.GetBinCenter(i);
                        sum_spin_X=sum_spin_X+mss*iweight;
                }
               Spin_Obs_X3=sum_spin_X/sum_weight;
                
                
                
        }

Spin_Obs_X4=-999;
        if (hist_dummy_spin_Obs_X4.Integral()>0 && hist_dummy_spin_Obs_X4.GetMaximumBin()<hist_dummy_spin_Obs_X4.GetNbinsX()-1){
                double sum_spin_X=0,sum_weight=0;
                int maxbin_high=(hist_dummy_spin_Obs_X4.GetMaximumBin()>hist_dummy_spin_Obs_X4.GetNbinsX()-2)?hist_dummy_spin_Obs_X4.GetMaximumBin()+1:hist_dummy_spin_Obs_X4.GetMaximumBin()+3;
                
                for (int i=hist_dummy_spin_Obs_X4.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=hist_dummy_spin_Obs_X4.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=hist_dummy_spin_Obs_X4.GetBinCenter(i);
                        sum_spin_X=sum_spin_X+mss*iweight;
                }
               Spin_Obs_X4=sum_spin_X/sum_weight;
                
                
                
        }

tau_pt=-999;
        if (taupt_dummy_all.Integral()>0 && taupt_dummy_all.GetMaximumBin()<taupt_dummy_all.GetNbinsX()-1){
                double sum_taupt=0,sum_weight=0;
                int maxbin_high=(taupt_dummy_all.GetMaximumBin()>taupt_dummy_all.GetNbinsX()-2)?taupt_dummy_all.GetMaximumBin()+1:taupt_dummy_all.GetMaximumBin()+3;

                for (int i=taupt_dummy_all.GetMaximumBin()-2;i<maxbin_high;i++){
                        double iweight=taupt_dummy_all.GetBinContent(i);
                        sum_weight=sum_weight+iweight;
                        double mss=taupt_dummy_all.GetBinCenter(i);
                        sum_taupt=sum_taupt+mss*iweight;
                }
                tau_pt=sum_taupt/sum_weight;

        }

//std::cout<<"......................."<<Diff_Spin_obs_X1.at(tt).GetName()<<std::endl;
if ( mass_dum_count%2==0 && mass_dum_count_sw){
                TString tmp_str;

		TFile fil("mass_Z.root","UPDATE");
		tmp_str=mass_dummy_all.GetName();tmp_str+=mass_dum_count;
		mass_dummy_all.SetName(tmp_str); mass_dummy_all.Write();

		tmp_str=hist_dummy_spin_Obs_X1.GetName();tmp_str+=mass_dum_count;
		hist_dummy_spin_Obs_X1.SetName(tmp_str);hist_dummy_spin_Obs_X1.Write();

		tmp_str=hist_dummy_spin_Obs_X2.GetName();tmp_str+=mass_dum_count;
		hist_dummy_spin_Obs_X2.SetName(tmp_str);hist_dummy_spin_Obs_X2.Write();

		tmp_str=hist_dummy_spin_Obs_X3.GetName();tmp_str+=mass_dum_count;
		hist_dummy_spin_Obs_X3.SetName(tmp_str);hist_dummy_spin_Obs_X3.Write();

		tmp_str=hist_dummy_spin_Obs_X4.GetName();tmp_str+=mass_dum_count;
		hist_dummy_spin_Obs_X4.SetName(tmp_str);hist_dummy_spin_Obs_X4.Write();

fil.Delete("mmc_polarizedtau_default_SPIN_OBS_X1MC_tautau_DY;1");
fil.Delete("mmc_polarizedtau_default_SPIN_OBS_X2MC_tautau_DY;1");
fil.Delete("mmc_polarizedtau_default_SPIN_OBS_X3MC_tautau_DY;1");

fil.Delete("mmc_polarizedtau_default-diff_Spin_Obs_X1MC_tautau_DY;1");
fil.Delete("mmc_polarizedtau_default-diff_Spin_Obs_X2MC_tautau_DY;1");
fil.Delete("mmc_polarizedtau_default-diff_Spin_Obs_X3MC_tautau_DY;1");

H_SPIN_OBS_X1.at(tt).Write();
H_SPIN_OBS_X2.at(tt).Write();
H_SPIN_OBS_X3.at(tt).Write();

Diff_Spin_obs_X1.at(tt).Write();
Diff_Spin_obs_X2.at(tt).Write();
Diff_Spin_obs_X3.at(tt).Write();


/*
		//tmp_str=Diff_Spin_obs_X1.at(t).GetName();tmp_str+=mass_dum_count;
		//Diff_Spin_obs_X1.at(t).SetName(tmp_str);
		Diff_Spin_obs_X1.at(tt).Write();
		//tmp_str=Diff_Spin_obs_X2.at(t).GetName();tmp_str+=mass_dum_count;
		//Diff_Spin_obs_X2.at(t).SetName(tmp_str);
		Diff_Spin_obs_X2.at(tt).Write();
		//tmp_str=Diff_Spin_obs_X3.at(t).GetName();tmp_str+=mass_dum_count;
		//Diff_Spin_obs_X3.at(t).SetName(tmp_str);
		Diff_Spin_obs_X3.at(tt).Write();

		//tmp_str=H_SPIN_OBS_X1.at(t).GetName();tmp_str+=mass_dum_count;
		//H_SPIN_OBS_X1.at(t).SetName(tmp_str);
		H_SPIN_OBS_X1.at(tt).Write();

		//tmp_str=H_SPIN_OBS_X2.at(t).GetName();tmp_str+=mass_dum_count;
	//	H_SPIN_OBS_X2.at(t).SetName(tmp_str);
	        H_SPIN_OBS_X2.at(tt).Write();
		//tmp_str=H_SPIN_OBS_X3.at(t).GetName();tmp_str+=mass_dum_count;
	//	H_SPIN_OBS_X3.at(t).SetName(tmp_str);
	         H_SPIN_OBS_X3.at(tt).Write();
	//	tmp_str=H_SPIN_OBS_X4.at(t).GetName();tmp_str+=mass_dum_count;
	//	H_SPIN_OBS_X4.at(t).SetName(tmp_str);
	         H_SPIN_OBS_X4.at(tt).Write();
*/
fil.Close();
}
  return check;
}
void  MMC_PolarizedTau::Finish() {
  
  Selection::Finish();
}







     






















