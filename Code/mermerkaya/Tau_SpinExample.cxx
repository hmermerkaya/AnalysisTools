#include "Tau_SpinExample.h"  
#include "TLorentzVector.h"   
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "TF1.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSpinerInterface.h"
                 
Tau_SpinExample::Tau_SpinExample(TString Name_, TString id_):
  Selection(Name_,id_)
{   
  //verbose=true;
}

Tau_SpinExample::~Tau_SpinExample() {
  for(unsigned int j=0; j<Npassed.size(); j++) {
    std::cout << "Tau_SpinExample::~Tau_SpinExample Selection Summary before: " 
	      << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	      << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Tau_SpinExample::~Tau_SpinExample()" << std::endl;
}
      
void  Tau_SpinExample::Configure() {
  // Setup Cut Values
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
  //////////////////////////////////////
  // Gen_Infos /////// Lab_Frame ///////    
  // X1 = Pi_Pt / Tau_Pt
  pi_PtRatio=HConfig.GetTH1D(Name+"_pi_PtRatio","pi_PtRatio",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}","Events");
  pi_PtRatio_hplus=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus","pi_PtRatio_hplus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}","Events");
  pi_PtRatio_hminus=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus","pi_PtRatio_hminus",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{-}}","Events");

  pi_PtRatio_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_Spin","pi_PtRatio_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_spin","Events");
  pi_PtRatio_hplus_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_hplus_Spin","pi_PtRatio_hplus_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{+}}|_spin","Events");
  pi_PtRatio_hminus_Spin=HConfig.GetTH1D(Name+"_pi_PtRatio_hminus_Spin","pi_PtRatio_hminus_Spin",20,0.0,1.0,"P_{t,#pi}/P_{t,#tau}|_{h^{-}}|_spin","Events");
  // X2 = Pi_Pt / Z_Mt
  pi_PtoverZmt=HConfig.GetTH1D(Name+"_pi_PtoverZmt","pi_PtoverZmt",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}","Events");
  pi_PtoverZmt_hplus=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hplus","pi_PtoverZmt_hplus",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{+}}","Events");
  pi_PtoverZmt_hminus=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hminus","pi_PtoverZmt_hminus",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{-}}","Events");
  pi_PtoverZmt_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_Spin","pi_PtoverZmt_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_spin","Events");
  pi_PtoverZmt_hplus_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hplus_Spin","pi_PtoverZmt_hplus_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{+}}|_spin","Events");
  pi_PtoverZmt_hminus_Spin=HConfig.GetTH1D(Name+"_pi_PtoverZmt_hminus_Spin","pi_PtoverZmt_hminus_Spin",20,0.0,1.0,"P_{t,#pi}/M_{t,Z}|_{h^{-}}|_spin","Events");
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  pi_EoverEtau=HConfig.GetTH1D(Name+"_pi_EoverEtau","pi_EoverEtau",20,0.0,1.0,"E_{#pi}/E_{#tau}","Events");
  pi_EoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_EoverEtau_hplus","pi_EoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_EoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_EoverEtau_hminus","pi_EoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  pi_EoverEtau_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_Spin","pi_EoverEtau_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_spin","Events");
  pi_EoverEtau_hplus_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_hplus_Spin","pi_EoverEtau_hplus_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}|_spin","Events");
  pi_EoverEtau_hminus_Spin=HConfig.GetTH1D(Name+"_pi_EoverEtau_hminus_Spin","pi_EoverEtau_hminus_Spin",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}|_spin","Events");
  // X4 = Pi_E / Z_M
  pi_EoverZm=HConfig.GetTH1D(Name+"_pi_EoverZm","pi_EoverZm",20,0.0,1.0,"E_{#pi}/M_{Z}","Events");
  pi_EoverZm_hplus=HConfig.GetTH1D(Name+"_pi_EoverZm_hplus","pi_EoverZm_hplus",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{+}}","Events");
  pi_EoverZm_hminus=HConfig.GetTH1D(Name+"_pi_EoverZm_hminus","pi_EoverZm_hminus",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{-}}","Events");
  pi_EoverZm_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_Spin","pi_EoverZm_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_spin","Events");
  pi_EoverZm_hplus_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_hplus_Spin","pi_EoverZm_hplus_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{+}}|_spin","Events");
  pi_EoverZm_hminus_Spin=HConfig.GetTH1D(Name+"_pi_EoverZm_hminus_Spin","pi_EoverZm_hminus_Spin",20,0.0,1.0,"E_{#pi}/M_{Z}|_{h^{-}}|_spin","Events");
  
  ////////////////////////////////////// 
  // Rec_Infos /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  pirec_PtRatio=HConfig.GetTH1D(Name+"_pirec_PtRatio","pirec_PtRatio",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}","Events");
  pirec_PtRatio_hplus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hplus","pirec_PtRatio_hplus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}","Events");
  pirec_PtRatio_hminus=HConfig.GetTH1D(Name+"_pirec_PtRatio_hminus","pirec_PtRatio_hminus",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{-}}","Events");
  pirec_PtRatio_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_Spin","pirec_PtRatio_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_spin","Events");
  pirec_PtRatio_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_hplus_Spin","pirec_PtRatio_hplus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_PtRatio_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_PtRatio_hminus_Spin","pirec_PtRatio_hminus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/P_{t,#tau_{rec}}|_{h^{-}}|_spin","Events");
  
  // X2 = Pirec_Pt / Zrec_Mt
  pirec_PtoverZmt=HConfig.GetTH1D(Name+"_pirec_PtoverZmt","pirec_PtoverZmt",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}","Events");
  pirec_PtoverZmt_hplus=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hplus","pirec_PtoverZmt_hplus",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{+}}","Events");
  pirec_PtoverZmt_hminus=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hminus","pirec_PtoverZmt_hminus",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{-}}","Events");
  pirec_PtoverZmt_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_Spin","pirec_PtoverZmt_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_spin","Events");
  pirec_PtoverZmt_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hplus_Spin","pirec_PtoverZmt_hplus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{+}}|_spin","Events");
  pirec_PtoverZmt_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_PtoverZmt_hminus_Spin","pirec_PtoverZmt_hminus_Spin",20,0.0,1.0,"P_{t,#pi_{rec}}/M_{t,Z_{rec}}|_{h^{-}}|_spin","Events");
  
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  pirec_EoverEtau=HConfig.GetTH1D(Name+"_pirec_EoverEtau","pirec_EoverEtau",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}","Events");
  pirec_EoverEtau_hplus=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hplus","pirec_EoverEtau_hplus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}","Events");
  pirec_EoverEtau_hminus=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hminus","pirec_EoverEtau_hminus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}","Events");
  pirec_EoverEtau_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_Spin","pirec_EoverEtau_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_spin","Events");
  pirec_EoverEtau_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hplus_Spin","pirec_EoverEtau_hplus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_EoverEtau_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverEtau_hminus_Spin","pirec_EoverEtau_hminus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}|_spin","Events");
  
  // X4 = Pirec_E / Zrec_M
  pirec_EoverZm=HConfig.GetTH1D(Name+"_pirec_EoverZm","pirec_EoverZm",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}","Events");
  pirec_EoverZm_hplus=HConfig.GetTH1D(Name+"_pirec_EoverZm_hplus","pirec_EoverZm_hplus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}","Events");
  pirec_EoverZm_hminus=HConfig.GetTH1D(Name+"_pirec_EoverZm_hminus","pirec_EoverZm_hminus",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}","Events");
  pirec_EoverZm_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_Spin","pirec_EoverZm_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_spin","Events");
  pirec_EoverZm_hplus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_hplus_Spin","pirec_EoverZm_hplus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{+}}|_spin","Events");
  pirec_EoverZm_hminus_Spin=HConfig.GetTH1D(Name+"_pirec_EoverZm_hminus_Spin","pirec_EoverZm_hminus_Spin",20,0.0,1.0,"E_{#pi_{rec}}/E_{#tau_{rec}}|_{h^{-}}|_spin","Events");
  
  
  
  // Setup Extra Histograms
  // Reco Objects
  //Delta_R(Mu_rec,Mu_gen), Delta_R(Pi_rec,Pi_gen)
  dR_recMu_genMu=HConfig.GetTH1D(Name+"_dR_recMu_genMu","dR_recMu_genMu",100,0.,8.,"dR_recMu_genMu","Events");
  mindR_recMu_genMu=HConfig.GetTH1D(Name+"_mindR_recMu_genMu","mindR_recMu_genMu",100,0.,1.,"mindR_recMu_genMu","Events");
  dR_recPi_genPi=HConfig.GetTH1D(Name+"_dR_recPi_genPi","dR_recPi_genPi",100,0.,8.,"dR_recPi_genPi","Events");
  mindR_recPi_genPi=HConfig.GetTH1D(Name+"_mindR_recPi_genPi","mindR_recPi_genPi",100,0.,1.,"mindR_recPi_genPi","Events");
  //Z---> TauPirec
  M_TauPi_rec=HConfig.GetTH1D(Name+"_M_TauPi_rec","M_TauPi_rec",180,0.,4.,"M_TauPi_rec","Events");
  Pt_TauPi_rec=HConfig.GetTH1D(Name+"_Pt_TauPi_rec","Pt_TauPi_rec",30,0.,70.,"Pt_TauPi_rec","Events");
  Eta_TauPi_rec=HConfig.GetTH1D(Name+"_Eta_TauPi_rec","Eta_TauPi_rec",50,-2.5,2.5,"Eta_TauPi_rec","Events");
  Phi_TauPi_rec=HConfig.GetTH1D(Name+"_Phi_TauPi_rec","Phi_TauPi_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauPi_rec","Events");
  //Z---> TauMurec
  M_TauMu_rec=HConfig.GetTH1D(Name+"_M_TauMu_rec","M_TauMu_rec",180,0.,4.,"M_TauMu_rec","Events");
  Pt_TauMu_rec=HConfig.GetTH1D(Name+"_Pt_TauMu_rec","Pt_TauMu_rec",30,0.,70.,"Pt_TauMu_rec","Events");
  Eta_TauMu_rec=HConfig.GetTH1D(Name+"_Eta_TauMu_rec","Eta_TauMu_rec",50,-2.5,2.5,"Eta_TauMu_rec","Events");
  Phi_TauMu_rec=HConfig.GetTH1D(Name+"_Phi_TauMu_rec","Phi_TauMu_rec",32,-TMath::Pi(),TMath::Pi(),"Phi_TauMu_rec","Events");

  // extra Histos
  mu_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{+}}","Events");
  mu_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_mu_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#mu}/E_{#tau}|_{h^{-}}","Events");
  
  pi_ExoverEtau_hplus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hplus","ExoverEtau_hplus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{+}}","Events");
  pi_ExoverEtau_hminus=HConfig.GetTH1D(Name+"_pi_ExoverEtau_hminus","ExoverEtau_hminus",20,0.0,1.0,"E_{#pi}/E_{#tau}|_{h^{-}}","Events");
  


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
}

void  Tau_SpinExample::Store_ExtraDist() {
  ////////////////////////////////////// 
  // Gen_Infos /////// Lab_Frame ///////
  // X1 = Pi_Pt / Tau_Pt     
  Extradist1d.push_back(&pi_PtRatio);
  Extradist1d.push_back(&pi_PtRatio_hplus);
  Extradist1d.push_back(&pi_PtRatio_hminus);
 
  Extradist1d.push_back(&pi_PtRatio_Spin);
  Extradist1d.push_back(&pi_PtRatio_hplus_Spin);
  Extradist1d.push_back(&pi_PtRatio_hminus_Spin);
  // X2 = Pi_Pt / Z_Mt
  Extradist1d.push_back(&pi_PtoverZmt);
  Extradist1d.push_back(&pi_PtoverZmt_hplus);
  Extradist1d.push_back(&pi_PtoverZmt_hminus);
  Extradist1d.push_back(&pi_PtoverZmt_Spin);
  Extradist1d.push_back(&pi_PtoverZmt_hplus_Spin);
  Extradist1d.push_back(&pi_PtoverZmt_hminus_Spin);
  ////// Z_Rest_Frame //////
  // X3 = Pi_E / Tau_E
  Extradist1d.push_back(&pi_EoverEtau);
  Extradist1d.push_back(&pi_EoverEtau_hplus);
  Extradist1d.push_back(&pi_EoverEtau_hminus);
  Extradist1d.push_back(&pi_EoverEtau_Spin);
  Extradist1d.push_back(&pi_EoverEtau_hplus_Spin);
  Extradist1d.push_back(&pi_EoverEtau_hminus_Spin);
  // X4 = Pi_E / Z_M
  Extradist1d.push_back(&pi_EoverZm);
  Extradist1d.push_back(&pi_EoverZm_hplus);
  Extradist1d.push_back(&pi_EoverZm_hminus);
  Extradist1d.push_back(&pi_EoverZm_Spin);
  Extradist1d.push_back(&pi_EoverZm_hplus_Spin);
  Extradist1d.push_back(&pi_EoverZm_hminus_Spin);
  
  ////////////////////////////////////// 
  // Rec_Infos /////// Lab_Frame ///////
  // X1 = Pirec_Pt / Taurec_Pt
  Extradist1d.push_back(&pirec_PtRatio);
  Extradist1d.push_back(&pirec_PtRatio_hplus);
  Extradist1d.push_back(&pirec_PtRatio_hminus);
  Extradist1d.push_back(&pirec_PtRatio_Spin);
  Extradist1d.push_back(&pirec_PtRatio_hplus_Spin);
  Extradist1d.push_back(&pirec_PtRatio_hminus_Spin);
  
  // X2 = Extradist1d.Push_Back(&Pirec_Pt / Zrec_Mt
  Extradist1d.push_back(&pirec_PtoverZmt);
  Extradist1d.push_back(&pirec_PtoverZmt_hplus);
  Extradist1d.push_back(&pirec_PtoverZmt_hminus);
  Extradist1d.push_back(&pirec_PtoverZmt_Spin);
  Extradist1d.push_back(&pirec_PtoverZmt_hplus_Spin);
  Extradist1d.push_back(&pirec_PtoverZmt_hminus_Spin);
  
  ////// Z_Rest_Frame //////
  // X3 = Pirec_E / Taurec_E
  Extradist1d.push_back(&pirec_EoverEtau);
  Extradist1d.push_back(&pirec_EoverEtau_hplus);
  Extradist1d.push_back(&pirec_EoverEtau_hminus);
  Extradist1d.push_back(&pirec_EoverEtau_Spin);
  Extradist1d.push_back(&pirec_EoverEtau_hplus_Spin);
  Extradist1d.push_back(&pirec_EoverEtau_hminus_Spin);
  
  // X4 = Extradist1d.Push_Back(&Pirec_E / Zrec_M
  Extradist1d.push_back(&pirec_EoverZm);
  Extradist1d.push_back(&pirec_EoverZm_hplus);
  Extradist1d.push_back(&pirec_EoverZm_hminus);
  Extradist1d.push_back(&pirec_EoverZm_Spin);
  Extradist1d.push_back(&pirec_EoverZm_hplus_Spin);
  Extradist1d.push_back(&pirec_EoverZm_hminus_Spin);
  
  //Reco objects
  Extradist1d.push_back(&M_TauPi_rec);        Extradist1d.push_back(&Pt_TauPi_rec);
  Extradist1d.push_back(&Eta_TauPi_rec);      Extradist1d.push_back(&Phi_TauPi_rec);
  Extradist1d.push_back(&M_TauMu_rec);        Extradist1d.push_back(&Pt_TauMu_rec);
  Extradist1d.push_back(&Eta_TauMu_rec);      Extradist1d.push_back(&Phi_TauMu_rec);        
  
  Extradist1d.push_back(&dR_recMu_genMu);     Extradist1d.push_back(&dR_recPi_genPi);
  Extradist1d.push_back(&mindR_recMu_genMu);  Extradist1d.push_back(&mindR_recPi_genPi);

  // Extra Histos
  Extradist1d.push_back(&mu_ExoverEtau_hplus);Extradist1d.push_back(&mu_ExoverEtau_hminus);
  Extradist1d.push_back(&pi_ExoverEtau_hplus);Extradist1d.push_back(&pi_ExoverEtau_hminus);

}

void  Tau_SpinExample::doEvent() {
  unsigned int t(0);
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)) {t=2;}// std::cout << "failed to find id" <<std::endl; return;}
  unsigned int Boson_idx,Boson3_idx,Boson4_idx,tau1_idx(0),tau2_idx(0),tau3_idx(0),tau4_idx(0);
  value.at(isZtautautopimu)=0;
  pass.at(isZtautautopimu) = true;
  if(pass.at(isZtautautopimu))value.at(isZtautautopimu)=1;
  
  double wobs=1;
  double w=1;
  
  //  if(verbose)  std::cout << "MCID=="<< Ntp->GetMCID() << "||Size==" << Npassed.size() << "||t== " << t << "||BosonIdx== " << Boson_idx << "||Tau1== " << tau1_idx << "||Tau2== " << tau1_idx << "||Ntaus== " << Ntp->NMCTaus() << std::endl;
  
  bool status=AnalysisCuts(t,w,wobs); 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status) {
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
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
    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    std::cout << "Spin_WT ==" << Spin_WT << std::endl;
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=1-hplus;//Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
    std::cout << "hplus " << hplus << " hminus " << hminus << std::endl;
   
  if (isnan(Spin_WT)) return;


 
    if(verbose)std::cout<< "A " <<std::endl;
    ////////////////////////////////////////////////
    //
    // Spin Validation
    //  
    ////////////////////////////////////////////////
    
    /////////////////////////////////    
    //
    //  Muon && Pion  Channel
    //
    /////////////////////////////////
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson3_idx,TauDecay::JAK_MUON,tau3_idx) && Ntp->hasSignalTauDecay(PDGInfo::Z0,Boson4_idx,TauDecay::JAK_PION,tau4_idx)) { //JAK_MUON && JAK_PION
      if(verbose)std::cout<< "pi-mu " <<std::endl;
      
      TLorentzVector GenZ=Ntp->MCSignalParticle_p4(Boson3_idx);
      TLorentzVector GenTauMu(0,0,0,0);
      TLorentzVector GenMu(0,0,0,0);  
      TLorentzVector GenNum(0,0,0,0);     
      TLorentzVector GenNutm(0,0,0,0);
      TLorentzVector Muon(0,0,0,0);
      TLorentzVector TauMuRec(0,0,0,0);
      TLorentzVector MuonRec(0,0,0,0);
      TLorentzVector Muon_mdR(0,0,0,0);
      
      TLorentzVector GenZ1=Ntp->MCSignalParticle_p4(Boson4_idx);
      TLorentzVector GenTauPi(0,0,0,0);
      TLorentzVector GenPi(0,0,0,0);
      TLorentzVector GenNutp(0,0,0,0);
      TLorentzVector Pion(0,0,0,0);
      TLorentzVector Zrec(0,0,0,0);
      TLorentzVector TauPiRec(0,0,0,0);
      TLorentzVector PionRec(0,0,0,0);
      TLorentzVector Pion_mdR(0,0,0,0);  
      
      double dR_Muons = -1;   double dR_Pions = -1;
      double mindRMu = 999;   double mindRPi=999; //dR_Muon < 0.1 && dR_Pion<0.4
      int mindRMu_Index = -1; int mindRPi_Index = -1;
      double Zrec_Mt = -1;    double Zgen_Mt = -1;
      
      //Gen Infos      
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
      } //end for loop 
      //Reco Muon
      for(unsigned int i=0; i < Ntp->NMuons(); i++) {
	dR_Muons = Ntp->Muon_p4(i).DeltaR(GenMu) ;
	dR_recMu_genMu.at(t).Fill(dR_Muons,w);  	  
	if(dR_Muons < mindRMu) { mindRMu = dR_Muons;  mindRMu_Index = i; }
	MuonRec = Ntp->Muon_p4(mindRMu_Index);
      }
      if(mindRMu<0.1 && mindRMu_Index>-1) { 
	mindR_recMu_genMu.at(t).Fill(mindRMu,w);   
	std::cout<< "Moun_minDR==" << mindRMu <<std::endl;
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
	dR_Pions = Ntp->PFTau_p4(i).DeltaR(GenPi) ;
	dR_recPi_genPi.at(t).Fill(dR_Pions,w);
	if(dR_Pions < mindRPi) { mindRPi = dR_Pions;  mindRPi_Index = i; }
	PionRec = Ntp->PFTau_p4(mindRPi_Index);
      }
      if(mindRPi<0.4 && mindRPi_Index>-1) { 
	mindR_recPi_genPi.at(t).Fill(mindRPi,w);
	Pion_mdR=PionRec;
      }     
      
      // Transvers mass of Z boson
      //GenZMt= sqrt(2*pttaum*pttaupi*(1-cos(Phitaum -Phipi )));
      //MT=sqrt(2*(Ntp->MET_CorrMVA_et())*Ntp->Muon_p4(mu_idx.at(0)).Pt()*fabs(1-cos(Ntp->Muon_p4(mu_idx.at(0)).Phi()-Ntp->MET_CorrMVA_phi())));
      
      /////////////////////////////////
      //// Tau_gen-XPlots ////
      /////////////////////////////////
      if(GenTauPi.E()>0) {
	//Lab-Frame
	// X1 = Pi_Pt / Tau_Pt  
	std::cout<<" SpinWt=="<< Spin_WT<< std::endl;
	pi_PtRatio.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w);
      	pi_PtRatio_hplus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus);
        pi_PtRatio_hminus.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus);

	pi_PtRatio_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*Spin_WT);
      	pi_PtRatio_hplus_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hplus*Spin_WT);
        pi_PtRatio_hminus_Spin.at(t).Fill(GenPi.Pt()/GenTauPi.Pt(),w*hminus*Spin_WT);
	
	// X2 = Pi_Pt / Z_Mt 
	Zgen_Mt = sqrt(2*GenTauMu.Pt()*GenTauPi.Pt()*(1-cos(GenTauMu.Phi() - GenTauPi.Phi())));	
	std::cout<<"Z_Mt=="<< Zgen_Mt<<std::endl;
	//std::cout<<"Z_Mt1=="<< GenZ1.Mt()<<std::endl;
	pi_PtoverZmt.at(t).Fill(GenPi.Pt()/Zgen_Mt,w);
	pi_PtoverZmt_hplus.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hplus);
	pi_PtoverZmt_hminus.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hminus); 

	pi_PtoverZmt_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*Spin_WT);
	pi_PtoverZmt_hplus_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hplus*Spin_WT);
	pi_PtoverZmt_hminus_Spin.at(t).Fill(GenPi.Pt()/Zgen_Mt,w*hminus*Spin_WT);
	// Z-Rest-Frame
	GenTauPi.Boost(-GenZ1.BoostVector());
        GenPi.Boost(-GenZ1.BoostVector());
        GenZ1.Boost(-GenZ1.BoostVector());
        // Now fill results                                   
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
	// X3 = Pi_E / Tau_E 
      	pi_EoverEtau.at(t).Fill(GenPi.E()/GenTauPi.E(),w);
	pi_EoverEtau_hplus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus);
      	pi_EoverEtau_hminus.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus);
	pi_EoverEtau_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*Spin_WT);
	pi_EoverEtau_hplus_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hplus*Spin_WT);
      	pi_EoverEtau_hminus_Spin.at(t).Fill(GenPi.E()/GenTauPi.E(),w*hminus*Spin_WT);
	// X4 = Pi_E / Z_M
	std::cout<<"Z_M =="<< GenZ1.M()<<std::endl;
	pi_EoverZm.at(t).Fill(GenPi.Pt()/GenZ1.M(),w);
	pi_EoverZm_hplus.at(t).Fill(GenPi.Pt()/GenZ1.M(),w*hplus);
	pi_EoverZm_hminus.at(t).Fill(GenPi.Pt()/GenZ1.M(),w*hminus); 
	pi_EoverZm_Spin.at(t).Fill(GenPi.Pt()/GenZ1.M(),w*Spin_WT);
	pi_EoverZm_hplus_Spin.at(t).Fill(GenPi.Pt()/GenZ1.M(),w*hplus*Spin_WT);
	pi_EoverZm_hminus_Spin.at(t).Fill(GenPi.Pt()/GenZ1.M(),w*hminus*Spin_WT);    
      }
      
      /////////////////////////////////
      //// Tau_rec  reconstruction ////
      /////////////////////////////////
      if(mindRMu < 0.1 && mindRPi < 0.4 && mindRMu_Index >-1 && mindRPi_Index >-1) {
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
      }
      if(TauPiRec.E()>0) {  
	//Lab-Frame
	// X1 = Pirec_Pt / Tau_Pt  
	std::cout<<" SpinWt(rec)=="<< Spin_WT<< std::endl;
	pirec_PtRatio.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w);
      	pirec_PtRatio_hplus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hplus);
        pirec_PtRatio_hminus.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hminus);
	pirec_PtRatio_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*Spin_WT);
      	pirec_PtRatio_hplus_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hplus*Spin_WT);
        pirec_PtRatio_hminus_Spin.at(t).Fill(Pion_mdR.Pt()/TauPiRec.Pt(),w*hminus*Spin_WT);
	// X2 = Pirec_Pt / Z_Mt 
	Zrec_Mt = sqrt(2.*TauMuRec.Pt()*TauPiRec.Pt()*(1-cos(TauMuRec.Phi() - TauPiRec.Phi())));
	std::cout<<"Zrec_Mt=="<< Zrec_Mt<< std::endl;
	//std::cout<<"Zrec_Mt=="<< Zrec.Mt()<<std::endl;
	
	pirec_PtoverZmt.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w);
	pirec_PtoverZmt_hplus.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hplus);
	pirec_PtoverZmt_hminus.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hminus); 
	pirec_PtoverZmt_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*Spin_WT);
	pirec_PtoverZmt_hplus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hplus*Spin_WT);
	pirec_PtoverZmt_hminus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec_Mt,w*hminus*Spin_WT);
	
	// Z-Rest-Frame
	TauPiRec.Boost(-Zrec.BoostVector());
        Pion_mdR.Boost(-Zrec.BoostVector());
        Zrec.Boost(-Zrec.BoostVector());
        // Now fill results                                   
        Ntp->TauSpinerSetSignal(Ntp->MCTau_charge(tau4_idx));
	// X3 = Pirec_E / Tau_E 
      	pirec_EoverEtau.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w);
	pirec_EoverEtau_hplus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hplus);
      	pirec_EoverEtau_hminus.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hminus);
	pirec_EoverEtau_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*Spin_WT);
	pirec_EoverEtau_hplus_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hplus*Spin_WT);
      	pirec_EoverEtau_hminus_Spin.at(t).Fill(Pion_mdR.E()/TauPiRec.E(),w*hminus*Spin_WT);
	// X4 = Pirec_E / Z_M
	//double Z1rec_M = TauMuRec.M() + TauPiRec.M();
	//std::cout<<"Z1rec_M =="<< Z1rec_M<<std::endl;
	std::cout<<"Zrec_M =="<< Zrec.M()<< std::endl;
	
	pirec_EoverZm.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w);
	pirec_EoverZm_hplus.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w*hplus);
	pirec_EoverZm_hminus.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w*hminus); 
	pirec_EoverZm_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w*Spin_WT);
	pirec_EoverZm_hplus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w*hplus*Spin_WT);
	pirec_EoverZm_hminus_Spin.at(t).Fill(Pion_mdR.Pt()/Zrec.M(),w*hminus*Spin_WT);    
      }
      
      
    } //end JAK_MUON && JAK_PION   ////////////   end Muon && Pion    ////////////    
  } //end if(Status)
} //end doEvent()




void  Tau_SpinExample::Finish() {
  
  Selection::Finish();
}







     






















