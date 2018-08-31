#ifndef _topjetshape_h_
#define _topjetshape_h_

#include "TLorentzVector.h"
#include "TopLJets2015/TopAnalysis/interface/ObjectTools.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

void RunTopJetShape(TString filename,
		    TString outname,
		    Int_t channelSelection, 
		    Int_t chargeSelection, 
		    SelectionTool::FlavourSplitting flavourSplitting,
		    TH1F *normH, 
		    Bool_t runSysts,
		    std::string systVar,
		    TString era,
		    Bool_t debug=false);

struct TopJetShapeEvent_t
{
  Int_t period, nw, nl, ngl, nj, ngj;
  Int_t gen_sel, reco_sel;
  
  Float_t weight[1000];
  
  Float_t l_pt[5], l_eta[5], l_phi[5], l_m[5];
  Int_t l_id[5];
  
  Float_t gl_pt[5], gl_eta[5], gl_phi[5], gl_m[5];
  Int_t gl_id[5];
  
  Float_t j_pt[50], j_eta[50], j_phi[50], j_m[50];
  Int_t j_flavor[50], j_partonflavor[50], j_overlap[50], j_gj[50];
  
  Float_t j_charged_hadron_fraction[50];
  Float_t j_neutral_hadron_fraction[50];
  Float_t j_photon_fraction[50];
  Float_t j_lepton_fraction[50];
  
  Float_t gj_pt[50], gj_eta[50], gj_phi[50], gj_m[50];
  Int_t gj_flavor[50], gj_partonflavor[50], gj_overlap[50], gj_j[50];
  
  Float_t met_pt,met_phi;
  
  Float_t rho_ch, grho_ch;
  
  Float_t j_mult_charged[50];
  Float_t j_width_charged[50];
  Float_t j_ptd_charged[50];
  Float_t j_ptds_charged[50];
  Float_t j_ecc_charged[50];
  Float_t j_tau21_charged[50];
  Float_t j_tau32_charged[50];
  Float_t j_tau43_charged[50];
  Float_t j_zg_charged[50];
  Float_t j_zgxdr_charged[50];
  Float_t j_zgdr_charged[50];
  Float_t j_ga_width_charged[50];
  Float_t j_ga_lha_charged[50];
  Float_t j_ga_thrust_charged[50];
  Float_t j_c1_02_charged[50];
  Float_t j_c1_00_charged[50];
  Float_t j_c1_05_charged[50];
  Float_t j_c1_10_charged[50];
  Float_t j_c1_20_charged[50];
  Float_t j_c2_02_charged[50];
  Float_t j_c2_00_charged[50];
  Float_t j_c2_05_charged[50];
  Float_t j_c2_10_charged[50];
  Float_t j_c2_20_charged[50];
  Float_t j_c3_00_charged[50];
  Float_t j_c3_02_charged[50];
  Float_t j_c3_05_charged[50];
  Float_t j_c3_10_charged[50];
  Float_t j_c3_20_charged[50];
  Float_t j_m2_b1_charged[50];
  Float_t j_n2_b1_charged[50];
  Float_t j_n3_b1_charged[50];
  Float_t j_m2_b2_charged[50];
  Float_t j_n2_b2_charged[50];
  Float_t j_n3_b2_charged[50];
  Float_t j_nsd_charged[50];
  Float_t j_mult_all[50];
  Float_t j_width_all[50];
  Float_t j_ptd_all[50];
  Float_t j_ptds_all[50];
  Float_t j_ecc_all[50];
  Float_t j_tau21_all[50];
  Float_t j_tau32_all[50];
  Float_t j_tau43_all[50];
  Float_t j_zg_all[50];
  Float_t j_zgxdr_all[50];
  Float_t j_zgdr_all[50];
  Float_t j_ga_width_all[50];
  Float_t j_ga_lha_all[50];
  Float_t j_ga_thrust_all[50];
  Float_t j_c1_02_all[50];
  Float_t j_c1_00_all[50];
  Float_t j_c1_05_all[50];
  Float_t j_c1_10_all[50];
  Float_t j_c1_20_all[50];
  Float_t j_c2_02_all[50];
  Float_t j_c2_00_all[50];
  Float_t j_c2_05_all[50];
  Float_t j_c2_10_all[50];
  Float_t j_c2_20_all[50];
  Float_t j_c3_00_all[50];
  Float_t j_c3_02_all[50];
  Float_t j_c3_05_all[50];
  Float_t j_c3_10_all[50];
  Float_t j_c3_20_all[50];
  Float_t j_m2_b1_all[50];
  Float_t j_n2_b1_all[50];
  Float_t j_n3_b1_all[50];
  Float_t j_m2_b2_all[50];
  Float_t j_n2_b2_all[50];
  Float_t j_n3_b2_all[50];
  Float_t j_nsd_all[50];


  Float_t gj_mult_charged[50];
  Float_t gj_width_charged[50];
  Float_t gj_ptd_charged[50];
  Float_t gj_ptds_charged[50];
  Float_t gj_ecc_charged[50];
  Float_t gj_tau21_charged[50];
  Float_t gj_tau32_charged[50];
  Float_t gj_tau43_charged[50];
  Float_t gj_zg_charged[50];
  Float_t gj_zgxdr_charged[50];
  Float_t gj_zgdr_charged[50];
  Float_t gj_ga_width_charged[50];
  Float_t gj_ga_lha_charged[50];
  Float_t gj_ga_thrust_charged[50];
  Float_t gj_c1_02_charged[50];
  Float_t gj_c1_00_charged[50];
  Float_t gj_c1_05_charged[50];
  Float_t gj_c1_10_charged[50];
  Float_t gj_c1_20_charged[50];
  Float_t gj_c2_02_charged[50];
  Float_t gj_c2_00_charged[50];
  Float_t gj_c2_05_charged[50];
  Float_t gj_c2_10_charged[50];
  Float_t gj_c2_20_charged[50];
  Float_t gj_c3_00_charged[50];
  Float_t gj_c3_02_charged[50];
  Float_t gj_c3_05_charged[50];
  Float_t gj_c3_10_charged[50];
  Float_t gj_c3_20_charged[50];
  Float_t gj_m2_b1_charged[50];
  Float_t gj_n2_b1_charged[50];
  Float_t gj_n3_b1_charged[50];
  Float_t gj_m2_b2_charged[50];
  Float_t gj_n2_b2_charged[50];
  Float_t gj_n3_b2_charged[50];
  Float_t gj_nsd_charged[50];
  Float_t gj_mult_all[50];
  Float_t gj_width_all[50];
  Float_t gj_ptd_all[50];
  Float_t gj_ptds_all[50];
  Float_t gj_ecc_all[50];
  Float_t gj_tau21_all[50];
  Float_t gj_tau32_all[50];
  Float_t gj_tau43_all[50];
  Float_t gj_zg_all[50];
  Float_t gj_zgxdr_all[50];
  Float_t gj_zgdr_all[50];
  Float_t gj_ga_width_all[50];
  Float_t gj_ga_lha_all[50];
  Float_t gj_ga_thrust_all[50];
  Float_t gj_c1_02_all[50];
  Float_t gj_c1_00_all[50];
  Float_t gj_c1_05_all[50];
  Float_t gj_c1_10_all[50];
  Float_t gj_c1_20_all[50];
  Float_t gj_c2_02_all[50];
  Float_t gj_c2_00_all[50];
  Float_t gj_c2_05_all[50];
  Float_t gj_c2_10_all[50];
  Float_t gj_c2_20_all[50];
  Float_t gj_c3_00_all[50];
  Float_t gj_c3_02_all[50];
  Float_t gj_c3_05_all[50];
  Float_t gj_c3_10_all[50];
  Float_t gj_c3_20_all[50];
  Float_t gj_m2_b1_all[50];
  Float_t gj_n2_b1_all[50];
  Float_t gj_n3_b1_all[50];
  Float_t gj_m2_b2_all[50];
  Float_t gj_n2_b2_all[50];
  Float_t gj_n3_b2_all[50];
  Float_t gj_nsd_all[50];

};

double calcGA(double beta, double kappa, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev);
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev);

double getMult(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getPtD(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getPtDs(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getWidth(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getEcc(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getTau(int N, int M, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getC(int N, double beta, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
std::vector<double> getZg(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
std::map<TString,double> getECF(Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);
double getNSD(double zcut, double beta, Jet jet, bool includeNeutrals = false, bool usePuppi = false, double ptcut = 1.0);

double getPFFraction(std::vector<int> pids, Jet jet);

double deltaR(TLorentzVector v1, TLorentzVector v2);
double mapAngleMPiToPi(double phi);
#endif
