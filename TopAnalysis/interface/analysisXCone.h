#ifndef _analysisXCone_h_
#define _analysisXCone_h_

#include "TH1.h"
#include "TString.h"
#include "TFile.h"

#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"

#include <map>

struct analysisXCone_t
{
  Int_t run,event,lumi,cat;
  Int_t nvtx;
  Int_t nw;
  Int_t NJettiness;
  Int_t NJettSlope;
  Float_t weight[400];
  Int_t passSel,gen_passSel;
  Int_t nj[11],  nb[11];
  Int_t gen_nj, gen_nb;
  Float_t ptpos[11],  phipos[11],  ptll[11],  phill[11],  mll[11],  sumpt[11],  dphill[11];
  Float_t gen_ptpos, gen_phipos, gen_ptll, gen_phill, gen_mll, gen_sumpt, gen_dphill;
  Float_t ptttbar[11],     phittbar[11];  
  Float_t gen_ptttbar,    gen_phittbar;  
  Int_t n,gen_n,id[1000];
  Float_t pt[1000],     eta[1000],     phi[1000]; 
  Int_t gen_id[1000],gen_rec[1000]; 
  Float_t gen_pt[1000], gen_eta[1000], gen_phi[1000]; 
};

void createanalysisXConeTree(TTree *t,analysisXCone_t &tue);
void resetanalysisXCone(analysisXCone_t &tue);
void RunanalysisXCone(TString filename,
	      TString outname,
	      Int_t channelSelection, 
	      Int_t chargeSelection, 
	      SelectionTool::FlavourSplitting flavourSplitting,
	      TH1F *normH, 
	      Bool_t runSysts,
	      TString era,
	      Bool_t debug);
#endif
