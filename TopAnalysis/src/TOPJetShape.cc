#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/TOP-16-006.h"
#include "TopLJets2015/TopAnalysis/interface/TOPJetShape.h"
#include "TopLJets2015/TopAnalysis/interface/BtagUncertaintyComputer.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
using namespace fastjet;
using namespace fastjet::contrib;

#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
using Rivet::Matrix;
using Rivet::EigenSystem;

using namespace std;


//
void RunTopJetShape(TString filename,
                 TString outname,
                 Int_t channelSelection, 
                 Int_t chargeSelection, 
                 FlavourSplitting flavourSplitting,
                 TH1F *normH, 
                 Bool_t runSysts)
{

  bool isTTbar( filename.Contains("_TTJets") );

  //prepare output
  TopJetShapeEvent_t tjsev;
  TString baseName=gSystem->BaseName(outname); 
  TString dirName=gSystem->DirName(outname);
  TFile *fOut=TFile::Open(dirName+"/"+baseName,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tjsev","tjsev");
  createTopJetShapeEventTree(outT,tjsev);
  outT->SetDirectory(fOut);

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *puTrue=(TH1 *)f->Get("analysis/putrue");
  puTrue->SetDirectory(0);
  puTrue->Scale(1./puTrue->Integral());
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev, true);
  Int_t nentries(t->GetEntriesFast());
  t->GetEntry(0);
  //bool requireEtriggerOnly(false);
  //if(ev.isData && filename.Contains("SingleElectron")) requireEtriggerOnly=true;

  cout << "...producing " << outname << " from " << nentries << " events" << endl;

  //PILEUP WEIGHTING
  std::vector<TGraph *>puWgtGr;
  if(!ev.isData)
    {
      TString puWgtUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/pileupWgts.root");
      gSystem->ExpandPathName(puWgtUrl);
      TFile *fIn=TFile::Open(puWgtUrl);
      TGraph *puData=(TGraph *)fIn->Get("pu_nom");
      Float_t totalData=puData->Integral();
      TH1 *tmp=(TH1 *)puTrue->Clone("tmp");
      for(Int_t xbin=1; xbin<=tmp->GetXaxis()->GetNbins(); xbin++)
        {
          Float_t yexp=puTrue->GetBinContent(xbin);
          Double_t xobs,yobs;
          puData->GetPoint(xbin-1,xobs,yobs);
          tmp->SetBinContent(xbin, yexp>0 ? yobs/(totalData*yexp) : 0. );
        }
      TGraph *gr=new TGraph(tmp);
      gr->SetName("puwgts_nom");
      puWgtGr.push_back( gr );
      tmp->Delete();
    }
    
  //LEPTON EFFICIENCIES
  TString lepEffUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/leptonEfficiencies.root");
  gSystem->ExpandPathName(lepEffUrl);
  std::map<TString,TH2 *> lepEffH;
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["m_sel"]=(TH2 *)fIn->Get("m_sel");
      lepEffH["m_trig"]=(TH2 *)fIn->Get("m_trig");      
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  lepEffUrl="${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root";
  gSystem->ExpandPathName(lepEffUrl);
  if(!ev.isData)
    {
      TFile *fIn=TFile::Open(lepEffUrl);
      lepEffH["e_sel"]=(TH2 *)fIn->Get("EGamma_SF2D");
      for(auto& it : lepEffH) it.second->SetDirectory(0);
      fIn->Close();
    }

  //B-TAG CALIBRATION
  TString btagUncUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/CSVv2.csv");
  gSystem->ExpandPathName(btagUncUrl);
  std::vector<BTagCalibrationReader *> sfbReaders, sflReaders;
  TString btagEffExpUrl("${CMSSW_BASE}/src/TopLJets2015/TopAnalysis/data/expTageff.root");
  gSystem->ExpandPathName(btagEffExpUrl);
  std::map<TString, TGraphAsymmErrors *> expBtagEff, expBtagEffPy8;
  BTagSFUtil myBTagSFUtil;
  if(!ev.isData)
    {
      BTagCalibration btvcalib("csvv2", btagUncUrl.Data());
      sfbReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "mujets", "central") );
      sflReaders.push_back( new BTagCalibrationReader(&btvcalib, BTagEntry::OP_MEDIUM, "incl", "central") );

      TFile *beffIn=TFile::Open(btagEffExpUrl);
      expBtagEffPy8["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEffPy8["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEffPy8["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
      
      TString btagExpPostFix("");
      if(isTTbar)
        {
          if(filename.Contains("_herwig")) btagExpPostFix="_herwig";
          if(filename.Contains("_scaleup")) btagExpPostFix="_scaleup";
          if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
        }
      btagEffExpUrl.ReplaceAll(".root",btagExpPostFix+".root");
      beffIn=TFile::Open(btagEffExpUrl);
      expBtagEff["b"]=(TGraphAsymmErrors *)beffIn->Get("b");
      expBtagEff["c"]=(TGraphAsymmErrors *)beffIn->Get("c");
      expBtagEff["udsg"]=(TGraphAsymmErrors *)beffIn->Get("udsg");
      beffIn->Close();
    }

  //BOOK HISTOGRAMS
  std::map<TString, TH1 *> allPlots;
  std::map<TString, TH2 *> all2dPlots;
  allPlots["puwgtctr"] = new TH1F("puwgtctr","Weight sums",2,0,2);

  std::vector<TString> lfsVec = { "E", "EE", "EM", "MM", "M" }; 
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      allPlots["nvtx_"+tag]  = new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",30,0,30);
      allPlots["njets_"+tag]  = new TH1F("njets_"+tag,";Jet multiplicity;Events",5,0,5);
      for(int i=0; i<2; i++)
        {
          TString pf(Form("l%d",i));          
          allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Lepton p_{t} [GeV];Events",50,0,250);
          allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Lepton pseudo-rapidity;Events",50,0,2.5);
        }
      for(int i=0; i<6; i++)
        {
          TString pf(Form("j%d",i));
          allPlots[pf+"pt_"+tag]  = new TH1F(pf+"pt_"+tag,";Jet transverse momentum [GeV];Events",50,0,250);
          allPlots[pf+"eta_"+tag]  = new TH1F(pf+"eta_"+tag,";Jet pseudo-rapidity;Events",50,0,4.7);
        }
    }
  for (auto& it : allPlots)   { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : all2dPlots) { it.second->Sumw2(); it.second->SetDirectory(0); }

  //nentries = 10000;
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopJetShapeEvent(tjsev);
      if(iev%100==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      
      //GENERATOR LEVEL SELECTION
      
      int sel_ngleptons = 0;
      int vet_ngleptons = 0;
      int sel_ngbjets   = 0;
      int sel_ngwcand   = 0;
      std::vector<Jet> genJets;
      
      //loop over gen jets and leptons from pseudotop producer
      //leptons
      for (int i = 0; i < ev.ng; i++) {
        if (abs(ev.g_id[i])>10) {
          //store to tree
          tjsev.gl_id [tjsev.ngl] = ev.g_id [i];
          tjsev.gl_pt [tjsev.ngl] = ev.g_pt [i];
          tjsev.gl_eta[tjsev.ngl] = ev.g_eta[i];
          tjsev.gl_phi[tjsev.ngl] = ev.g_phi[i];
          tjsev.gl_m  [tjsev.ngl] = ev.g_m  [i];
          tjsev.gl_id [tjsev.ngl] = ev.g_id [i];
          //count
          tjsev.ngl++;
          if      (ev.g_pt[i]>30 && abs(ev.g_eta[i])<2.1) sel_ngleptons++;
          else if (ev.g_pt[i]>15 && abs(ev.g_eta[i])<2.5) vet_ngleptons++;
        }
      }
      //jets
      for (int i = 0; i < ev.ng; i++) {
        if (ev.g_pt[i]>30 && abs(ev.g_eta[i])<2.4 && abs(ev.g_id[i])<10) {
          TLorentzVector jp4; jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
          //check overlap with leptons
          bool lepJetOverlap = false;
          for (int l = 0; l < tjsev.ngl; l++) {
            if (tjsev.gl_pt[l] < 30) continue;
            TLorentzVector lp4; lp4.SetPtEtaPhiM(tjsev.gl_pt[l],tjsev.gl_eta[l],tjsev.gl_phi[l],tjsev.gl_m[l]);
            if (jp4.DeltaR(lp4) < 0.4) lepJetOverlap = true;
          }
          if (lepJetOverlap) continue;
          genJets.push_back(Jet(jp4, abs(ev.g_id[i]), i));
          //count
          tjsev.ngj++;
          if (abs(ev.g_id[i])==5) sel_ngbjets++;
        }
      }
      
      //flag non-b jets as part of W boson candidates: flavor 0->1
      //TODO: matched events at 8 TeV: mu=84.23, sigma=12.39 GeV. Correction needed?
      for (int i = 0; i < tjsev.ngj; i++) {
        for (int j = i+1; j < tjsev.ngj; j++) {
          if (genJets[i].p4.DeltaR(genJets[j].p4) < 0.8) {
            genJets[i].overlap = 1;
            genJets[j].overlap = 1;
          }
          if (genJets[i].flavor==5 or genJets[j].flavor==5) continue;
          TLorentzVector wCand = genJets[i].p4 + genJets[j].p4;
          if (abs(wCand.M()-80.4) < 15.) {
            genJets[i].flavor = 1;
            genJets[j].flavor = 1;
            sel_ngwcand++;
          }
        }
      }
      
      //event selected on generator level?
      if (sel_ngbjets==2 && sel_ngwcand>0 && 
          sel_ngleptons==1 && vet_ngleptons==0) tjsev.gen_sel = 1;
      

      //RECO LEVEL SELECTION
      
      int sel_nbjets = 0;
      int sel_nwcand = 0;
      std::vector<Jet> jets;
      
      //account for pu weights and effect on normalization
      float puWeight(1.0);
      if(!ev.isData) 
        {
          puWeight=puWgtGr[0]->Eval(ev.putrue);  
          allPlots["puwgtctr"]->Fill(0.,1.0);
          allPlots["puwgtctr"]->Fill(1.,puWeight);
        }

      //select 1 good lepton
      std::vector<int> tightLeptons,looseLeptons;
      for(int il=0; il<ev.nl; il++)
        {
          bool passTightKin(ev.l_pt[il]>30 && fabs(ev.l_eta[il])<2.1);
          bool passLooseKin(ev.l_pt[il]>15 && fabs(ev.l_eta[il])<2.5);
          bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
          float relIso(ev.l_relIso[il]);
          bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1);
          bool passLooseIso( ev.l_id[il]==13 ? relIso<0.25 : (ev.l_pid[il]   )&0x1);
          if(passTightKin && passTightId && passTightIso) tightLeptons.push_back(il);
          else if(passLooseKin && passLooseIso)           looseLeptons.push_back(il);
        }
      
      //check if triggers have fired
      //bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      //bool hasEleTrigger((ev.elTrigger & 0x1)!=0);

      //decide the channel
      std::vector<int> selLeptons;
      selLeptons.insert(selLeptons.end(), tightLeptons.begin(), tightLeptons.end());
      selLeptons.insert(selLeptons.end(), looseLeptons.begin(), looseLeptons.end());
      /*
      if(tightLeptons.size()==1 && looseLeptons.size()==0)
        {
          selLeptons.push_back( tightLeptons[0] );
        }
      if(tightLeptons.size()==0 && looseLeptons.size()==1)
        {
          selLeptons.push_back( looseLeptons[0] );
        }
      if(tightLeptons.size()>=2)
        {
          selLeptons.push_back(tightLeptons[0]);
          selLeptons.push_back(tightLeptons[1]);
        }
      if(tightLeptons.size()==1 && looseLeptons.size()>=1)
        {
          selLeptons.push_back(tightLeptons[0]);
          selLeptons.push_back(looseLeptons[0]);          
        }
      */

      //save lepton kinematics
      std::vector<TLorentzVector> leptons;
      for(size_t il=0; il<selLeptons.size(); il++)
        {
          int lepIdx=selLeptons[il];
          TLorentzVector lp4;
          lp4.SetPtEtaPhiM(ev.l_pt[lepIdx],ev.l_eta[lepIdx],ev.l_phi[lepIdx],ev.l_mass[lepIdx]);
          leptons.push_back(lp4);
        }

      //select jets
      for (int k=0; k<ev.nj;k++)
        {
          //check kinematics
          TLorentzVector jp4;
          jp4.SetPtEtaPhiM(ev.j_pt[k],ev.j_eta[k],ev.j_phi[k],ev.j_mass[k]);

          //cross clean with leptons
          bool overlapsWithLepton(false);
          for(size_t il=0; il<leptons.size(); il++)
            {
              if(jp4.DeltaR(leptons[il])>0.4) continue;
              overlapsWithLepton=true;
            }
          if(overlapsWithLepton) continue;

          //smear jet energy resolution for MC
          float genJet_pt(0);
          if(ev.j_g[k]>-1) genJet_pt=ev.g_pt[ ev.j_g[k] ];
          if(!ev.isData && genJet_pt>0) 
            {
              float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Eta(),genJet_pt)[0];
              jp4 *= jerSmear;
            }
          
          //jet kinematic selection
          if(jp4.Pt() < 30 || abs(jp4.Eta()) > 2.4) continue;
          
          //b-tag
          float csv = ev.j_csv[k];          
          bool isBTagged(csv>0.800);
          if(!ev.isData)
            {
              float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
              float expEff(1.0), jetBtagSF(1.0);
              if(abs(ev.j_hadflav[k])==4) 
                {         
                  expEff    = expBtagEff["c"]->Eval(jptForBtag); 
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_C, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? expBtagEffPy8["c"]->Eval(jptForBtag)/expBtagEff["c"]->Eval(jptForBtag) : 0.;
                }
              else if(abs(ev.j_hadflav[k])==5) 
                { 
                  expEff    = expBtagEff["b"]->Eval(jptForBtag); 
                  jetBtagSF = sfbReaders[0]->eval( BTagEntry::FLAV_B, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff>0 ? expBtagEffPy8["b"]->Eval(jptForBtag)/expBtagEff["b"]->Eval(jptForBtag) : 0.;
                }
              else
                {
                  expEff    = expBtagEff["udsg"]->Eval(jptForBtag);
                  jetBtagSF = sflReaders[0]->eval( BTagEntry::FLAV_UDSG, jetaForBtag, jptForBtag);
                  jetBtagSF *= expEff> 0 ? expBtagEffPy8["udsg"]->Eval(jptForBtag)/expBtagEff["udsg"]->Eval(jptForBtag) : 0.;
                }
              
              //updated b-tagging decision with the data/MC scale factor
              myBTagSFUtil.modifyBTagsWithSF(isBTagged,    jetBtagSF,     expEff);
            }

          //flavor based on b tagging
          int flavor = 0;
          if (isBTagged) {
            flavor = 5;
            sel_nbjets++;
          }
          
          //store to jet analysis objects
          jets.push_back(Jet(jp4, flavor, k));
        }
      
      tjsev.nj=jets.size();
      
      //flag non-b jets as part of W boson candidates: flavor 0->1
      for (int i = 0; i < tjsev.nj; i++) {
        for (int j = i+1; j < tjsev.nj; j++) {
          if (jets[i].p4.DeltaR(jets[j].p4) < 0.8) {
            jets[i].overlap = 1;
            jets[j].overlap = 1;
          }
          if (jets[i].flavor==5 or jets[j].flavor==5) continue;
          TLorentzVector wCand = jets[i].p4 + jets[j].p4;
          if (abs(wCand.M()-80.4) < 15.) {
            jets[i].flavor = 1;
            jets[j].flavor = 1;
            sel_nwcand++;
          }
        }
      }
      
      //event selected on reco level?
      if (sel_nbjets==2 && sel_nwcand>0 && 
          tightLeptons.size()==1 && looseLeptons.size()==0) tjsev.reco_sel = 1;
      
      //proceed only if event is selected on gen or reco level
      if (tjsev.gen_sel + tjsev.reco_sel == -2) continue;
      
      
      //GEN LEVEL ANALYSIS
      
      //fill jet constituents
      for (int i = 0; i < tjsev.ngj; i++) {
        int idx = genJets[i].oldidx;
        for (int p = 0; p < ev.ngpf; p++) {
          if (ev.gpf_g[p] == idx) {
            TLorentzVector pp4;
            pp4.SetPtEtaPhiM(ev.gpf_pt[p],ev.gpf_eta[p],ev.gpf_phi[p],ev.gpf_m[p]);
            genJets[i].particles.push_back(Particle(pp4, ev.gpf_c[p], 1.));
          }
        }
      }
      
      //store jets to tree, including generalized angularities (jet shapes)
      for (int i = 0; i < tjsev.ngj; i++) {
        tjsev.gj_pt  [i] = genJets[i].p4.Pt();
        tjsev.gj_eta [i] = genJets[i].p4.Eta();
        tjsev.gj_phi [i] = genJets[i].p4.Phi();
        tjsev.gj_m   [i] = genJets[i].p4.M();
        tjsev.gj_flavor[i] = genJets[i].flavor;
        tjsev.gj_overlap[i] = genJets[i].overlap;
        //for (int i1 = 0; i1 < 3; ++i1) for (int i2 = 0; i2 < 3; ++i2) for (int i3 = 0; i3 < 3; ++i3) for (int i4 = 0; i4 < 3; ++i4)
        //tjsev.gj_ga[i][i1][i2][i3][i4] = calcGA(genJets[i], i1, i2, i3, i4);
        
        tjsev.gj_mult_charged[i] = getMult(genJets[i]);
        tjsev.gj_mult_all[i]     = getMult(genJets[i], true);
        tjsev.gj_mult_puppi[i]   = getMult(genJets[i], true, true);
        
        tjsev.gj_ptd_charged[i] = getPtD(genJets[i]);
        tjsev.gj_ptd_all[i]     = getPtD(genJets[i], true);
        tjsev.gj_ptd_puppi[i]   = getPtD(genJets[i], true, true);
        
        tjsev.gj_width_charged[i] = getWidth(genJets[i]);
        tjsev.gj_width_all[i]     = getWidth(genJets[i], true);
        tjsev.gj_width_puppi[i]   = getWidth(genJets[i], true, true);
        
        tjsev.gj_ecc_charged[i] = getEcc(genJets[i]);
        tjsev.gj_ecc_all[i]     = getEcc(genJets[i], true);
        tjsev.gj_ecc_puppi[i]   = getEcc(genJets[i], true, true);
        
        tjsev.gj_tau21_charged[i] = getTau(2, 1, genJets[i]);
        tjsev.gj_tau21_all[i]     = getTau(2, 1, genJets[i], true);
        tjsev.gj_tau21_puppi[i]   = getTau(2, 1, genJets[i], true, true);
        
        tjsev.gj_tau32_charged[i] = getTau(3, 2, genJets[i]);
        tjsev.gj_tau32_all[i]     = getTau(3, 2, genJets[i], true);
        tjsev.gj_tau32_puppi[i]   = getTau(3, 2, genJets[i], true, true);
        
        tjsev.gj_tau43_charged[i] = getTau(4, 3, genJets[i]);
        tjsev.gj_tau43_all[i]     = getTau(4, 3, genJets[i], true);
        tjsev.gj_tau43_puppi[i]   = getTau(4, 3, genJets[i], true, true);
        
        std::vector<double> zgResult_charged = getZg(genJets[i]);
        std::vector<double> zgResult_all     = getZg(genJets[i], true);
        std::vector<double> zgResult_puppi   = getZg(genJets[i], true, true);
        
        tjsev.gj_zg_charged[i]    = zgResult_charged[0];
        tjsev.gj_zg_all[i]        = zgResult_all[0];
        tjsev.gj_zg_puppi[i]      = zgResult_puppi[0];
        
        tjsev.gj_zgxdr_charged[i] = zgResult_charged[1];
        tjsev.gj_zgxdr_all[i]     = zgResult_all[1];
        tjsev.gj_zgxdr_puppi[i]   = zgResult_puppi[1];
        
        tjsev.gj_zgdr_charged[i]  = zgResult_charged[2];
        tjsev.gj_zgdr_all[i]      = zgResult_all[2];
        tjsev.gj_zgdr_puppi[i]    = zgResult_puppi[2];
        
        tjsev.gj_ga_width_charged[i]  = calcGA(1., 1., genJets[i]);
        tjsev.gj_ga_width_all[i]      = calcGA(1., 1., genJets[i], true);
        tjsev.gj_ga_width_puppi[i]    = calcGA(1., 1., genJets[i], true, true);
        
        tjsev.gj_ga_lha_charged[i]    = calcGA(0.5, 1., genJets[i]);
        tjsev.gj_ga_lha_all[i]        = calcGA(0.5, 1., genJets[i], true);
        tjsev.gj_ga_lha_puppi[i]      = calcGA(0.5, 1., genJets[i], true, true);
        
        tjsev.gj_ga_thrust_charged[i] = calcGA(2., 1., genJets[i]);
        tjsev.gj_ga_thrust_all[i]     = calcGA(2., 1., genJets[i], true);
        tjsev.gj_ga_thrust_puppi[i]   = calcGA(2., 1., genJets[i], true, true);
        
        tjsev.gj_c1_02_charged[i] = getC(1, 0.2, genJets[i]);
        tjsev.gj_c1_02_all[i]     = getC(1, 0.2, genJets[i], true);
        tjsev.gj_c1_02_puppi[i]   = getC(1, 0.2, genJets[i], true, true);
        tjsev.gj_c1_05_charged[i] = getC(1, 0.5, genJets[i]);
        tjsev.gj_c1_05_all[i]     = getC(1, 0.5, genJets[i], true);
        tjsev.gj_c1_05_puppi[i]   = getC(1, 0.5, genJets[i], true, true);
        tjsev.gj_c1_10_charged[i] = getC(1, 1.0, genJets[i]);
        tjsev.gj_c1_10_all[i]     = getC(1, 1.0, genJets[i], true);
        tjsev.gj_c1_10_puppi[i]   = getC(1, 1.0, genJets[i], true, true);
        tjsev.gj_c1_20_charged[i] = getC(1, 2.0, genJets[i]);
        tjsev.gj_c1_20_all[i]     = getC(1, 2.0, genJets[i], true);
        tjsev.gj_c1_20_puppi[i]   = getC(1, 2.0, genJets[i], true, true);
        tjsev.gj_c2_02_charged[i] = getC(2, 0.2, genJets[i]);
        tjsev.gj_c2_02_all[i]     = getC(2, 0.2, genJets[i], true);
        tjsev.gj_c2_02_puppi[i]   = getC(2, 0.2, genJets[i], true, true);
        tjsev.gj_c2_05_charged[i] = getC(2, 0.5, genJets[i]);
        tjsev.gj_c2_05_all[i]     = getC(2, 0.5, genJets[i], true);
        tjsev.gj_c2_05_puppi[i]   = getC(2, 0.5, genJets[i], true, true);
        tjsev.gj_c2_10_charged[i] = getC(2, 1.0, genJets[i]);
        tjsev.gj_c2_10_all[i]     = getC(2, 1.0, genJets[i], true);
        tjsev.gj_c2_10_puppi[i]   = getC(2, 1.0, genJets[i], true, true);
        tjsev.gj_c2_20_charged[i] = getC(2, 2.0, genJets[i]);
        tjsev.gj_c2_20_all[i]     = getC(2, 2.0, genJets[i], true);
        tjsev.gj_c2_20_puppi[i]   = getC(2, 2.0, genJets[i], true, true);
        tjsev.gj_c3_02_charged[i] = getC(3, 0.2, genJets[i]);
        tjsev.gj_c3_02_all[i]     = getC(3, 0.2, genJets[i], true);
        tjsev.gj_c3_02_puppi[i]   = getC(3, 0.2, genJets[i], true, true);
        tjsev.gj_c3_05_charged[i] = getC(3, 0.5, genJets[i]);
        tjsev.gj_c3_05_all[i]     = getC(3, 0.5, genJets[i], true);
        tjsev.gj_c3_05_puppi[i]   = getC(3, 0.5, genJets[i], true, true);
        tjsev.gj_c3_10_charged[i] = getC(3, 1.0, genJets[i]);
        tjsev.gj_c3_10_all[i]     = getC(3, 1.0, genJets[i], true);
        tjsev.gj_c3_10_puppi[i]   = getC(3, 1.0, genJets[i], true, true);
        tjsev.gj_c3_20_charged[i] = getC(3, 2.0, genJets[i]);
        tjsev.gj_c3_20_all[i]     = getC(3, 2.0, genJets[i], true);
        tjsev.gj_c3_20_puppi[i]   = getC(3, 2.0, genJets[i], true, true);

      }

      
      //RECO LEVEL ANALYSIS
      
      //event weight
      float wgt(1.0);
      if(!ev.isData)
        {
          //update lepton selection scale factors, if found
          float lepTriggerSF(1.0),lepSelSF(1.0);
          for(UInt_t il=0; il<leptons.size(); il++)
            {
              TString prefix(abs(ev.l_id[il])==11 ? "e" : "m");
              float minEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmin() ), maxEtaForEff( lepEffH[prefix+"_sel"]->GetXaxis()->GetXmax()-0.01 );
              float etaForEff=TMath::Max(TMath::Min(float(fabs(leptons[il].Eta())),maxEtaForEff),minEtaForEff);
              Int_t etaBinForEff=lepEffH[prefix+"_sel"]->GetXaxis()->FindBin(etaForEff);
              
              float minPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmin() ), maxPtForEff( lepEffH[prefix+"_sel"]->GetYaxis()->GetXmax()-0.01 );
              float ptForEff=TMath::Max(TMath::Min(float(leptons[il].Pt()),maxPtForEff),minPtForEff);
              Int_t ptBinForEff=lepEffH[prefix+"_sel"]->GetYaxis()->FindBin(ptForEff);
                                    
              lepSelSF=(lepEffH[prefix+"_sel"]->GetBinContent(etaBinForEff,ptBinForEff));
              
              if(il!=0) continue;
              if(prefix=="m")
                {
                  lepTriggerSF=(lepEffH[prefix+"_trig"]->GetBinContent(etaBinForEff,ptBinForEff));
                
                }
            }

          //update nominal event weight
          float norm( normH ? normH->GetBinContent(1) : 1.0);
          wgt=lepTriggerSF*lepSelSF*puWeight*norm;
          if(ev.ttbar_nw>0) wgt*=ev.ttbar_w[0];
        }
      
      //nominal selection control histograms
      //TODO allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);

      //fill leptons
      tjsev.nl=leptons.size();
      for(int il=0; il<(int)leptons.size(); il++)
        {
          tjsev.l_pt[il]=leptons[il].Pt();
          tjsev.l_eta[il]=leptons[il].Eta();
          tjsev.l_phi[il]=leptons[il].Phi();
          tjsev.l_m[il]=leptons[il].M();
          tjsev.l_id[il]=ev.l_id[ selLeptons[il] ];
        }
      
      //fill jet constituents
      for (int i = 0; i < tjsev.nj; i++) {
        int idx = jets[i].oldidx;
        for (int p = 0; p < ev.npf; p++) {
          if (ev.pf_j[p] == idx) {
            TLorentzVector pp4;
            pp4.SetPtEtaPhiM(ev.pf_pt[p],ev.pf_eta[p],ev.pf_phi[p],ev.pf_m[p]);
            jets[i].particles.push_back(Particle(pp4, ev.pf_c[p], ev.pf_puppiWgt[p]));
          }
        }
      }
      
      //fill jets (with jet shapes)
      for(int ij=0; ij<(int)jets.size(); ij++)
        {
          tjsev.j_pt[ij]   = jets[ij].p4.Pt();
          tjsev.j_eta[ij]  = jets[ij].p4.Eta();
          tjsev.j_phi[ij]  = jets[ij].p4.Phi();
          tjsev.j_m[ij]    = jets[ij].p4.M(); 
          tjsev.j_flavor[ij] = jets[ij].flavor;
          tjsev.j_overlap[ij] = jets[ij].overlap;
          //for (int i1 = 0; i1 < 3; ++i1) for (int i2 = 0; i2 < 3; ++i2) for (int i3 = 0; i3 < 3; ++i3) for (int i4 = 0; i4 < 3; ++i4)
          //tjsev.j_ga[ij][i1][i2][i3][i4] = calcGA(jets[ij], i1, i2, i3, i4);
          
          tjsev.j_mult_charged[ij] = getMult(jets[ij]);
          tjsev.j_mult_all[ij]     = getMult(jets[ij], true);
          tjsev.j_mult_puppi[ij]   = getMult(jets[ij], true, true);
          
          tjsev.j_ptd_charged[ij] = getPtD(jets[ij]);
          tjsev.j_ptd_all[ij]     = getPtD(jets[ij], true);
          tjsev.j_ptd_puppi[ij]   = getPtD(jets[ij], true, true);
          
          tjsev.j_width_charged[ij] = getWidth(jets[ij]);
          tjsev.j_width_all[ij]     = getWidth(jets[ij], true);
          tjsev.j_width_puppi[ij]   = getWidth(jets[ij], true, true);
          
          tjsev.j_ecc_charged[ij] = getEcc(jets[ij]);
          tjsev.j_ecc_all[ij]     = getEcc(jets[ij], true);
          tjsev.j_ecc_puppi[ij]   = getEcc(jets[ij], true, true);
          
          tjsev.j_tau21_charged[ij] = getTau(2, 1, jets[ij]);
          tjsev.j_tau21_all[ij]     = getTau(2, 1, jets[ij], true);
          tjsev.j_tau21_puppi[ij]   = getTau(2, 1, jets[ij], true, true);
          
          tjsev.j_tau32_charged[ij] = getTau(3, 2, jets[ij]);
          tjsev.j_tau32_all[ij]     = getTau(3, 2, jets[ij], true);
          tjsev.j_tau32_puppi[ij]   = getTau(3, 2, jets[ij], true, true);
          
          tjsev.j_tau43_charged[ij] = getTau(4, 3, jets[ij]);
          tjsev.j_tau43_all[ij]     = getTau(4, 3, jets[ij], true);
          tjsev.j_tau43_puppi[ij]   = getTau(4, 3, jets[ij], true, true);
          
          std::vector<double> zgResult_charged = getZg(jets[ij]);
          std::vector<double> zgResult_all     = getZg(jets[ij], true);
          std::vector<double> zgResult_puppi   = getZg(jets[ij], true, true);
          
          tjsev.j_zg_charged[ij]    = zgResult_charged[0];
          tjsev.j_zg_all[ij]        = zgResult_all[0];
          tjsev.j_zg_puppi[ij]      = zgResult_puppi[0];
          
          tjsev.j_zgxdr_charged[ij] = zgResult_charged[1];
          tjsev.j_zgxdr_all[ij]     = zgResult_all[1];
          tjsev.j_zgxdr_puppi[ij]   = zgResult_puppi[1];
          
          tjsev.j_zgdr_charged[ij]  = zgResult_charged[2];
          tjsev.j_zgdr_all[ij]      = zgResult_all[2];
          tjsev.j_zgdr_puppi[ij]    = zgResult_puppi[2];
          
          tjsev.j_ga_width_charged[ij]  = calcGA(1., 1., jets[ij]);
          tjsev.j_ga_width_all[ij]      = calcGA(1., 1., jets[ij], true);
          tjsev.j_ga_width_puppi[ij]    = calcGA(1., 1., jets[ij], true, true);
          
          tjsev.j_ga_lha_charged[ij]    = calcGA(0.5, 1., jets[ij]);
          tjsev.j_ga_lha_all[ij]        = calcGA(0.5, 1., jets[ij], true);
          tjsev.j_ga_lha_puppi[ij]      = calcGA(0.5, 1., jets[ij], true, true);
          
          tjsev.j_ga_thrust_charged[ij] = calcGA(2., 1., jets[ij]);
          tjsev.j_ga_thrust_all[ij]     = calcGA(2., 1., jets[ij], true);
          tjsev.j_ga_thrust_puppi[ij]   = calcGA(2., 1., jets[ij], true, true);
          
          tjsev.j_c1_02_charged[ij] = getC(1, 0.2, jets[ij]);
          tjsev.j_c1_02_all[ij]     = getC(1, 0.2, jets[ij], true);
          tjsev.j_c1_02_puppi[ij]   = getC(1, 0.2, jets[ij], true, true);
          tjsev.j_c1_05_charged[ij] = getC(1, 0.5, jets[ij]);
          tjsev.j_c1_05_all[ij]     = getC(1, 0.5, jets[ij], true);
          tjsev.j_c1_05_puppi[ij]   = getC(1, 0.5, jets[ij], true, true);
          tjsev.j_c1_10_charged[ij] = getC(1, 1.0, jets[ij]);
          tjsev.j_c1_10_all[ij]     = getC(1, 1.0, jets[ij], true);
          tjsev.j_c1_10_puppi[ij]   = getC(1, 1.0, jets[ij], true, true);
          tjsev.j_c1_20_charged[ij] = getC(1, 2.0, jets[ij]);
          tjsev.j_c1_20_all[ij]     = getC(1, 2.0, jets[ij], true);
          tjsev.j_c1_20_puppi[ij]   = getC(1, 2.0, jets[ij], true, true);
          tjsev.j_c2_02_charged[ij] = getC(2, 0.2, jets[ij]);
          tjsev.j_c2_02_all[ij]     = getC(2, 0.2, jets[ij], true);
          tjsev.j_c2_02_puppi[ij]   = getC(2, 0.2, jets[ij], true, true);
          tjsev.j_c2_05_charged[ij] = getC(2, 0.5, jets[ij]);
          tjsev.j_c2_05_all[ij]     = getC(2, 0.5, jets[ij], true);
          tjsev.j_c2_05_puppi[ij]   = getC(2, 0.5, jets[ij], true, true);
          tjsev.j_c2_10_charged[ij] = getC(2, 1.0, jets[ij]);
          tjsev.j_c2_10_all[ij]     = getC(2, 1.0, jets[ij], true);
          tjsev.j_c2_10_puppi[ij]   = getC(2, 1.0, jets[ij], true, true);
          tjsev.j_c2_20_charged[ij] = getC(2, 2.0, jets[ij]);
          tjsev.j_c2_20_all[ij]     = getC(2, 2.0, jets[ij], true);
          tjsev.j_c2_20_puppi[ij]   = getC(2, 2.0, jets[ij], true, true);
          tjsev.j_c3_02_charged[ij] = getC(3, 0.2, jets[ij]);
          tjsev.j_c3_02_all[ij]     = getC(3, 0.2, jets[ij], true);
          tjsev.j_c3_02_puppi[ij]   = getC(3, 0.2, jets[ij], true, true);
          tjsev.j_c3_05_charged[ij] = getC(3, 0.5, jets[ij]);
          tjsev.j_c3_05_all[ij]     = getC(3, 0.5, jets[ij], true);
          tjsev.j_c3_05_puppi[ij]   = getC(3, 0.5, jets[ij], true, true);
          tjsev.j_c3_10_charged[ij] = getC(3, 1.0, jets[ij]);
          tjsev.j_c3_10_all[ij]     = getC(3, 1.0, jets[ij], true);
          tjsev.j_c3_10_puppi[ij]   = getC(3, 1.0, jets[ij], true, true);
          tjsev.j_c3_20_charged[ij] = getC(3, 2.0, jets[ij]);
          tjsev.j_c3_20_all[ij]     = getC(3, 2.0, jets[ij], true);
          tjsev.j_c3_20_puppi[ij]   = getC(3, 2.0, jets[ij], true, true);

          //matching to gen jet
          for(int ig=0; ig<(int)genJets.size(); ig++) {
	          if(jets[ij].p4.DeltaR(genJets[ig].p4)>0.4) continue;
	          tjsev.j_gj[ij] = ig;
	          tjsev.gj_j[ig] = ij;
	          break;
	        }
        }

      tjsev.nw=1;
      tjsev.weight[0]=wgt;
      tjsev.met_pt=ev.met_pt[0];
      tjsev.met_phi=ev.met_phi[0];
      
      outT->Fill();
    }
  
  //close input file
  f->Close();

  //save histos to file  
  fOut->cd();
  outT->Write();
  for (auto& it : allPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : all2dPlots)  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  fOut->Close();
}

//Calculate generalized angularities
/*
  Definition from arXiv:1408.3122
  (beta,kappa) =
    (0,0) -> multiplicity
    (0,2) -> pt dispersion
    (1,1) -> broadening/width/girth
    (2,1) -> thrust (m^2/e)
*/
double calcGA(double beta, double kappa, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  int mult = 0;
  double sumpt = 0.;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    sumpt  += p.pt()*weight;
  }
  if (mult < 1) return -1.;
  //std::cout << "sumpt" << beta << kappa << iptcut << icharge << ": " << sumpt << std::endl;
  
  double ga = 0.;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    ga += weight * pow(p.p4.Pt()/sumpt, kappa) * pow(jet.p4.DeltaR(p.p4)/0.4, beta);
  }
  
  //std::cout << "ga" << beta << kappa << iptcut << icharge << ": " << ga << std::endl;
  
  /*
  if (ga > 1. && beta == 1 && kappa == 1 && iptcut == 0 && icharge == 0) {
    std::cout << "ga(1,1) = " << ga << std::endl;
    std::cout << "pt/sum(pt) deltaR" << std::endl;
    double sumpt_test = 0;
    for (auto p : jet.particles) {
      if (p.p4.Pt() < (iptcut+1)*0.500) continue;
      if(p.charge==0) continue;
      std::cout << p.p4.Pt()/sumpt << " " << jet.p4.DeltaR(p.p4) << std::endl;
      sumpt_test += p.p4.Pt()/sumpt;
    }
    std::cout << "sumpt_test = " << sumpt_test << "\n" << std::endl;
  }
  */
  
  return ga;
}

double getMult(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getMult()" << std::endl;
  double sumWeight = 0.;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    sumWeight += weight;
  }
  return sumWeight;
}

double getPtD(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getPtD()" << std::endl;
  int mult = 0;
  double sumpt  = 0.;
  double sumpt2 = 0.;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    sumpt  += p.pt()*weight;
    sumpt2 += pow(p.pt()*weight, 2);
  }
  if (mult < 2) return -1.;
  double ptd = sqrt(sumpt2)/sumpt;
  return ptd;
}

//TODO: use version of jet with particle cuts applied, and recalculated 4vec?
double getWidth(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getWidth()" << std::endl;
  int mult = 0;
  double sumpt   = 0.;
  double sumptdr = 0.;
  Particle axis(TLorentzVector(), 0., 0.);
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    axis.p4 += p.momentum()*weight;
  }
  if (mult < 2) return -1.;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    sumpt   += p.pt()*weight;
    sumptdr += p.pt()*weight * deltaR(axis.momentum(), p.momentum());
  }
  double width = sumptdr/sumpt;
  return width;
}

double getEcc(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getEcc()" << std::endl;
  // Get mean axis
  int mult = 0;
  Particle axis(TLorentzVector(), 0., 0.);
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    axis.p4 += p.momentum()*weight;
  }
  if (mult < 4) return -1.;
  // Covariance matrix
  Matrix<2> M;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    Matrix<2> MPart;
    MPart.set(0, 0, (p.eta() - axis.eta()) * (p.eta() - axis.eta()));
    MPart.set(0, 1, (p.eta() - axis.eta()) * mapAngleMPiToPi(p.phi() - axis.phi()));
    MPart.set(1, 0, mapAngleMPiToPi(p.phi() - axis.phi()) * (p.eta() - axis.eta()));
    MPart.set(1, 1, mapAngleMPiToPi(p.phi() - axis.phi()) * mapAngleMPiToPi(p.phi() - axis.phi()));
    M += MPart * p.energy() * weight;
  }
  // Calculate eccentricity from eigenvalues
  const EigenSystem<2> eigen = diagonalize(M);
  double ecc = 1. - eigen.getEigenValues()[1]/eigen.getEigenValues()[0];
  
  return ecc;
}

double getTau(int N, int M, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < N+1) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  NsubjettinessRatio tau_ratio(N, M, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.0));
  
  return tau_ratio(cajet);
}

double getC(int N, double beta, Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < N+1) return -1.;
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  
  EnergyCorrelatorDoubleRatio C(N, beta);
  
  return C(cajet);
}

std::vector<double> getZg(Jet jet, bool includeNeutrals, bool usePuppi, double ptcut) {
  //std::cout << "getZg()" << std::endl;
  // Recluster constituents with CA
  int mult = 0.;
  vector<PseudoJet> particles;
  for (auto p : jet.particles) {
    if (not includeNeutrals and p.charge == 0) continue;
    if (ptcut > p.pt()) continue;
    double weight = usePuppi ? p.puppi : 1.;
    if (weight > 0.) ++mult;
    particles.push_back( PseudoJet(p.px(), p.py(), p.pz(), p.e())*weight );
  }
  if (mult < 2) return {-1., -1.};
  
  JetDefinition jet_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
  
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(1));
  
  PseudoJet cajet = jets[0];
  PseudoJet cajetp1, cajetp2;
  double zg = 0.;
  while (zg < 0.1 and cajet.has_parents(cajetp1, cajetp2)) {
    zg    = cajetp2.pt()/cajet.pt();
    cajet = cajetp1;
  }
  //std::cout << "zg = " << zg << std::endl;
  std::vector<double> results;
  results.push_back(zg);
  results.push_back(zg*cajetp1.delta_R(cajetp2));
  results.push_back(cajetp1.delta_R(cajetp2));
  return results;
}

//
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev)
{
  //event weights
  t->Branch("nw",  &tjsev.nw, "nw/I");
  t->Branch("weight",  tjsev.weight, "weight[nw]/F");

  //met
  t->Branch("met_pt",  &tjsev.met_pt, "met_pt/F");
  t->Branch("met_phi",  &tjsev.met_phi, "met_phi/F");

  //leptons
  t->Branch("nl",  &tjsev.nl, "nl/I");
  t->Branch("l_pt",  tjsev.l_pt ,  "l_pt[nl]/F");
  t->Branch("l_eta", tjsev.l_eta , "l_eta[nl]/F");
  t->Branch("l_phi", tjsev.l_phi , "l_phi[nl]/F");
  t->Branch("l_m",   tjsev.l_m ,   "l_m[nl]/F");
  t->Branch("l_id",   tjsev.l_id ,   "l_id[nl]/I");
  
  //gen leptons
  t->Branch("ngl",  &tjsev.ngl, "ngl/I");
  t->Branch("gl_pt",  tjsev.gl_pt ,  "gl_pt[ngl]/F");
  t->Branch("gl_eta", tjsev.gl_eta , "gl_eta[ngl]/F");
  t->Branch("gl_phi", tjsev.gl_phi , "gl_phi[ngl]/F");
  t->Branch("gl_m",   tjsev.gl_m ,   "gl_m[ngl]/F");
  t->Branch("gl_id",  tjsev.gl_id ,  "gl_id[ngl]/I");

  //jets
  t->Branch("nj",  &tjsev.nj, "nj/I");
  t->Branch("j_pt",  tjsev.j_pt ,  "j_pt[nj]/F");
  t->Branch("j_eta", tjsev.j_eta , "j_eta[nj]/F");
  t->Branch("j_phi", tjsev.j_phi , "j_phi[nj]/F");
  t->Branch("j_m",   tjsev.j_m ,   "j_m[nj]/F");
  t->Branch("j_flavor",  tjsev.j_flavor ,  "j_flavor[nj]/I");
  t->Branch("j_overlap",  tjsev.j_overlap ,  "j_overlap[nj]/I");
  t->Branch("j_gj",  tjsev.j_gj ,  "j_gj[nj]/I");
  
  //gen jets
  t->Branch("ngj",  &tjsev.ngj, "ngj/I");
  t->Branch("gj_pt",  tjsev.gj_pt ,  "gj_pt[ngj]/F");
  t->Branch("gj_eta", tjsev.gj_eta , "gj_eta[ngj]/F");
  t->Branch("gj_phi", tjsev.gj_phi , "gj_phi[ngj]/F");
  t->Branch("gj_m",   tjsev.gj_m ,   "gj_m[ngj]/F");
  t->Branch("gj_flavor",  tjsev.gj_flavor ,  "gj_flavor[ngj]/I");
  t->Branch("gj_overlap",  tjsev.gj_overlap ,  "gj_overlap[ngj]/I");
  t->Branch("gj_j",  tjsev.gj_j ,  "gj_j[ngj]/I");
  
  t->Branch("gen_sel", &tjsev.gen_sel ,  "gen_sel/I");
  t->Branch("reco_sel", &tjsev.reco_sel ,  "reco_sel/I");
  
  //observables
  t->Branch("j_mult_charged",   tjsev.j_mult_charged,   "j_mult_charged[nj]/F");
  t->Branch("j_mult_puppi",   tjsev.j_mult_puppi,   "j_mult_puppi[nj]/F");
  t->Branch("j_mult_all",   tjsev.j_mult_all,   "j_mult_all[nj]/F");
  t->Branch("j_width_charged",   tjsev.j_width_charged,   "j_width_charged[nj]/F");
  t->Branch("j_width_puppi",   tjsev.j_width_puppi,   "j_width_puppi[nj]/F");
  t->Branch("j_width_all",   tjsev.j_width_all,   "j_width_all[nj]/F");
  t->Branch("j_ptd_charged",   tjsev.j_ptd_charged,   "j_ptd_charged[nj]/F");
  t->Branch("j_ptd_puppi",   tjsev.j_ptd_puppi,   "j_ptd_puppi[nj]/F");
  t->Branch("j_ptd_all",   tjsev.j_ptd_all,   "j_ptd_all[nj]/F");
  t->Branch("j_ecc_charged",   tjsev.j_ecc_charged,   "j_ecc_charged[nj]/F");
  t->Branch("j_ecc_puppi",   tjsev.j_ecc_puppi,   "j_ecc_puppi[nj]/F");
  t->Branch("j_ecc_all",   tjsev.j_ecc_all,   "j_ecc_all[nj]/F");
  t->Branch("j_tau21_charged",   tjsev.j_tau21_charged,   "j_tau21_charged[nj]/F");
  t->Branch("j_tau21_puppi",   tjsev.j_tau21_puppi,   "j_tau21_puppi[nj]/F");
  t->Branch("j_tau21_all",   tjsev.j_tau21_all,   "j_tau21_all[nj]/F");
  t->Branch("j_tau32_charged",   tjsev.j_tau32_charged,   "j_tau32_charged[nj]/F");
  t->Branch("j_tau32_puppi",   tjsev.j_tau32_puppi,   "j_tau32_puppi[nj]/F");
  t->Branch("j_tau32_all",   tjsev.j_tau32_all,   "j_tau32_all[nj]/F");
  t->Branch("j_tau43_charged",   tjsev.j_tau43_charged,   "j_tau43_charged[nj]/F");
  t->Branch("j_tau43_puppi",   tjsev.j_tau43_puppi,   "j_tau43_puppi[nj]/F");
  t->Branch("j_tau43_all",   tjsev.j_tau43_all,   "j_tau43_all[nj]/F");
  t->Branch("j_zg_charged",   tjsev.j_zg_charged,   "j_zg_charged[nj]/F");
  t->Branch("j_zg_puppi",   tjsev.j_zg_puppi,   "j_zg_puppi[nj]/F");
  t->Branch("j_zg_all",   tjsev.j_zg_all,   "j_zg_all[nj]/F");
  t->Branch("j_zgxdr_charged",   tjsev.j_zgxdr_charged,   "j_zgxdr_charged[nj]/F");
  t->Branch("j_zgxdr_puppi",   tjsev.j_zgxdr_puppi,   "j_zgxdr_puppi[nj]/F");
  t->Branch("j_zgxdr_all",   tjsev.j_zgxdr_all,   "j_zgxdr_all[nj]/F");
  t->Branch("j_zgdr_charged",   tjsev.j_zgdr_charged,   "j_zgdr_charged[nj]/F");
  t->Branch("j_zgdr_puppi",   tjsev.j_zgdr_puppi,   "j_zgdr_puppi[nj]/F");
  t->Branch("j_zgdr_all",   tjsev.j_zgdr_all,   "j_zgdr_all[nj]/F");
  t->Branch("j_ga_width_charged",   tjsev.j_ga_width_charged,   "j_ga_width_charged[nj]/F");
  t->Branch("j_ga_width_puppi",   tjsev.j_ga_width_puppi,   "j_ga_width_puppi[nj]/F");
  t->Branch("j_ga_width_all",   tjsev.j_ga_width_all,   "j_ga_width_all[nj]/F");
  t->Branch("j_ga_lha_charged",   tjsev.j_ga_lha_charged,   "j_ga_lha_charged[nj]/F");
  t->Branch("j_ga_lha_puppi",   tjsev.j_ga_lha_puppi,   "j_ga_lha_puppi[nj]/F");
  t->Branch("j_ga_lha_all",   tjsev.j_ga_lha_all,   "j_ga_lha_all[nj]/F");
  t->Branch("j_ga_thrust_charged",   tjsev.j_ga_thrust_charged,   "j_ga_thrust_charged[nj]/F");
  t->Branch("j_ga_thrust_puppi",   tjsev.j_ga_thrust_puppi,   "j_ga_thrust_puppi[nj]/F");
  t->Branch("j_ga_thrust_all",   tjsev.j_ga_thrust_all,   "j_ga_thrust_all[nj]/F");
  t->Branch("j_c1_02_charged",   tjsev.j_c1_02_charged,   "j_c1_02_charged[nj]/F");
  t->Branch("j_c1_02_puppi",   tjsev.j_c1_02_puppi,   "j_c1_02_puppi[nj]/F");
  t->Branch("j_c1_02_all",   tjsev.j_c1_02_all,   "j_c1_02_all[nj]/F");
  t->Branch("j_c1_05_charged",   tjsev.j_c1_05_charged,   "j_c1_05_charged[nj]/F");
  t->Branch("j_c1_05_puppi",   tjsev.j_c1_05_puppi,   "j_c1_05_puppi[nj]/F");
  t->Branch("j_c1_05_all",   tjsev.j_c1_05_all,   "j_c1_05_all[nj]/F");
  t->Branch("j_c1_10_charged",   tjsev.j_c1_10_charged,   "j_c1_10_charged[nj]/F");
  t->Branch("j_c1_10_puppi",   tjsev.j_c1_10_puppi,   "j_c1_10_puppi[nj]/F");
  t->Branch("j_c1_10_all",   tjsev.j_c1_10_all,   "j_c1_10_all[nj]/F");
  t->Branch("j_c1_20_charged",   tjsev.j_c1_20_charged,   "j_c1_20_charged[nj]/F");
  t->Branch("j_c1_20_puppi",   tjsev.j_c1_20_puppi,   "j_c1_20_puppi[nj]/F");
  t->Branch("j_c1_20_all",   tjsev.j_c1_20_all,   "j_c1_20_all[nj]/F");
  t->Branch("j_c2_02_charged",   tjsev.j_c2_02_charged,   "j_c2_02_charged[nj]/F");
  t->Branch("j_c2_02_puppi",   tjsev.j_c2_02_puppi,   "j_c2_02_puppi[nj]/F");
  t->Branch("j_c2_02_all",   tjsev.j_c2_02_all,   "j_c2_02_all[nj]/F");
  t->Branch("j_c2_05_charged",   tjsev.j_c2_05_charged,   "j_c2_05_charged[nj]/F");
  t->Branch("j_c2_05_puppi",   tjsev.j_c2_05_puppi,   "j_c2_05_puppi[nj]/F");
  t->Branch("j_c2_05_all",   tjsev.j_c2_05_all,   "j_c2_05_all[nj]/F");
  t->Branch("j_c2_10_charged",   tjsev.j_c2_10_charged,   "j_c2_10_charged[nj]/F");
  t->Branch("j_c2_10_puppi",   tjsev.j_c2_10_puppi,   "j_c2_10_puppi[nj]/F");
  t->Branch("j_c2_10_all",   tjsev.j_c2_10_all,   "j_c2_10_all[nj]/F");
  t->Branch("j_c2_20_charged",   tjsev.j_c2_20_charged,   "j_c2_20_charged[nj]/F");
  t->Branch("j_c2_20_puppi",   tjsev.j_c2_20_puppi,   "j_c2_20_puppi[nj]/F");
  t->Branch("j_c2_20_all",   tjsev.j_c2_20_all,   "j_c2_20_all[nj]/F");
  t->Branch("j_c3_02_charged",   tjsev.j_c3_02_charged,   "j_c3_02_charged[nj]/F");
  t->Branch("j_c3_02_puppi",   tjsev.j_c3_02_puppi,   "j_c3_02_puppi[nj]/F");
  t->Branch("j_c3_02_all",   tjsev.j_c3_02_all,   "j_c3_02_all[nj]/F");
  t->Branch("j_c3_05_charged",   tjsev.j_c3_05_charged,   "j_c3_05_charged[nj]/F");
  t->Branch("j_c3_05_puppi",   tjsev.j_c3_05_puppi,   "j_c3_05_puppi[nj]/F");
  t->Branch("j_c3_05_all",   tjsev.j_c3_05_all,   "j_c3_05_all[nj]/F");
  t->Branch("j_c3_10_charged",   tjsev.j_c3_10_charged,   "j_c3_10_charged[nj]/F");
  t->Branch("j_c3_10_puppi",   tjsev.j_c3_10_puppi,   "j_c3_10_puppi[nj]/F");
  t->Branch("j_c3_10_all",   tjsev.j_c3_10_all,   "j_c3_10_all[nj]/F");
  t->Branch("j_c3_20_charged",   tjsev.j_c3_20_charged,   "j_c3_20_charged[nj]/F");
  t->Branch("j_c3_20_puppi",   tjsev.j_c3_20_puppi,   "j_c3_20_puppi[nj]/F");
  t->Branch("j_c3_20_all",   tjsev.j_c3_20_all,   "j_c3_20_all[nj]/F");


  t->Branch("gj_mult_charged",   tjsev.gj_mult_charged,   "gj_mult_charged[ngj]/F");
  t->Branch("gj_mult_puppi",   tjsev.gj_mult_puppi,   "gj_mult_puppi[ngj]/F");
  t->Branch("gj_mult_all",   tjsev.gj_mult_all,   "gj_mult_all[ngj]/F");
  t->Branch("gj_width_charged",   tjsev.gj_width_charged,   "gj_width_charged[ngj]/F");
  t->Branch("gj_width_puppi",   tjsev.gj_width_puppi,   "gj_width_puppi[ngj]/F");
  t->Branch("gj_width_all",   tjsev.gj_width_all,   "gj_width_all[ngj]/F");
  t->Branch("gj_ptd_charged",   tjsev.gj_ptd_charged,   "gj_ptd_charged[ngj]/F");
  t->Branch("gj_ptd_puppi",   tjsev.gj_ptd_puppi,   "gj_ptd_puppi[ngj]/F");
  t->Branch("gj_ptd_all",   tjsev.gj_ptd_all,   "gj_ptd_all[ngj]/F");
  t->Branch("gj_ecc_charged",   tjsev.gj_ecc_charged,   "gj_ecc_charged[ngj]/F");
  t->Branch("gj_ecc_puppi",   tjsev.gj_ecc_puppi,   "gj_ecc_puppi[ngj]/F");
  t->Branch("gj_ecc_all",   tjsev.gj_ecc_all,   "gj_ecc_all[ngj]/F");
  t->Branch("gj_tau21_charged",   tjsev.gj_tau21_charged,   "gj_tau21_charged[ngj]/F");
  t->Branch("gj_tau21_puppi",   tjsev.gj_tau21_puppi,   "gj_tau21_puppi[ngj]/F");
  t->Branch("gj_tau21_all",   tjsev.gj_tau21_all,   "gj_tau21_all[ngj]/F");
  t->Branch("gj_tau32_charged",   tjsev.gj_tau32_charged,   "gj_tau32_charged[ngj]/F");
  t->Branch("gj_tau32_puppi",   tjsev.gj_tau32_puppi,   "gj_tau32_puppi[ngj]/F");
  t->Branch("gj_tau32_all",   tjsev.gj_tau32_all,   "gj_tau32_all[ngj]/F");
  t->Branch("gj_tau43_charged",   tjsev.gj_tau43_charged,   "gj_tau43_charged[ngj]/F");
  t->Branch("gj_tau43_puppi",   tjsev.gj_tau43_puppi,   "gj_tau43_puppi[ngj]/F");
  t->Branch("gj_tau43_all",   tjsev.gj_tau43_all,   "gj_tau43_all[ngj]/F");
  t->Branch("gj_zg_charged",   tjsev.gj_zg_charged,   "gj_zg_charged[ngj]/F");
  t->Branch("gj_zg_puppi",   tjsev.gj_zg_puppi,   "gj_zg_puppi[ngj]/F");
  t->Branch("gj_zg_all",   tjsev.gj_zg_all,   "gj_zg_all[ngj]/F");
  t->Branch("gj_zgxdr_charged",   tjsev.gj_zgxdr_charged,   "gj_zgxdr_charged[ngj]/F");
  t->Branch("gj_zgxdr_puppi",   tjsev.gj_zgxdr_puppi,   "gj_zgxdr_puppi[ngj]/F");
  t->Branch("gj_zgxdr_all",   tjsev.gj_zgxdr_all,   "gj_zgxdr_all[ngj]/F");
  t->Branch("gj_zgdr_charged",   tjsev.gj_zgdr_charged,   "gj_zgdr_charged[ngj]/F");
  t->Branch("gj_zgdr_puppi",   tjsev.gj_zgdr_puppi,   "gj_zgdr_puppi[ngj]/F");
  t->Branch("gj_zgdr_all",   tjsev.gj_zgdr_all,   "gj_zgdr_all[ngj]/F");
  t->Branch("gj_ga_width_charged",   tjsev.gj_ga_width_charged,   "gj_ga_width_charged[ngj]/F");
  t->Branch("gj_ga_width_puppi",   tjsev.gj_ga_width_puppi,   "gj_ga_width_puppi[ngj]/F");
  t->Branch("gj_ga_width_all",   tjsev.gj_ga_width_all,   "gj_ga_width_all[ngj]/F");
  t->Branch("gj_ga_lha_charged",   tjsev.gj_ga_lha_charged,   "gj_ga_lha_charged[ngj]/F");
  t->Branch("gj_ga_lha_puppi",   tjsev.gj_ga_lha_puppi,   "gj_ga_lha_puppi[ngj]/F");
  t->Branch("gj_ga_lha_all",   tjsev.gj_ga_lha_all,   "gj_ga_lha_all[ngj]/F");
  t->Branch("gj_ga_thrust_charged",   tjsev.gj_ga_thrust_charged,   "gj_ga_thrust_charged[ngj]/F");
  t->Branch("gj_ga_thrust_puppi",   tjsev.gj_ga_thrust_puppi,   "gj_ga_thrust_puppi[ngj]/F");
  t->Branch("gj_ga_thrust_all",   tjsev.gj_ga_thrust_all,   "gj_ga_thrust_all[ngj]/F");
  t->Branch("gj_c1_02_charged",   tjsev.gj_c1_02_charged,   "gj_c1_02_charged[ngj]/F");
  t->Branch("gj_c1_02_puppi",   tjsev.gj_c1_02_puppi,   "gj_c1_02_puppi[ngj]/F");
  t->Branch("gj_c1_02_all",   tjsev.gj_c1_02_all,   "gj_c1_02_all[ngj]/F");
  t->Branch("gj_c1_05_charged",   tjsev.gj_c1_05_charged,   "gj_c1_05_charged[ngj]/F");
  t->Branch("gj_c1_05_puppi",   tjsev.gj_c1_05_puppi,   "gj_c1_05_puppi[ngj]/F");
  t->Branch("gj_c1_05_all",   tjsev.gj_c1_05_all,   "gj_c1_05_all[ngj]/F");
  t->Branch("gj_c1_10_charged",   tjsev.gj_c1_10_charged,   "gj_c1_10_charged[ngj]/F");
  t->Branch("gj_c1_10_puppi",   tjsev.gj_c1_10_puppi,   "gj_c1_10_puppi[ngj]/F");
  t->Branch("gj_c1_10_all",   tjsev.gj_c1_10_all,   "gj_c1_10_all[ngj]/F");
  t->Branch("gj_c1_20_charged",   tjsev.gj_c1_20_charged,   "gj_c1_20_charged[ngj]/F");
  t->Branch("gj_c1_20_puppi",   tjsev.gj_c1_20_puppi,   "gj_c1_20_puppi[ngj]/F");
  t->Branch("gj_c1_20_all",   tjsev.gj_c1_20_all,   "gj_c1_20_all[ngj]/F");
  t->Branch("gj_c2_02_charged",   tjsev.gj_c2_02_charged,   "gj_c2_02_charged[ngj]/F");
  t->Branch("gj_c2_02_puppi",   tjsev.gj_c2_02_puppi,   "gj_c2_02_puppi[ngj]/F");
  t->Branch("gj_c2_02_all",   tjsev.gj_c2_02_all,   "gj_c2_02_all[ngj]/F");
  t->Branch("gj_c2_05_charged",   tjsev.gj_c2_05_charged,   "gj_c2_05_charged[ngj]/F");
  t->Branch("gj_c2_05_puppi",   tjsev.gj_c2_05_puppi,   "gj_c2_05_puppi[ngj]/F");
  t->Branch("gj_c2_05_all",   tjsev.gj_c2_05_all,   "gj_c2_05_all[ngj]/F");
  t->Branch("gj_c2_10_charged",   tjsev.gj_c2_10_charged,   "gj_c2_10_charged[ngj]/F");
  t->Branch("gj_c2_10_puppi",   tjsev.gj_c2_10_puppi,   "gj_c2_10_puppi[ngj]/F");
  t->Branch("gj_c2_10_all",   tjsev.gj_c2_10_all,   "gj_c2_10_all[ngj]/F");
  t->Branch("gj_c2_20_charged",   tjsev.gj_c2_20_charged,   "gj_c2_20_charged[ngj]/F");
  t->Branch("gj_c2_20_puppi",   tjsev.gj_c2_20_puppi,   "gj_c2_20_puppi[ngj]/F");
  t->Branch("gj_c2_20_all",   tjsev.gj_c2_20_all,   "gj_c2_20_all[ngj]/F");
  t->Branch("gj_c3_02_charged",   tjsev.gj_c3_02_charged,   "gj_c3_02_charged[ngj]/F");
  t->Branch("gj_c3_02_puppi",   tjsev.gj_c3_02_puppi,   "gj_c3_02_puppi[ngj]/F");
  t->Branch("gj_c3_02_all",   tjsev.gj_c3_02_all,   "gj_c3_02_all[ngj]/F");
  t->Branch("gj_c3_05_charged",   tjsev.gj_c3_05_charged,   "gj_c3_05_charged[ngj]/F");
  t->Branch("gj_c3_05_puppi",   tjsev.gj_c3_05_puppi,   "gj_c3_05_puppi[ngj]/F");
  t->Branch("gj_c3_05_all",   tjsev.gj_c3_05_all,   "gj_c3_05_all[ngj]/F");
  t->Branch("gj_c3_10_charged",   tjsev.gj_c3_10_charged,   "gj_c3_10_charged[ngj]/F");
  t->Branch("gj_c3_10_puppi",   tjsev.gj_c3_10_puppi,   "gj_c3_10_puppi[ngj]/F");
  t->Branch("gj_c3_10_all",   tjsev.gj_c3_10_all,   "gj_c3_10_all[ngj]/F");
  t->Branch("gj_c3_20_charged",   tjsev.gj_c3_20_charged,   "gj_c3_20_charged[ngj]/F");
  t->Branch("gj_c3_20_puppi",   tjsev.gj_c3_20_puppi,   "gj_c3_20_puppi[ngj]/F");
  t->Branch("gj_c3_20_all",   tjsev.gj_c3_20_all,   "gj_c3_20_all[ngj]/F");

}

//
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev)
{
  tjsev.nw=0;   tjsev.nl=0;   tjsev.nj=0;   tjsev.ngj=0;   tjsev.ngl=0;   tjsev.met_pt=0; tjsev.met_phi=0;
  for(int i=0; i<10; i++) tjsev.weight[i]=0;
  for(int i=0; i<5; i++) { tjsev.l_pt[i]=0;   tjsev.l_eta[i]=0;   tjsev.l_phi[i]=0;   tjsev.l_m[i]=0; tjsev.l_id[i]=0; tjsev.gl_pt[i]=0;   tjsev.gl_eta[i]=0;   tjsev.gl_phi[i]=0;   tjsev.gl_m[i]=0; tjsev.gl_id[i]=0; }
  for(int i=0; i<50; i++) {
    tjsev.j_pt[i]=0;   tjsev.j_eta[i]=0;   tjsev.j_phi[i]=0;   tjsev.j_m[i]=0; tjsev.j_flavor[i]=0; tjsev.j_overlap[i]=0; tjsev.j_gj[i]=-1;
    tjsev.gj_pt[i]=0;   tjsev.gj_eta[i]=0;   tjsev.gj_phi[i]=0;   tjsev.gj_m[i]=0; tjsev.gj_flavor[i]=0; tjsev.gj_overlap[i]=0; tjsev.gj_j[i]=-1;
    
  tjsev.j_mult_charged[i]=-99;
  tjsev.j_mult_puppi[i]=-99;
  tjsev.j_mult_all[i]=-99;
  tjsev.j_width_charged[i]=-99;
  tjsev.j_width_puppi[i]=-99;
  tjsev.j_width_all[i]=-99;
  tjsev.j_ptd_charged[i]=-99;
  tjsev.j_ptd_puppi[i]=-99;
  tjsev.j_ptd_all[i]=-99;
  tjsev.j_ecc_charged[i]=-99;
  tjsev.j_ecc_puppi[i]=-99;
  tjsev.j_ecc_all[i]=-99;
  tjsev.j_tau21_charged[i]=-99;
  tjsev.j_tau21_puppi[i]=-99;
  tjsev.j_tau21_all[i]=-99;
  tjsev.j_tau32_charged[i]=-99;
  tjsev.j_tau32_puppi[i]=-99;
  tjsev.j_tau32_all[i]=-99;
  tjsev.j_tau43_charged[i]=-99;
  tjsev.j_tau43_puppi[i]=-99;
  tjsev.j_tau43_all[i]=-99;
  tjsev.j_zg_charged[i]=-99;
  tjsev.j_zg_puppi[i]=-99;
  tjsev.j_zg_all[i]=-99;
  tjsev.j_zgxdr_charged[i]=-99;
  tjsev.j_zgxdr_puppi[i]=-99;
  tjsev.j_zgxdr_all[i]=-99;
  tjsev.j_zgdr_charged[i]=-99;
  tjsev.j_zgdr_puppi[i]=-99;
  tjsev.j_zgdr_all[i]=-99;
  tjsev.j_ga_width_charged[i]=-99;
  tjsev.j_ga_width_puppi[i]=-99;
  tjsev.j_ga_width_all[i]=-99;
  tjsev.j_ga_lha_charged[i]=-99;
  tjsev.j_ga_lha_puppi[i]=-99;
  tjsev.j_ga_lha_all[i]=-99;
  tjsev.j_ga_thrust_charged[i]=-99;
  tjsev.j_ga_thrust_puppi[i]=-99;
  tjsev.j_ga_thrust_all[i]=-99;
  tjsev.j_c1_02_charged[i]=-99;
  tjsev.j_c1_02_puppi[i]=-99;
  tjsev.j_c1_02_all[i]=-99;
  tjsev.j_c1_05_charged[i]=-99;
  tjsev.j_c1_05_puppi[i]=-99;
  tjsev.j_c1_05_all[i]=-99;
  tjsev.j_c1_10_charged[i]=-99;
  tjsev.j_c1_10_puppi[i]=-99;
  tjsev.j_c1_10_all[i]=-99;
  tjsev.j_c1_20_charged[i]=-99;
  tjsev.j_c1_20_puppi[i]=-99;
  tjsev.j_c1_20_all[i]=-99;
  tjsev.j_c2_02_charged[i]=-99;
  tjsev.j_c2_02_puppi[i]=-99;
  tjsev.j_c2_02_all[i]=-99;
  tjsev.j_c2_05_charged[i]=-99;
  tjsev.j_c2_05_puppi[i]=-99;
  tjsev.j_c2_05_all[i]=-99;
  tjsev.j_c2_10_charged[i]=-99;
  tjsev.j_c2_10_puppi[i]=-99;
  tjsev.j_c2_10_all[i]=-99;
  tjsev.j_c2_20_charged[i]=-99;
  tjsev.j_c2_20_puppi[i]=-99;
  tjsev.j_c2_20_all[i]=-99;
  tjsev.j_c3_02_charged[i]=-99;
  tjsev.j_c3_02_puppi[i]=-99;
  tjsev.j_c3_02_all[i]=-99;
  tjsev.j_c3_05_charged[i]=-99;
  tjsev.j_c3_05_puppi[i]=-99;
  tjsev.j_c3_05_all[i]=-99;
  tjsev.j_c3_10_charged[i]=-99;
  tjsev.j_c3_10_puppi[i]=-99;
  tjsev.j_c3_10_all[i]=-99;
  tjsev.j_c3_20_charged[i]=-99;
  tjsev.j_c3_20_puppi[i]=-99;
  tjsev.j_c3_20_all[i]=-99;

  tjsev.gj_mult_charged[i]=-99;
  tjsev.gj_mult_puppi[i]=-99;
  tjsev.gj_mult_all[i]=-99;
  tjsev.gj_width_charged[i]=-99;
  tjsev.gj_width_puppi[i]=-99;
  tjsev.gj_width_all[i]=-99;
  tjsev.gj_ptd_charged[i]=-99;
  tjsev.gj_ptd_puppi[i]=-99;
  tjsev.gj_ptd_all[i]=-99;
  tjsev.gj_ecc_charged[i]=-99;
  tjsev.gj_ecc_puppi[i]=-99;
  tjsev.gj_ecc_all[i]=-99;
  tjsev.gj_tau21_charged[i]=-99;
  tjsev.gj_tau21_puppi[i]=-99;
  tjsev.gj_tau21_all[i]=-99;
  tjsev.gj_tau32_charged[i]=-99;
  tjsev.gj_tau32_puppi[i]=-99;
  tjsev.gj_tau32_all[i]=-99;
  tjsev.gj_tau43_charged[i]=-99;
  tjsev.gj_tau43_puppi[i]=-99;
  tjsev.gj_tau43_all[i]=-99;
  tjsev.gj_zg_charged[i]=-99;
  tjsev.gj_zg_puppi[i]=-99;
  tjsev.gj_zg_all[i]=-99;
  tjsev.gj_zgxdr_charged[i]=-99;
  tjsev.gj_zgxdr_puppi[i]=-99;
  tjsev.gj_zgxdr_all[i]=-99;
  tjsev.gj_zgdr_charged[i]=-99;
  tjsev.gj_zgdr_puppi[i]=-99;
  tjsev.gj_zgdr_all[i]=-99;
  tjsev.gj_ga_width_charged[i]=-99;
  tjsev.gj_ga_width_puppi[i]=-99;
  tjsev.gj_ga_width_all[i]=-99;
  tjsev.gj_ga_lha_charged[i]=-99;
  tjsev.gj_ga_lha_puppi[i]=-99;
  tjsev.gj_ga_lha_all[i]=-99;
  tjsev.gj_ga_thrust_charged[i]=-99;
  tjsev.gj_ga_thrust_puppi[i]=-99;
  tjsev.gj_ga_thrust_all[i]=-99;
  tjsev.gj_c1_02_charged[i]=-99;
  tjsev.gj_c1_02_puppi[i]=-99;
  tjsev.gj_c1_02_all[i]=-99;
  tjsev.gj_c1_05_charged[i]=-99;
  tjsev.gj_c1_05_puppi[i]=-99;
  tjsev.gj_c1_05_all[i]=-99;
  tjsev.gj_c1_10_charged[i]=-99;
  tjsev.gj_c1_10_puppi[i]=-99;
  tjsev.gj_c1_10_all[i]=-99;
  tjsev.gj_c1_20_charged[i]=-99;
  tjsev.gj_c1_20_puppi[i]=-99;
  tjsev.gj_c1_20_all[i]=-99;
  tjsev.gj_c2_02_charged[i]=-99;
  tjsev.gj_c2_02_puppi[i]=-99;
  tjsev.gj_c2_02_all[i]=-99;
  tjsev.gj_c2_05_charged[i]=-99;
  tjsev.gj_c2_05_puppi[i]=-99;
  tjsev.gj_c2_05_all[i]=-99;
  tjsev.gj_c2_10_charged[i]=-99;
  tjsev.gj_c2_10_puppi[i]=-99;
  tjsev.gj_c2_10_all[i]=-99;
  tjsev.gj_c2_20_charged[i]=-99;
  tjsev.gj_c2_20_puppi[i]=-99;
  tjsev.gj_c2_20_all[i]=-99;
  tjsev.gj_c3_02_charged[i]=-99;
  tjsev.gj_c3_02_puppi[i]=-99;
  tjsev.gj_c3_02_all[i]=-99;
  tjsev.gj_c3_05_charged[i]=-99;
  tjsev.gj_c3_05_puppi[i]=-99;
  tjsev.gj_c3_05_all[i]=-99;
  tjsev.gj_c3_10_charged[i]=-99;
  tjsev.gj_c3_10_puppi[i]=-99;
  tjsev.gj_c3_10_all[i]=-99;
  tjsev.gj_c3_20_charged[i]=-99;
  tjsev.gj_c3_20_puppi[i]=-99;
  tjsev.gj_c3_20_all[i]=-99;

  } 
  
  tjsev.gen_sel = -1; tjsev.reco_sel = -1;
}
