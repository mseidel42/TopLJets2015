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
  bool requireEtriggerOnly(false);
  if(ev.isData && filename.Contains("SingleElectron")) requireEtriggerOnly=true;

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

  nentries = 100000;
  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetTopJetShapeEvent(tjsev);
      if(iev%100==0) printf ("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      
      //GENERATOR LEVEL
      
      int sel_ngleptons = 0;
      int vet_ngleptons = 0;
      int sel_ngbjets   = 0;
      int sel_ngjets    = 0;
      std::vector<Jet> genJets;
      
      //loop over gen objects from pseudotop producer
      for (int i = 0; i < ev.ng; i++) {
        //jets
        if (ev.g_pt[i]>30 && abs(ev.g_eta[i])<2.4 && abs(ev.g_id[i])<10) {
          //store to jet analysis objects
          TLorentzVector jp4; jp4.SetPtEtaPhiM(ev.g_pt[i],ev.g_eta[i],ev.g_phi[i],ev.g_m[i]);
          genJets.push_back(Jet(jp4, abs(ev.g_id[i]), i));
          //count
          tjsev.ngj++;
          if      (abs(ev.g_id[i]) <5) sel_ngjets++;
          else if (abs(ev.g_id[i])==5) sel_ngbjets++;
        }
        //leptons
        else if (abs(ev.g_id[i])>10) {
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
          else if (ev.g_pt[i]>10 && abs(ev.g_eta[i])<2.5) vet_ngleptons++;
        }
      }
      
      //event selected on generator level?
      if (sel_ngbjets==2 && sel_ngjets+sel_ngbjets>=4 &&
          sel_ngleptons==1 && vet_ngleptons==0) tjsev.gen_sel = 1;
      
      //flag jets as part of W boson candidates
      for (int i = 0; i < tjsev.ngj; i++) {
        if (genJets[i].flav==5) continue;
        for (int j = i+1; j < tjsev.ngj; j++) {
          if (genJets[j].flav==5) continue;
          TLorentzVector wCand = genJets[i].p4 + genJets[j].p4;
          if (abs(wCand.M()-80.4) < 10.) {
            genJets[i].wcandidate++;
            genJets[j].wcandidate++;
          }
        }
      }
      
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
      
      //calculate generalized angularities (jet shapes)
      for (int i = 0; i < tjsev.ngj; i++) {
        for (int i1 = 0; i1 < 3; ++i1) for (int i2 = 0; i2 < 3; ++i2) for (int i3 = 0; i3 < 3; ++i3) for (int i4 = 0; i4 < 3; ++i4)
        tjsev.gj_ga[i][i1][i2][i3][i4] = calcGA(genJets[i], i1, i2, i3, i4);
      }
      
      //store analyzed jets to tree
      for (int i = 0; i < tjsev.ngj; i++) {
        tjsev.gj_pt  [i] = genJets[i].p4.Pt();
        tjsev.gj_eta [i] = genJets[i].p4.Eta();
        tjsev.gj_phi [i] = genJets[i].p4.Phi();
        tjsev.gj_m   [i] = genJets[i].p4.M();
        tjsev.gj_flav[i] = genJets[i].flav;
        tjsev.gj_w   [i] = genJets[i].wcandidate;
      }

      //RECO LEVEL
      /*
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
          bool passLooseKin(ev.l_pt[il]>10 && fabs(ev.l_eta[il])<2.5);
          bool passTightId(ev.l_id[il]==13 ? (ev.l_pid[il]>>1)&0x1  : (ev.l_pid[il]>>2)&0x1);
          float relIso(ev.l_relIso[il]);
          bool passTightIso( ev.l_id[il]==13 ? relIso<0.15 : (ev.l_pid[il]>>1)&0x1);
          bool passLooseIso( ev.l_id[il]==13 ? relIso<0.25 : (ev.l_pid[il]   )&0x1);
          if(passTightKin && passTightId && passTightIso) tightLeptons.push_back(il);
          else if(passLooseKin && passLooseIso)           looseLeptons.push_back(il);
        }
      
      //check if triggers have fired
      bool hasMuTrigger((ev.muTrigger & 0x3)!=0);
      bool hasEleTrigger((ev.elTrigger & 0x1)!=0);

      //decide the channel
      TString chTag("");
      std::vector<int> selLeptons;
      if(tightLeptons.size()==1 && looseLeptons.size()==0)
        {
          selLeptons.push_back( tightLeptons[0] );
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
      if(selLeptons.size()==0) continue;
      if(selLeptons.size()==1)
        {
          if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="E";
          if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger)  chTag="M";
        }
      if(selLeptons.size()==2)
        {
          if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*11 && hasEleTrigger) chTag="EE";
          if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==13*13 && hasMuTrigger) chTag="MM";
          if(abs(ev.l_id[ selLeptons[0] ]*ev.l_id[ selLeptons[1] ])==11*13)
            {
              if(tightLeptons.size()>=2 && (hasEleTrigger || hasMuTrigger)) chTag="EM";
              if(tightLeptons.size()==1)
                {
                  if(abs(ev.l_id[ selLeptons[0] ])==11 && hasEleTrigger) chTag="EM";
                  if(abs(ev.l_id[ selLeptons[0] ])==13 && hasMuTrigger) chTag="EM";
                }
            }
          if(hasMuTrigger && requireEtriggerOnly) chTag="";
        }
      if(chTag=="") continue;

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
      std::vector<int> genBJetsFlav,genOtherjetsFlav,genBJetsHadFlav,genOtherjetsHadFlav;
      std::vector<TLorentzVector> bJets,otherjets, genBJets, genOtherjets;
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
              float jerSmear=getJetResolutionScales(jp4.Pt(),jp4.Pt(),genJet_pt)[0];
              jp4 *= jerSmear;
            }

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

          if(jp4.Pt()<30) continue;

          Int_t hadFlav=ev.j_hadflav[k];
          Int_t flav=ev.j_flav[k];
          TLorentzVector gjp4(0,0,0,0);
          if(ev.j_g[k]>-1)
            {
              int gidx=ev.j_g[k];
              gjp4.SetPtEtaPhiM( ev.g_pt[gidx], ev.g_eta[gidx], ev.g_phi[gidx], ev.g_m[gidx] );
            }

          //save jets 
          if(isBTagged && fabs(jp4.Eta()) < 2.4 && bJets.size()<2) 
            { bJets.push_back(jp4);     genBJets.push_back(gjp4);     genBJetsFlav.push_back(flav);     genBJetsHadFlav.push_back(hadFlav);    }
          else 
            { otherjets.push_back(jp4); genOtherjets.push_back(gjp4); genOtherjetsFlav.push_back(flav); genOtherjetsHadFlav.push_back(hadFlav); }
        }

      //2 b-jets are required 
      //notice more can be b-tagged but they are put into the other jets category
      if(bJets.size()!=2) continue;

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
      allPlots["nvtx_"+chTag]->Fill(ev.nvtx,wgt);
      allPlots["njets_"+chTag]->Fill(otherjets.size(),wgt);

      tjsev.nl=leptons.size();
      for(int il=0; il<(int)leptons.size(); il++)
        {
          TString pf(Form("l%d",il));
          allPlots[pf+"pt_"+chTag]->Fill(leptons[il].Pt(),wgt);
          allPlots[pf+"eta_"+chTag]->Fill(fabs(leptons[il].Eta()),wgt);
          tjsev.l_pt[il]=leptons[il].Pt();
          tjsev.l_eta[il]=leptons[il].Eta();
          tjsev.l_phi[il]=leptons[il].Phi();
          tjsev.l_m[il]=leptons[il].M();
          tjsev.l_id[il]=ev.l_id[ selLeptons[il] ];
          for(Int_t ig=0; ig<ev.ng; ig++)
            {
              if(abs(ev.g_id[ig])!=ev.l_id[ selLeptons[il] ]) continue;
              TLorentzVector glp4;
              glp4.SetPtEtaPhiM( ev.g_pt[ig], ev.g_eta[ig], ev.g_phi[ig], ev.g_m[ig]);
              if(glp4.DeltaR( leptons[il] ) > 0.3) continue;
              tjsev.gl_id[il]=ev.g_id[ig];
              tjsev.gl_pt[il]=ev.g_pt[ig];
              tjsev.gl_eta[il]=ev.g_eta[ig];
              tjsev.gl_phi[il]=ev.g_phi[ig];
              tjsev.gl_m[il]=ev.g_m[ig];
            }
        }
      
      tjsev.nj=bJets.size();
      for(int ij=0; ij<(int)bJets.size(); ij++)
        {
          TString pf(Form("j%d",ij));
          allPlots[pf+"pt_"+chTag]->Fill(bJets[ij].Pt(),wgt);
          allPlots[pf+"eta_"+chTag]->Fill(fabs(bJets[ij].Eta()),wgt);
          tjsev.j_pt[ij]=bJets[ij].Pt();
          tjsev.j_eta[ij]=bJets[ij].Eta();
          tjsev.j_phi[ij]=bJets[ij].Phi();
          tjsev.j_m[ij]=bJets[ij].M();
          tjsev.gj_flav[ij]=genBJetsFlav[ij];
          tjsev.gj_hadflav[ij]=genBJetsHadFlav[ij];
          tjsev.gj_pt[ij]=genBJets[ij].Pt();
          tjsev.gj_eta[ij]=genBJets[ij].Eta();
          tjsev.gj_phi[ij]=genBJets[ij].Phi();
          tjsev.gj_m[ij]=genBJets[ij].M();          
        }
      
      tjsev.nj+=otherjets.size();
      for(int ij=0; ij<(int)otherjets.size(); ij++)
        {
          if(ij+2<6)
            {
              TString pf(Form("j%d",ij+2));
              allPlots[pf+"pt_"+chTag]->Fill(otherjets[ij].Pt(),wgt);
              allPlots[pf+"eta_"+chTag]->Fill(fabs(otherjets[ij].Eta()),wgt);
            }
          tjsev.j_pt[ij+2]=otherjets[ij].Pt();
          tjsev.j_eta[ij+2]=otherjets[ij].Eta();
          tjsev.j_phi[ij+2]=otherjets[ij].Phi();
          tjsev.j_m[ij+2]=otherjets[ij].M();
          tjsev.gj_flav[ij+2]=genOtherjetsFlav[ij];
          tjsev.gj_hadflav[ij+2]=genOtherjetsHadFlav[ij];
          tjsev.gj_pt[ij+2]=genOtherjets[ij].Pt();
          tjsev.gj_eta[ij+2]=genOtherjets[ij].Eta();
          tjsev.gj_phi[ij+2]=genOtherjets[ij].Phi();
          tjsev.gj_m[ij+2]=genOtherjets[ij].M();          
        }

      tjsev.cat=11;
      if(chTag=="M") tjsev.cat=13;
      if(chTag=="MM") tjsev.cat=13*13;
      if(chTag=="EM") tjsev.cat=11*13;
      if(chTag=="EE") tjsev.cat=11*11;
      tjsev.nw=1;
      tjsev.weight[0]=wgt;
      tjsev.nt=0;
      tjsev.met_pt=ev.met_pt[0];
      tjsev.met_phi=ev.met_phi[0];
      if(ev.ngtop>0)
        {
          for(int i=0; i<ev.ngtop; i++)
            {
              if(abs(ev.gtop_id[i])!=6) continue;
              tjsev.t_pt[tjsev.nt]=ev.gtop_pt[i];
              tjsev.t_eta[tjsev.nt]=ev.gtop_eta[i];
              tjsev.t_phi[tjsev.nt]=ev.gtop_phi[i];
              tjsev.t_m[tjsev.nt]=ev.gtop_m[i];
              tjsev.t_id[tjsev.nt]=ev.gtop_id[i];
              tjsev.nt++;
              if(tjsev.nt>4) break;
            }
        }
      */
      int reco_sel = -1;
      if (tjsev.gen_sel + reco_sel > -1) outT->Fill();
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
double calcGA(Jet jet, int beta, int kappa, int iptcut, int icharge) {
  double sumpt = 0;
  for (auto p : jet.particles) {
    if (p.p4.Pt() < iptcut*0.500) continue;
    if      (icharge == 0 && p.charge!=0) sumpt += p.p4.Pt();
    else if (icharge == 1) sumpt += p.p4.Pt();
    else if (icharge == 2) sumpt += p.p4.Pt()*p.puppi;
  }
  //std::cout << "sumpt" << beta << kappa << iptcut << icharge << ": " << sumpt << std::endl;
  
  double ga = 0;
  for (auto p : jet.particles) {
    if (p.p4.Pt() < iptcut*0.500) continue;
    if(icharge == 0 && p.charge==0) continue;
    double weight = 1.;
    if (icharge == 2) weight = p.puppi;
    ga += weight * pow(p.p4.Pt()/sumpt, kappa) * pow(jet.p4.DeltaR(p.p4)/0.4, beta);
  }
  
  //std::cout << "ga" << beta << kappa << iptcut << icharge << ": " << ga << std::endl;
  return ga;
}

//
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev)
{
  //event category
  t->Branch("cat", &tjsev.cat,"cat/I");

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
  
  //gen jets
  t->Branch("ngj",  &tjsev.ngj, "ngj/I");
  t->Branch("gj_pt",  tjsev.gj_pt ,  "gj_pt[ngj]/F");
  t->Branch("gj_eta", tjsev.gj_eta , "gj_eta[ngj]/F");
  t->Branch("gj_phi", tjsev.gj_phi , "gj_phi[ngj]/F");
  t->Branch("gj_m",   tjsev.gj_m ,   "gj_m[ngj]/F");
  t->Branch("gj_flav",  tjsev.gj_flav ,  "gj_flav[ngj]/I");
  t->Branch("gj_hadflav",  tjsev.gj_hadflav ,  "gj_hadflav[ngj]/I");
  t->Branch("gj_w",  tjsev.gj_w ,  "gj_w[ngj]/I");
  t->Branch("gj_ga",  tjsev.gj_ga ,  "gj_ga[ngj][3][3][3][3]/F");

  //mc truth
  t->Branch("nt",  &tjsev.nt, "nt/I");
  t->Branch("t_pt",  tjsev.t_pt ,  "t_pt[nt]/F");
  t->Branch("t_eta", tjsev.t_eta , "t_eta[nt]/F");
  t->Branch("t_phi", tjsev.t_phi , "t_phi[nt]/F");
  t->Branch("t_m",   tjsev.t_m ,   "t_m[nt]/F");
  t->Branch("t_id",  tjsev.t_id ,  "t_id[nt]/I");
  
  t->Branch("gen_sel", &tjsev.gen_sel ,  "gen_sel/I");
}

//
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev)
{
  tjsev.cat=0;   tjsev.nw=0;   tjsev.nl=0;   tjsev.nj=0;   tjsev.ngj=0;   tjsev.ngl=0;   tjsev.nt=0;
  tjsev.met_pt=0; tjsev.met_phi=0;
  for(int i=0; i<10; i++) tjsev.weight[i]=0;
  for(int i=0; i<5; i++) { tjsev.l_pt[i]=0;   tjsev.l_eta[i]=0;   tjsev.l_phi[i]=0;   tjsev.l_m[i]=0; tjsev.l_id[i]=0; tjsev.gl_pt[i]=0;   tjsev.gl_eta[i]=0;   tjsev.gl_phi[i]=0;   tjsev.gl_m[i]=0; tjsev.gl_id[i]=0; }
  for(int i=0; i<50; i++) {
    tjsev.j_pt[i]=0;   tjsev.j_eta[i]=0;   tjsev.j_phi[i]=0;   tjsev.j_m[i]=0; tjsev.gj_pt[i]=0;   tjsev.gj_eta[i]=0;   tjsev.gj_phi[i]=0;   tjsev.gj_m[i]=0; tjsev.gj_flav[i]=0; tjsev.gj_hadflav[i]=0; tjsev.gj_w[i]=0;
    for (int i1 = 0; i1 < 3; ++i1) for (int i2 = 0; i2 < 3; ++i2) for (int i3 = 0; i3 < 3; ++i3) for (int i4 = 0; i4 < 3; ++i4) tjsev.gj_ga[i][i1][i2][i3][i4] = 0;
  } 
  for(int i=0; i<4; i++) { tjsev.t_pt[i]=0;   tjsev.t_eta[i]=0;   tjsev.t_phi[i]=0;   tjsev.t_m[i]=0; tjsev.t_id[i]=0; }
  tjsev.gen_sel = -1;
}
