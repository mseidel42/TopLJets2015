#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"
#include "TopLJets2015/TopAnalysis/interface/analysisXCone.h"
#include "TopLJets2015/TopAnalysis/interface/CorrectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/LeptonEfficiencyWrapper.h"
#include "TopLJets2015/TopAnalysis/interface/SelectionTools.h"
#include "TopLJets2015/TopAnalysis/interface/CommonTools.h"

#include "fastjet/tools/Recluster.hh"
#include "fastjet/contrib/Njettiness.hh"
using fastjet::contrib::Njettiness;
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh" 
using fastjet::contrib::Nsubjettiness;
using fastjet::contrib::OnePass_WTA_KT_Axes;
using fastjet::contrib::UnnormalizedMeasure;
using fastjet::contrib::XConeMeasure;


#include <vector>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TKey.h"

using namespace std;

//
void RunanalysisXCone(TString filename,
              TString outname,
              Int_t channelSelection, 
              Int_t chargeSelection, 
              SelectionTool::FlavourSplitting flavourSplitting,
              TH1F *normH, 
              Bool_t runSysts,
              std::string systVar,
              TString era,
              Bool_t debug)
{
  /////////////////////
  // INITIALIZATION //
  ///////////////////

  bool isTTbar( filename.Contains("_TTJets") or (normH and TString(normH->GetTitle()).Contains("_TTJets")));
  bool isData( filename.Contains("Data") );
  
  // explicit systematics
  std::vector<std::string> vSystVar;
  boost::split(vSystVar, systVar, boost::is_any_of("_"));

  //READ TREE FROM FILE
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TH1 *genPU=(TH1 *)f->Get("analysis/putrue");
  TH1 *triggerList=(TH1 *)f->Get("analysis/triggerList");
  TTree *t = (TTree*)f->Get("analysis/data");
  attachToMiniEventTree(t,ev,true);
  Int_t nentries(t->GetEntriesFast());
  if (debug) nentries = 10000;

  //EVENT SELECTION WRAPPER
  SelectionTool evsel(filename,false,triggerList);

  //PREPARE CORRECTIONS  
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  
  //LUMI
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isData )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }
  
  //PILEUP WEIGHTING
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isData ) puWgtGr=getPileupWeightsMap(era,genPU);
  
  //LEPTON EFFICIENCIES
  LeptonEfficiencyWrapper lepEffH(filename.Contains("Data13TeV"),era+"GH");

  //B-TAG CALIBRATION
  BTagSFUtil* myBTagSFUtil = new BTagSFUtil();
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  //dummy calls
  btvsfReaders["GH"][BTagEntry::FLAV_B]->eval_auto_bounds("central", BTagEntry::FLAV_B,   0., 30.);
  btvsfReaders["GH"][BTagEntry::FLAV_UDSG]->eval_auto_bounds("central", BTagEntry::FLAV_UDSG,   0., 30.);

  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *>    expBtagEffPy8 = readExpectedBtagEff(era);
  TString btagExpPostFix("");
  if(isTTbar) {
    if(filename.Contains("_herwig"))    btagExpPostFix="_herwig";
    if(filename.Contains("_scaleup"))   btagExpPostFix="_scaleup";
    if(filename.Contains("_scaledown")) btagExpPostFix="_scaledown";
  }
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff=readExpectedBtagEff(era,btagExpPostFix);
  
  //JET ENERGY UNCERTAINTIES
  std::string jecVar = "Total";
  if (vSystVar[0] == "jec") {
    if (vSystVar.size() != 3) {
      std::cout << "Wrong format of JEC uncertainty, expected jec_Source_up/down. Exiting..." << std::endl;
      return;
    }
    jecVar = vSystVar[1];
  }
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), jecVar);
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
  
  //JER SMEARING
  TString resolutionFile(era+"/Spring16_25nsV10_MC_PtResolution_AK4PFchs.txt");
  gSystem->ExpandPathName(resolutionFile);
  JME::JetResolution *jer = new JME::JetResolution(resolutionFile.Data());
  
  //BFRAG UNCERTAINTIES
  std::map<TString, TGraph*> bfrag = getBFragmentationWeights(era);
  std::map<TString, std::map<int, double> > semilepbr = getSemilepBRWeights(era);

  //PREPARE OUTPUT
  analysisXCone_t tue;
  TFile *fOut=TFile::Open(outname,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tue","tue");
  createanalysisXConeTree(outT,tue);
  outT->SetDirectory(fOut);

  //BOOK CONTROL HISTOGRAMS
  HistTool ht; 
  if (isData) ht.setNsyst(0);
  ht.addHist("puwgtctr", new TH1F("puwgtctr","Weight sums",4,0,4) );  
  std::vector<TString> lfsVec = { "EM" };
  for(size_t ilfs=0; ilfs<lfsVec.size(); ilfs++)   
    { 
      TString tag(lfsVec[ilfs]);
      if(ratevsrunH) ht.addHist("ratevsrun_"+tag, (TH1 *)ratevsrunH->Clone("ratevsrun_"+tag) );
      ht.addHist("nvtx_"+tag, new TH1F("nvtx_"+tag,";Vertex multiplicity;Events",40,0,40) );
      ht.addHist("rho_"+tag, new TH1F("rho_"+tag,";#rho;Events",40,0,40));
      for(size_t i=0; i<=2; i++)
        {
          TString subtag(tag);
          if(i<2) { subtag += i; subtag += "t"; }
          ht.addHist("mll_"+subtag,  new TH1F("mll_"+subtag,";Dilepton invariant mass [GeV];Events",50,0,400) );
          ht.addHist("met_"+subtag       , new TH1F("met_"+subtag,";Missing transverse momentum [GeV];Events",50,0,300) );
        }
      // Jet multiplicity!
      ht.addHist("njets_"+tag, new TH1F("njets_"+tag,";Extra jet multiplicity;Events",7,0,7) );
      ht.addHist("ptttbar_"+tag   , new TH1F("ptttbar_"+tag,";p_{T}(t#bar{t}) [GeV];Events",50,0,200) );
      ht.addHist("sumpt_"+tag     , new TH1F("sumpt_"+tag,";Transverse momentum sum [GeV];Events",50,40,300) );      
      ht.addHist("ptpos_"+tag     , new TH1F("ptpos_"+tag,";Lepton transverse momentum [GeV];Events",50,20,200) );
      ht.addHist("ptll_"+tag      , new TH1F("ptll_"+tag,";Dilepton transverse momentum [GeV];Events",50,0,200) );
      ht.addHist("nbtags_"+tag    , new TH1F("nbtags_"+tag,";b-tag multiplicity;Events",3,0,3) );
      ht.addHist("nch_"+tag      , new TH1F("nch_"+tag,";Charged particle multiplicity;Events",50,0,200) ); 
      ht.addHist("XConeCut_"+tag, new TH1F("XConeCut_"+tag,"; Jet XConeCut; Events",10,0,10));
      ht.addHist("XConeDiff_"+tag, new TH1F("XConeDiff_"+tag,"; Jet XConeDiff; Events",10,0,10));
      for (Int_t i = 2;  i<10; i++){
        TString ii = std::to_string(i);
        ht.addHist(ii+"_Jettiness"+tag, new TH1F(ii+"_Jettiness"+tag,";tau ;Events",50,0,200));
        ht.addHist(ii+"_Jettiness_slope"+tag, new TH1F(ii+"_Jettiness_slope"+tag,";d tau/d N ;Events",20,-0.1,0));
      }
    }
  for (auto& it : ht.getPlots() )     { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : ht.get2dPlots() )   { it.second->Sumw2(); it.second->SetDirectory(0); }

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
      t->GetEntry(iev);
      resetanalysisXCone(tue);
      if(iev%int(nentries/100)==0) printf ("[%3.0f%%] done\n", 100.*(float)iev/(float)nentries);
      
      //assign a run period and correct the event accordingly
      TString period = assignRunPeriod(runPeriods);
      
      //////////////////
      // CORRECTIONS //
      ////////////////
      
      double csvm = 0.8484;
      if (vSystVar[0] == "csv") {
          if (vSystVar[1] == "heavy") {
              //heavy flavor uncertainty +/-3.5%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, 0.8726, csvm);
              if (vSystVar[2] == "down") addBTagDecisions(ev, 0.8190, csvm);
          }
          if (vSystVar[1] == "light") {
              //light flavor uncertainty +/-10%
              if (vSystVar[2] == "up")   addBTagDecisions(ev, csvm, 0.8557);
              if (vSystVar[2] == "down") addBTagDecisions(ev, csvm, 0.8415);
          }
      }
      else addBTagDecisions(ev, csvm, csvm);
      
      if(!ev.isData) {
        //jec
        if (vSystVar[0] == "jec") {
          applyJetCorrectionUncertainty(ev, jecUnc, jecVar, vSystVar[2]);
        }
        //jer
        if (vSystVar[0] == "jer") {
          smearJetEnergies(ev, jer, vSystVar[1]);
        }
        else smearJetEnergies(ev, jer);
        //b tagging
        if (vSystVar[0] == "btag") {
          if (vSystVar[1] == "heavy") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,vSystVar[2],"central");
          }
          if (vSystVar[1] == "light") {
            updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil,"central",vSystVar[2]);
          }
        }
        else updateBTagDecisions(ev, btvsfReaders[period],expBtagEff,expBtagEffPy8,myBTagSFUtil);
        //tracking efficiency
        if (vSystVar[0] == "tracking") {
          // "up": no correction
          // "down": apply twice
          if (vSystVar[1] == "down") applyTrackingEfficiencySF(ev, pow(lepEffH.getTrackingCorrection(ev.nvtx, period).first,2));
        }
        else applyTrackingEfficiencySF(ev, lepEffH.getTrackingCorrection(ev.nvtx, period).first);
      }

      ///////////////////////////
      // RECO LEVEL SELECTION //
      /////////////////////////

      //evsel.setDebug(true);
      TString chTag = evsel.flagFinalState(ev);
      if(chTag=="EM" )
        {
          // leptons
          std::vector<Particle> &leptons = evsel.getSelLeptons();
          
          TLorentzVector l1(leptons[0].p4()), l2(leptons[1].p4());
          TLorentzVector dil(l1+l2);
          float mll = dil.M();
          bool passLepPresel(mll>12
                             && (l1.Pt()>25 || l2.Pt()>25)
                             && (fabs(l1.Eta())<2.5 && fabs(l2.Eta())<2.5) );

          // jets
          std::vector<Jet> jets = evsel.getGoodJets(ev,25.,2.4,leptons);

          int nb = 0;
          for (auto& jet : jets) {
            if (jet.flavor() == 5) ++nb;
          }
          
          // select events with 2 leptons and 2 b jets
          if (passLepPresel and nb >= 2) tue.reco_sel = 1;
          
          // store some selection observables to tree
          tue.nj     = jets.size();
          tue.nb     = nb;
          tue.mll    = mll;
          tue.ptpos  = leptons[0].charge()>0 ? l1.Pt() : l2.Pt();
          tue.phipos = leptons[0].charge()>0 ? l1.Phi() : l2.Phi();
          tue.ptll   = dil.Pt();
          tue.phill  = dil.Phi();
          tue.sumpt  = l1.Pt()+l2.Pt();
          tue.dphill = TMath::Abs(l1.DeltaPhi(l2));
          
          ////////////////////
          // EVENT WEIGHTS //
          //////////////////
          
          float wgt(1.0);
          // Pairs for systematic uncertainty weights
          // double: weight value (divided by central weight)
          // bool: use weight for plotting, otherwise just save to tree
          std::vector<std::pair<double, bool> > varweights;
          std::vector<double> plotwgts;
          
          if (!ev.isData) {
            // norm weight
            wgt  = (normH? normH->GetBinContent(1) : 1.0);
            
            // pu weight
            double puWgt(puWgtGr[period][0]->Eval(ev.g_pu));
            if (std::isnan(puWgt)) puWgt = 1.;
            if (puWgt == 0.) puWgt = 1.;
            wgt *= puWgt;
            double puWgt1(puWgtGr[period][1]->Eval(ev.g_pu));
            if (std::isnan(puWgt1)) puWgt1 = 1.;
            varweights.push_back(std::make_pair(puWgt1/puWgt, true)); // 1
            double puWgt2(puWgtGr[period][2]->Eval(ev.g_pu));
            if (std::isnan(puWgt2)) puWgt2 = 1.;
            varweights.push_back(std::make_pair(puWgt2/puWgt, true)); // 2
            
            // lepton trigger*selection weights
            EffCorrection_t selSF(1.0,0.0),trigSF(1.0,0.0);
            trigSF=lepEffH.getTriggerCorrection(leptons,period);
            varweights.push_back(std::make_pair(1.+trigSF.second, true)); // 3
            varweights.push_back(std::make_pair(1.-trigSF.second, true)); // 4
            for(size_t il=0; il<2; il++) {
              EffCorrection_t sf = lepEffH.getOfflineCorrection(leptons[il].id(),leptons[il].pt(),leptons[il].eta(),period);
              selSF.second = sqrt( pow(selSF.first*sf.second,2)+pow(selSF.second*sf.first,2));
              selSF.first *= sf.first;
            }
	          varweights.push_back(std::make_pair(1.+selSF.second, true));  // 5
            varweights.push_back(std::make_pair(1.-selSF.second, true));  // 6
            wgt *= trigSF.first*selSF.first;
            
            // bfrag weights
            varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["upFrag"]), true));       // 7
            varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["downFrag"]), true));     // 8
            varweights.push_back(std::make_pair(computeBFragmentationWeight(ev, bfrag["PetersonFrag"]), true)); // 9
            // weights for semilep BR
            // simultaneous variation for all hadrons, average over particle and antiparticle
            // divide by mean weight from 100k events to avoid change in cross section
            varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrUp"], 0, true)/1.00480, true)); // 10
            varweights.push_back(std::make_pair(computeSemilepBRWeight(ev, semilepbr["semilepbrDown"], 0, true)/0.992632, true)); // 11
            
        	  //top pt weighting
            double topptsf = 1.0;
            if(isTTbar) {
              for (int igen=0; igen<ev.ngtop; igen++) {
                if(abs(ev.gtop_id[igen])!=6) continue;
                topptsf *= TMath::Exp(0.0615-0.0005*ev.gtop_pt[igen]);
              }
            }
            varweights.push_back(std::make_pair(topptsf, true)); // 12
            
            // lhe weights
            wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
            
            std::set<std::string> scalesForPlotter = {
              "id1002muR1muF2hdampmt272.7225",    // 13
              "id1003muR1muF0.5hdampmt272.7225",  // 14
              "id1004muR2muF1hdampmt272.7225",    // 15
              "id1005muR2muF2hdampmt272.7225",    // 16
              "id1007muR0.5muF1hdampmt272.7225",  // 18
              "id1009muR0.5muF0.5hdampmt272.7225" // 20
            };
            for (int i = 1; i < ev.g_nw; ++i) {
              bool forPlotter = (normH and scalesForPlotter.count(normH->GetXaxis()->GetBinLabel(i)) != 0);
              varweights.push_back(std::make_pair(ev.g_w[i]/ev.g_w[0], forPlotter));
            }
            
            tue.nw = 1 + varweights.size();
            tue.weight[0]=wgt;
            for (unsigned int i = 0; i < varweights.size(); ++i) {
              tue.weight[i+1] = varweights[i].first;
            }
            plotwgts.push_back(wgt);
            for (auto vw : varweights)
              if (vw.second)
                plotwgts.push_back(vw.first);
          }
          else {
            tue.nw=1;
            tue.weight[0]=wgt;
            plotwgts.push_back(wgt);
          }

          //////////////////////
          // FILL HISTOGRAMS //
          ////////////////////
          
          ht.fill("nvtx_"+chTag,ev.nvtx,plotwgts);
          ht.fill("rho_"+chTag,ev.rho,plotwgts);
          if(nb<2)
            {
              TString subTag(chTag);
              if(nb==0) subTag += "0t";
              if(nb==1) subTag += "1t";              
              ht.fill("mll_"+subTag,mll,plotwgts);
              ht.fill("met_"+subTag,evsel.getMET().Pt(),plotwgts);              
            }
            
          ht.fill("nbtags_"+chTag,nb,plotwgts);
          if(nb>=2) ht.fill("mll_"+chTag,mll,plotwgts);
          
          if (tue.reco_sel == 1) {
            std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
            if(rIt!=lumiMap.end() && ratevsrunH) ht.getPlots()["ratevsrun_"+chTag]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
            ht.fill("met_"+chTag,evsel.getMET().Pt(),plotwgts);
            ht.fill("njets_"+chTag,jets.size(),plotwgts);
            ht.fill("sumpt_"+chTag,tue.sumpt,plotwgts);             
            ht.fill("ptpos_"+chTag,tue.ptpos,plotwgts);
            ht.fill("ptll_"+chTag,tue.ptll,plotwgts);
          
            //////////////////////////
            // RECO LEVEL ANALYSIS //
            ////////////////////////
            
            // Collect charged particles into PseudoJet
            std::vector<fastjet::PseudoJet> fj_particles;
            for(int ipf = 0; ipf < ev.npf; ipf++) {
              if(ev.pf_c[ipf]==0) continue;

              TLorentzVector tkP4(0,0,0,0);
              tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);

              //fiducial cuts
              bool passKin(ev.pf_pt[ipf]>0.9 && fabs(ev.pf_eta[ipf])<2.5);
              if (not passKin) continue;
              
              //matching to leptons
              Double_t relDpt2lep(9999.);
              for(int ilep=0; ilep<2; ilep++) {
                float dR=tkP4.DeltaR(leptons[ilep].p4());
                if(dR>0.05) continue;
                float relDpt=fabs(leptons[ilep].pt()-tkP4.Pt())/leptons[ilep].pt();
                if(relDpt>relDpt2lep) continue;
                relDpt2lep=relDpt;
              }
              bool matchedToLepton(relDpt2lep<0.05);
              if (matchedToLepton) continue;
              
              fj_particles.push_back( fastjet::PseudoJet( tkP4.Px(),tkP4.Py(), tkP4.Pz(), tkP4.E() ) );
            }
            
            ht.fill("nch_"+chTag, fj_particles.size(), plotwgts);
            
            // NJettiness variables
            // TODO: maybe move again to a place where it can be used by both reco and gen analysis
            double tau_cut = 30.0 ;
            double beta = 2.;
            double R = 0.4;
            double tau_temp_C = 1000.;
            double tau_temp_S = 1000.0;
            unsigned int N_temp_S = 1;
            unsigned int N_temp_C = 1;
            
            Njettiness* tauN = new fastjet::contrib::Njettiness(OnePass_WTA_KT_Axes(), XConeMeasure(beta,R));        
            for (int N = 2; N < 10; ++N) {
              double tau_C = tauN->getTau(N, fj_particles);
              double tau_S = 0;
              if (tauN->getTau(N-1,fj_particles) > 0.) {
                tau_S = (tauN->getTau(N,fj_particles) - tauN->getTau(N-1,fj_particles))/(tauN->getTau(N-1,fj_particles));
                if (tau_S == -1.0){
                  tau_S = 0.0;
                }
              }

              ht.fill(N+"_Jettiness"+chTag,tau_C,plotwgts);
              ht.fill(N+"_Jettiness_slope"+chTag,tau_S,plotwgts);

              if (tau_C < tau_temp_C and tau_C > tau_cut){
                tau_temp_C = tau_C;
                N_temp_C = N;
              }
              if (tau_S < tau_temp_S and tau_S !=-1.0){
                tau_temp_S = tau_S;
                N_temp_S = N;
              }
            }
            delete tauN;
            
            // Fill NJettiness results into tree and histograms
            
            tue.NJettiness = int( N_temp_C);
            tue.NJettSlope = int (N_temp_S);
            
            ht.fill("XConeCut_"+chTag,int(N_temp_C),plotwgts);
            ht.fill("XConeDiff_"+chTag,int(N_temp_S),plotwgts);
          }
        }

      //////////////////////////
      // GEN LEVEL SELECTION //
      ////////////////////////
      
      TString genChTag = evsel.flagGenFinalState(ev);
      if(isTTbar && genChTag=="EM")
        {
          // leptons
          std::vector<Particle> &leptons = evsel.getGenLeptons();
          TLorentzVector dil(leptons[0].p4()+leptons[1].p4());
          float mll = dil.M();
          bool passLepPresel(mll>12
                             && (leptons[0].pt()>25 || leptons[1].pt()>25)
                             && (fabs(leptons[0].eta())<2.5 && fabs(leptons[1].eta())<2.5) );
          
          // jets
          std::vector<Jet> &jets=evsel.getGenJets();

          int nb = 0;
          for (auto& jet : jets) {
            if (jet.flavor() == 5) ++nb;
          }
          
          // select events with 2 leptons and 2 b jets
          if (passLepPresel and nb >= 2) tue.gen_sel = 1;

          // fill tree with selection observables
          tue.gen_nj = jets.size();
          tue.gen_nb = nb;
          int posLepton( leptons[0].charge()>0 ? 0 : 1 );
          tue.gen_mll    = mll;
          tue.gen_ptpos  = leptons[posLepton].pt();
          tue.gen_phipos = leptons[posLepton].phi();
          tue.gen_ptll   = dil.Pt();
          tue.gen_phill  = dil.Phi();
          tue.gen_sumpt  = leptons[0].pt()+leptons[1].pt();
          tue.gen_dphill = TMath::Abs(leptons[0].p4().DeltaPhi(leptons[1].p4()));          
          
          /////////////////////////
          // GEN LEVEL ANALYSIS //
          ///////////////////////
          
          // Collect charged particles into PseudoJet
          std::vector<fastjet::PseudoJet> fj_particles;
          for(int ipf = 0; ipf < ev.ngpf; ipf++) {
            if(ev.gpf_c[ipf]==0) continue;

            TLorentzVector tkP4(0,0,0,0);
            tkP4.SetPtEtaPhiM(ev.gpf_pt[ipf],ev.gpf_eta[ipf],ev.gpf_phi[ipf],0.);

            //fiducial cuts
            bool passKin(ev.gpf_pt[ipf]>0.9 && fabs(ev.gpf_eta[ipf])<2.5);
            if (not passKin) continue;
         
            //matching to leptons
            Double_t relDpt2lep(9999.);
            for(int ilep=0; ilep<2; ilep++) {
              float dR=tkP4.DeltaR(leptons[ilep].p4());
              if(dR>0.05) continue;
              float relDpt=fabs(leptons[ilep].pt()-tkP4.Pt())/leptons[ilep].pt();
              if(relDpt>relDpt2lep) continue;
              relDpt2lep=relDpt;
            }
            bool matchedToLepton(relDpt2lep<0.05);
            if (matchedToLepton) continue;
            
            fj_particles.push_back( fastjet::PseudoJet( tkP4.Px(),tkP4.Py(), tkP4.Pz(), tkP4.E() ) );
          }
          
          // TODO: calculate NJettiness at gen level

        }
      
      //proceed only if event is selected on gen or reco level
      if (tue.gen_sel + tue.reco_sel == -2) continue;

      //finalize ntuple
      tue.cat=0;
      if(chTag=="MM") tue.cat=13*13;
      if(chTag=="EM") tue.cat=11*13;
      if(chTag=="EE") tue.cat=11*11;
      tue.run=ev.run;
      tue.event=ev.event;
      tue.lumi=ev.lumi;
      tue.nvtx=ev.nvtx;
   
      //all done
      if(tue.cat==11*13) outT->Fill();
    }
  f->Close();
  //save histos to file  
  fOut->cd();
  cout << "Selected " << outT->GetEntriesFast() << " events, saving to output" << endl;
  outT->Write();
  for (auto& it : ht.getPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  for (auto& it : ht.get2dPlots())  { 
    it.second->SetDirectory(fOut); it.second->Write(); 
  }
  cout << "Histograms were saved" << endl;
  fOut->Close();
}



//
void createanalysisXConeTree(TTree *t,analysisXCone_t &tue)
{
  //event category
  t->Branch("run",         &tue.run,         "run/I");
  t->Branch("event",       &tue.event,       "event/I");
  t->Branch("lumi",        &tue.lumi,        "lumi/I");
  t->Branch("cat",         &tue.cat,         "cat/I");
  t->Branch("gen_passSel", &tue.gen_passSel, "gen_passSel/I");
  t->Branch("passSel",     &tue.passSel,     "passSel/I");
  t->Branch("gen_nj",      &tue.gen_nj,      "gen_nj/I");
  t->Branch("gen_nb",      &tue.gen_nb,      "gen_nb/I");
  t->Branch("nj",          &tue.nj,          "nj/I");
  t->Branch("nb",          &tue.nb,          "nb/I");
  t->Branch("nvtx",        &tue.nvtx,        "nvtx/I");
  
  //NJettiness
  //TODO: consistent names of variables and branches ;)
  t->Branch("N_cut",       &tue.NJettiness,  "N_cut/I");
  t->Branch("N_slope",     &tue.NJettSlope,  "N_slope/I");

  //event weights
  t->Branch("nw",     &tue.nw, "nw/I");
  t->Branch("weight",  tue.weight, "weight[nw]/F");

  //leptonic quantities
  t->Branch("ptpos",      &tue.ptpos ,      "ptpos/F");
  t->Branch("phipos",     &tue.phipos ,     "phipos/F");
  t->Branch("ptll",       &tue.ptll ,       "ptll/F");
  t->Branch("phill",      &tue.phill ,      "phill/F");
  t->Branch("mll",        &tue.mll ,        "mll/F");
  t->Branch("sumpt",      &tue.sumpt ,      "sumpt/F");
  t->Branch("dphill",     &tue.dphill ,     "dphill/F");
  t->Branch("gen_ptpos",  &tue.gen_ptpos ,  "gen_ptpos/F");
  t->Branch("gen_phipos", &tue.gen_phipos , "gen_phipos/F");
  t->Branch("gen_ptll",   &tue.gen_ptll ,   "gen_ptll/F");
  t->Branch("gen_phill",  &tue.gen_phill ,  "gen_phill/F");
  t->Branch("gen_mll",    &tue.gen_mll ,    "gen_mll/F");
  t->Branch("gen_sumpt",  &tue.gen_sumpt ,  "gen_sumpt/F");
  t->Branch("gen_dphill", &tue.gen_dphill , "gen_dphill/F");
}

//
void resetanalysisXCone(analysisXCone_t &tue)
{
  //dummy event header
  tue.run=-1;  tue.lumi=0;  tue.event=0;              
  
  //reset selection flags
  tue.cat=0;   
  tue.gen_passSel=0;       
  tue.passSel=0; 

  //reset weights
  tue.nw=0;      
  
  //reset all MC truth
  tue.gen_ptpos=0;
  tue.gen_phipos=0;
  tue.gen_ptll=0;
  tue.gen_mll=0;
  tue.gen_sumpt=0;
  tue.gen_dphill=0;
  tue.gen_nj=0;             
  tue.gen_nb=0;
  
  //reset selection flags
  tue.gen_sel = -1; tue.reco_sel = -1;
}
