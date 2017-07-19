#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

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
	      TString era,
	      Bool_t debug)
{
  bool isDataFile(filename.Contains("Data"));
  cout<< " We are in the analysis file"<< endl;
  //check file type from name
  if(isDataFile) runSysts=false;
  bool isTTbar( filename.Contains("_TTJets") );

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
  cout<<"Corrections"<<endl;
  //CORRECTIONS
  std::vector<RunPeriod_t> runPeriods=getRunPeriods(era);
  cout<<"lumi"<<endl;
  //lumi
  TH1F *ratevsrunH=0;
  std::map<Int_t,Float_t> lumiMap;
  if( isDataFile )  
    {
      std::pair<std::map<Int_t,Float_t>, TH1F *> result=parseLumiInfo(era);
      lumiMap   = result.first;
      ratevsrunH = result.second;
    }
	cout<<"pileup"<<endl;
  //pileup
  std::map<TString, std::vector<TGraph *> > puWgtGr;
  if( !isDataFile ) puWgtGr=getPileupWeightsMap(era,genPU);
  cout<<"lepton eff."<<endl;
  //lepton efficiencies
  LeptonEfficiencyWrapper lepEffH(isDataFile,era);
	cout<<"btag"<<endl;
  //b-tagging
  BTagSFUtil myBTagSFUtil;
  std::map<TString, std::map<BTagEntry::JetFlavor, BTagCalibrationReader *> > btvsfReaders = getBTVcalibrationReadersMap(era, BTagEntry::OP_MEDIUM);
  std::map<BTagEntry::JetFlavor, TGraphAsymmErrors *> expBtagEff = readExpectedBtagEff(era);
  cout<<"jet energy uncertai"<<endl;
  //JET ENERGY UNCERTAINTIES    
  TString jecUncUrl(era+"/Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");
  gSystem->ExpandPathName(jecUncUrl);
  JetCorrectorParameters *jecParam = new JetCorrectorParameters(jecUncUrl.Data(), "Total");
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty( *jecParam );
cout<<"prepare out"<<endl;
  //PREPARE OUTPUT
  analysisXCone_t tue;
  TFile *fOut=TFile::Open(outname,"RECREATE");
  fOut->cd();
  TTree *outT=new TTree("tue","tue");
  createanalysisXConeTree(outT,tue);
  outT->SetDirectory(fOut);
cout<<"book histograms"<<endl;
  //BOOK CONTROL HISTOGRAMS
  HistTool ht; 
  ht.setNsyst(isDataFile ? 0 : 3);
  ht.addHist("puwgtctr", new TH1F("puwgtctr","Weight sums",4,0,4) );  
  std::vector<TString> lfsVec = { "EE", "EM", "MM" };
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
      ht.addHist("nchvsnvtx_"+tag , new TH2F("nchvsnvtx_"+tag,";Vertex multiplicity;Charged particle multiplicity;Events",10,0,40,50,0,100) );
      ht.addHist("nchvsrho_"+tag  , new TH2F("nchvsrho_"+tag,";#rho;Charged particle multiplicity;Events",10,0,40,50,0,100) );      
      ht.addHist("nch_"+tag      , new TH1F("nch_"+tag,";Charged particle multiplicity;Events",50,0,200) );      
      ht.addHist("chavgpt_"+tag  , new TH1F("chavgpt_"+tag,";Charged particle average p_{T} [GeV];Events",50,0,15) );      
      ht.addHist("chsumpt_"+tag   , new TH1F("chsumpt_"+tag,";Charged particle sum p_{T} [GeV];Events",50,0,400) );
      ht.addHist("chavgpz_"+tag   , new TH1F("chavgpz_"+tag,";Charged particle average p_{z} [GeV];Events",50,0,15) );      
      ht.addHist("chsumpz_"+tag   , new TH1F("chsumpz_"+tag,";Charged particle sum p_{z} [GeV];Events",50,0,400) );
	ht.addHist("XConeCut_"+tag, new TH1F("XConeCut_"+tag,"; Jet XConeCut; Events",10,0,10));
	ht.addHist("XConeDiff_"+tag, new TH1F("XConeDiff_"+tag,"; Jet XConeDiff; Events",10,0,10));
	for (Int_t ii = 2;  ii<10; ii++){
		ht.addHist(ii+"_Jettiness"+tag, new TH1F(ii+"_Jettiness"+tag,";tau ;Events",50,0,200));
		ht.addHist(ii+"_Jettiness_slope"+tag, new TH1F(ii+"_Jettiness_slope"+tag,";d tau/d N ;Events",20,-0.1,0));
      }
    }
  for (auto& it : ht.getPlots() )     { it.second->Sumw2(); it.second->SetDirectory(0); }
  for (auto& it : ht.get2dPlots() )   { it.second->Sumw2(); it.second->SetDirectory(0); }
cout<<"Njettiness"<<endl;
  //NJettiness variables 
  double tau_cut = 30.0 ;
  double beta = 2.;
  double R = 0.4;
 

  //LOOP OVER EVENTS
  for (Int_t iev=0;iev<nentries;iev++)
    {
	  
	   double tau_temp_C = 1000.;
       double tau_temp_S = 1000.0;
       unsigned int N_temp_S = 1;
       unsigned int N_temp_C = 1;
      TString period("");
      t->GetEntry(iev);
      resetanalysisXCone(tue);
      if(iev%1000==0) printf("\r [%3.0f/100] done",100.*(float)(iev)/(float)(nentries));
      //assign a run period and correct the event accordingly
      float puWgt(1.0),puWgtUp(1.0),puWgtDn(1.0),topptsf(1.0);
      if(!ev.isData)
	{
	  period=assignRunPeriod(runPeriods);
	  puWgt   = puWgtGr[period][0]->Eval(ev.g_pu);
	  puWgtUp = puWgtGr[period][1]->Eval(ev.g_pu);
	  puWgtDn = puWgtGr[period][2]->Eval(ev.g_pu);

	  //top pt weighting
	  if(isTTbar)
	    {
	      for(Int_t igen=0; igen<ev.ngtop; igen++)
		{
		  if(abs(ev.gtop_id[igen])!=6) continue;
		  topptsf *= TMath::Exp(0.156-0.00137*ev.gtop_pt[igen]);
		}
	    }
	}

      //the main objects of the analysis
      std::vector<Particle> recoTracks, genTracks;

      //,
      
      //RECO LEVEL analysis
      //

      //selection
      //evsel.setDebug(true);
      TString chTag = evsel.flagFinalState(ev);
      if(chTag=="EM" )
	{
	  //leptons
	  std::vector<Particle> &leptons=evsel.getSelLeptons();

	  //divide jets
	  std::vector<Jet> jets=evsel.getGoodJets(ev,25.,2.4,leptons);
	  sort(jets.begin(),jets.end(),Jet::sortJetsByCSV);

	  //select the charged PF candidates 
	  //veto neutrals
	  //veto if associated to the two leading b-jets		
	    
	  std::vector<fastjet::PseudoJet> fj_particles;
	  
	  std::vector<std::pair<int,TLorentzVector> > selTracks,hpTracks;

	  for(int ipf = 0; ipf < ev.npf; ipf++) 
	    {
	      if(ev.pf_c[ipf]==0) continue;

	      TLorentzVector tkP4(0,0,0,0);
	      tkP4.SetPtEtaPhiM(ev.pf_pt[ipf],ev.pf_eta[ipf],ev.pf_phi[ipf],0.);

	      //fiducial cuts
	      bool passKin(ev.pf_pt[ipf]>0.9 && fabs(ev.pf_eta[ipf])<2.5);
	   
	      //matching to leptons
	      Double_t relDpt2lep(9999.);
	      for(int ilep=0; ilep<2; ilep++)
		{
		  float dR=tkP4.DeltaR(leptons[ilep].p4());
		  if(dR>0.05) continue;
		  float relDpt=fabs(leptons[ilep].pt()-tkP4.Pt())/leptons[ilep].pt();
		  if(relDpt>relDpt2lep) continue;
		  relDpt2lep=relDpt;
		}
	      bool matchedToLepton(relDpt2lep<0.05);
		  if (matchedToLepton) continue;
			fj_particles.push_back( fastjet::PseudoJet( tkP4.Px(),tkP4.Py(), tkP4.Pz(), tkP4.E() ) );
			
	      //matching to leading CSV jet candidates
	      bool clusteredInBjet(false);
	      for(size_t ibj=0; ibj<min(jets.size(),size_t(2)); ibj++)
		{
		  std::vector<Particle> &pinJet=jets[ibj].particles();
		  for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
		    {
		      if(pinJet[ipinj].charge()==0) continue;
		      if(pinJet[ipinj].originalReference()!=ipf) continue;
		      clusteredInBjet=true;
		      break;
		    }
		  if(clusteredInBjet) break;
		}

	      recoTracks.push_back( Particle(tkP4, ev.pf_c[ipf], ev.pf_id[ipf],
					     (passKin | matchedToLepton <<1 | clusteredInBjet <<2),
					     ipf,
					     1) 
				    );
	    }

	
	  //save PF cands
	  tue.n=0;
	  float nch(0.),chSumPt(0.),chSumPz(0.);
	  for(auto p : recoTracks)
	    {
	      //check if it's not matched to hard process and charge is non-null
	      if(p.qualityFlags()!=1 || p.charge()==0)
		{
		  p.setOriginalReference(-1);
		  continue;
		}
	      p.setOriginalReference(tue.n); //set to the reference in the ntuple
	      tue.pt[tue.n]  = p.pt();
	      tue.eta[tue.n] = p.eta();
	      tue.phi[tue.n] = p.phi();
	      tue.id[tue.n]  = p.id();
	      tue.n++;
	      nch++;
	      chSumPt+=p.pt();
	      chSumPz+=fabs(p.pz());
	    }

	  //event weight
	  float wgt(1.0);
	  ht.getPlots()["puwgtctr"]->Fill(0.,1.0);
	  EffCorrection_t lepselSF(1.0,0.0),trigSF(1.0,0.0);
	  if(!ev.isData) 
	    {
	      ht.getPlots()["puwgtctr"]->Fill(1,puWgt);	      
	      ht.getPlots()["puwgtctr"]->Fill(2,puWgtUp);	      
	      ht.getPlots()["puwgtctr"]->Fill(3,puWgtDn);	      
	      for(size_t il=0; il<2; il++)
		{
                  EffCorrection_t sf=lepEffH.getOfflineCorrection(leptons[il].id(),leptons[il].pt(),leptons[il].eta(),period);
                  lepselSF.second = sqrt( pow(lepselSF.first*sf.second,2)+pow(lepselSF.second*sf.first,2));
                  lepselSF.first *= sf.first;
                }
	      trigSF=lepEffH.getTriggerCorrection(leptons,period);

	      wgt  = (normH? normH->GetBinContent(1) : 1.0);
	      wgt *= puWgt;
	      wgt *= trigSF.first;
	      wgt *= lepselSF.first;
	      wgt *= (ev.g_nw>0 ? ev.g_w[0] : 1.0);
	    }

	  //weights
	  tue.nw=1;
	  tue.weight[0]=wgt;
	  if(runSysts && isTTbar && ev.g_nw>1)
	    {
	      //pu{up,down}
	      tue.weight[1]=puWgt!=0 ? wgt*puWgtUp/puWgt : wgt;
	      tue.weight[2]=puWgt!=0 ? wgt*puWgtDn/puWgt : wgt;

	      //eff{up,down}
	      float effCen=trigSF.first*lepselSF.first;
	      float effUp=(trigSF.first+trigSF.second)*(lepselSF.first+lepselSF.second);
	      tue.weight[3]=wgt*effUp/effCen;
	      float effDn=(trigSF.first-trigSF.second)*(lepselSF.first-lepselSF.second);
	      tue.weight[4]=wgt*effDn/effCen;

	      //top pt
	      tue.weight[5]=wgt*topptsf;

	      //generator level weights
	      for(size_t iw=1; iw<=120; iw++)
		tue.weight[5+iw]= ev.g_w[0]!=0 ? wgt* ev.g_w[iw]/ev.g_w[0] : wgt;     

	      tue.nw=126;
	    }

	  std::vector<double>plotwgts(1,wgt);
	  if(!ev.isData)
	    {
	      plotwgts.push_back(puWgt!=0 ? puWgtUp/puWgt : 1.0);
	      plotwgts.push_back(puWgt!=0 ? puWgtDn/puWgt : 1.0);
	    }
	  

	  //check the selection for different experimental variations
	  //nom,b-tag{up,down},jes{up,down},jer{up,down},ees{up,down},mes{up,down}
	  tue.passSel=0;
	  std::vector<Bool_t> nomBtag;
	  for(size_t ivar=0; ivar<=(runSysts ? 10 : 0); ivar++)
	    {
	      //leptons
	      TLorentzVector l1(leptons[0].p4()), l2(leptons[1].p4());
	      float l1scaleUnc(ev.l_scaleUnc[leptons[0].originalReference()]/l1.E());
	      float l2scaleUnc(ev.l_scaleUnc[leptons[1].originalReference()]/l2.E());
	      if(ivar==7 || ivar==8 || ivar==9 || ivar==10)
		{
		  float sign((ivar==7||ivar==8) ? +1.0 : -1.0);
		  if(ivar==7 || ivar==9)
		    {
		      if( abs(leptons[0].id())==13 ) l1 *= (1+sign*l1scaleUnc);
		      else                           l2 *= (1+sign*l2scaleUnc);
		    }
		  if(ivar==8 || ivar==10)
		    {
		      if( abs(leptons[0].id())==11 ) l1 *= (1+sign*l1scaleUnc);
		      else                           l2 *= (1+sign*l2scaleUnc);
		    }
		}
	      TLorentzVector dil(l1+l2);
	      float mll=dil.M();
	      bool passLepPresel(mll>12
				 && (l1.Pt()>25 || l2.Pt()>25)
				 && (fabs(l1.Eta())<2.5 && fabs(l2.Eta())<2.5) );
	      
	      //jets
	      std::vector< std::pair<TLorentzVector,int> > lightJets,bJets;
	      for(size_t ij=0; ij<jets.size(); ij++)
		{
		  TLorentzVector jp4(jets[ij].p4());
		  int k=jets[ij].getJetIndex();
		  bool isBTagged( ev.j_csv[k]>0.8484);

		  if(!ev.isData) 
		    {
		      //JER
		      float genJet_pt(0);
		      if(ev.j_g[k]>-1) genJet_pt = ev.g_pt[ ev.j_g[k] ];
		      if(genJet_pt>0) {
			std::string option("");
			if(ivar==5) option="up";
			if(ivar==6) option="down";
			smearJetEnergy(jp4,genJet_pt,option);
		      }

		      //JES
		      if(ivar==3 || ivar==4)
			{
			  std::string option( ivar==3 ? "up" : "down" );
			  applyJetCorrectionUncertainty(jp4,jecUnc,option);
			}

		      //b-tagging
		      if(ivar<3)
			{
			  float jptForBtag(jp4.Pt()>1000. ? 999. : jp4.Pt()), jetaForBtag(fabs(jp4.Eta()));
			  BTagEntry::JetFlavor hadFlav=BTagEntry::FLAV_UDSG;
			  std::string option("central");
			  if(ivar==1) option="up";
			  if(ivar==2) option="down";
			  if(abs(ev.j_hadflav[k])==4) { hadFlav=BTagEntry::FLAV_C; option = "central"; }
			  if(abs(ev.j_hadflav[k])==5) { hadFlav=BTagEntry::FLAV_B; option = "central"; }
			  float expEff = expBtagEff[hadFlav]->Eval(jptForBtag);
			  float jetBtagSF = btvsfReaders[period][hadFlav]->eval_auto_bounds( option, hadFlav, jetaForBtag, jptForBtag);
			  myBTagSFUtil.modifyBTagsWithSF(isBTagged, jetBtagSF, expEff);
			}
		    }

		  //select jet
		  if(ivar==0) nomBtag.push_back(isBTagged);
		  if(ivar>=3) isBTagged=nomBtag[ij];
		  if(jp4.Pt()<30 || fabs(jp4.Eta())>2.5) continue;
		  if(isBTagged) bJets.push_back( std::pair<TLorentzVector,bool>(jp4,ij) );
		  else          lightJets.push_back( std::pair<TLorentzVector,bool>(jp4,ij) );
		}
	      
	      //require that the two leading CSV are b-tagged
	      int nb(0),nj(bJets.size()+lightJets.size());
	      if(bJets.size()>0) 
		{
		  if(bJets[0].second==0) nb++;
		  if(bJets.size()>1) if(bJets[1].second==1) nb++;
		}

	      nj-=nb;
	      	      
	      bool passPresel = passLepPresel && nb==2;
	      if(chTag=="EE" || chTag=="MM") passPresel &= fabs(mll-91)>15;
	      
	      //flag if passes selection
	      tue.passSel |= (passPresel << ivar);
	      tue.nj[ivar]=nj;
	      tue.nb[ivar]=nb;

	      TLorentzVector rec_tt(dil);
	      if(nb==2) rec_tt += bJets[0].first+bJets[1].first;
	      rec_tt += evsel.getMET();
	      tue.ptttbar[ivar]=rec_tt.Pt();
	      tue.phittbar[ivar]=rec_tt.Phi();	  
      
	      tue.mll[ivar]    = mll;
	      tue.ptpos[ivar]  = leptons[0].charge()>0 ? l1.Pt() : l2.Pt();
	      tue.phipos[ivar] = leptons[0].charge()>0 ? l1.Phi() : l2.Phi();
	      tue.ptll[ivar]   = dil.Pt();
	      tue.phill[ivar]  = dil.Phi();
	      tue.sumpt[ivar]  = l1.Pt()+l2.Pt();
	      tue.dphill[ivar] = TMath::Abs(l1.DeltaPhi(l2));

		
	      cout << N_temp_C << ", " << int(N_temp_C) << ", " << tue.NJettiness << std::endl;
	      //histogram filling for the nominal selection
	      if(ivar!=0) continue;
	      if(!passLepPresel) continue;
	      
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

	      ht.fill("nvtx_"+chTag,ev.nvtx,plotwgts);
	      ht.fill("rho_"+chTag,ev.rho,plotwgts);
	      tue.NJettiness= int( N_temp_C);
		tue.NJettSlope = int (N_temp_S);
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
	      if(!passPresel) continue;
	      std::map<Int_t,Float_t>::iterator rIt=lumiMap.find(ev.run);
              if(rIt!=lumiMap.end() && ratevsrunH) ht.getPlots()["ratevsrun_"+chTag]->Fill(std::distance(lumiMap.begin(),rIt),1./rIt->second);
              ht.fill("met_"+chTag,evsel.getMET().Pt(),plotwgts);	      
              ht.fill("ptttbar_"+chTag,tue.ptttbar[0],plotwgts);
	      ht.fill("njets_"+chTag,nj,plotwgts);
	      ht.fill("sumpt_"+chTag,tue.sumpt[0],plotwgts);	     
	      ht.fill("ptpos_"+chTag,tue.ptpos[0],plotwgts);
	      ht.fill("ptll_"+chTag,tue.ptll[0],plotwgts);	      
	      ht.fill("nch_"+chTag,nch,plotwgts);	  
	      ht.get2dPlots()["nchvsnvtx_"+chTag]->Fill(ev.nvtx,nch,wgt);
	      ht.get2dPlots()["nchvsrho_"+chTag]->Fill(ev.rho,nch,wgt);
	      ht.fill("chavgpt_"+chTag,nch>0 ? chSumPt/nch : -1,plotwgts);
	      ht.fill("chsumpt_"+chTag,chSumPt,plotwgts);
	      ht.fill("chavgpz_"+chTag,nch>0 ? chSumPz/nch : -1,plotwgts);
	      ht.fill("chsumpz_"+chTag,chSumPz,plotwgts);
		ht.fill("XConeCut_"+chTag,int(N_temp_C),plotwgts);
	      ht.fill("XConeDiff_"+chTag,int(N_temp_S),plotwgts);
	    
	    }
	}

      //
      // GEN LEVEL ANALYSIS
      //
      TString genChTag = evsel.flagGenFinalState(ev);
      if(!ev.isData && (genChTag=="EM" || genChTag=="EE" || genChTag=="MM"))
	{
	  std::vector<Particle> &leptons=evsel.getGenLeptons();
	  TLorentzVector dil(leptons[0].p4()+leptons[1].p4());
	  float mll=dil.M();
	  bool passLepPresel(mll>12
                             && (leptons[0].pt()>25 || leptons[1].pt()>25)
                             && (fabs(leptons[0].eta())<2.5 && fabs(leptons[1].eta())<2.5) );

	  //divide jets                                                                                                                              
	  std::vector<Jet> &jets=evsel.getGenJets() ;
	  std::vector<size_t> lightJetsIdx, bJetsIdx;
          for(size_t ij=0; ij<jets.size(); ij++)
            {
              if(abs(jets[ij].flavor())==5 && bJetsIdx.size()<2) bJetsIdx.push_back(ij);
              else                          lightJetsIdx.push_back(ij);
            }
          bool passPresel = passLepPresel && (bJetsIdx.size()>=2);
          if(genChTag=="EE" || genChTag=="MM") passPresel &= fabs(mll-91)>15;

	  //track selection
	  for(int ipf = 0; ipf < ev.ngpf; ipf++) 
	    {
	      if(ev.gpf_c[ipf]==0) continue;

	      TLorentzVector tkP4(0,0,0,0);
	      tkP4.SetPtEtaPhiM(ev.gpf_pt[ipf],ev.gpf_eta[ipf],ev.gpf_phi[ipf],0.);
	      
	      //fiducial cuts
              bool passKin(ev.gpf_pt[ipf]>0.9 && fabs(ev.gpf_eta[ipf])<2.5);
	      
	      //check if matching the leptons
	      bool matchedToLepton(false);
	      for(int ilep=0; ilep<2; ilep++)
		{
		  float dR=tkP4.DeltaR(leptons[ilep].p4());
		  if(dR<0.1 && leptons[ilep].id()==ev.gpf_id[ipf]) matchedToLepton=true;
		}

	      //matching to b-jet candidates
              bool clusteredInBjet(false);
	      for(size_t ibj=0; ibj<min(bJetsIdx.size(),size_t(2)); ibj++)
		{
		  std::vector<Particle> &pinJet=jets[ bJetsIdx[ibj] ].particles();
		  for(size_t ipinj=0; ipinj<pinJet.size(); ipinj++)
		    {
		      if(pinJet[ipinj].charge()==0) continue;
		      if(pinJet[ipinj].originalReference()!=ipf) continue;
		      clusteredInBjet=true;
		    }
		}
		
	      genTracks.push_back( Particle(tkP4,ev.gpf_c[ipf],ev.gpf_id[ipf],
					    (passKin | matchedToLepton <<1 | clusteredInBjet <<2),
					    ipf,
					    1)
				   );
	    }
	  
	  //save gen candidates
	  tue.gen_n=0;
	  for(auto p : genTracks)
	    {
	      if(p.qualityFlags()!=1) continue;
	      tue.gen_pt[tue.gen_n]  = p.pt();
	      tue.gen_eta[tue.gen_n] = p.eta();
	      tue.gen_phi[tue.gen_n] = p.phi();
	      tue.gen_id[tue.gen_n]  = p.id();

	      //match to reco level
	      Particle *recoMatch=0;
	      for(auto r : recoTracks)
		{
		  float dR=p.p4().DeltaR(r.p4());
		  if(dR>0.1) continue;
		  if(recoMatch==0) recoMatch=&p;
		  else
		    {
		      float dpt=fabs(p.pt()-r.pt());
		      float prev_dpt=fabs(p.pt()-recoMatch->pt());
		      if(dpt>prev_dpt) continue;
		      recoMatch=&p;
		    }	      
		}
	      tue.gen_rec[tue.gen_n] = recoMatch!=0 ? recoMatch->originalReference() : -1;

	      //check the unmatched, high pT
	      // if(passPresel && tue.passSel && p.pt()>10)
	      // 	{
	      // 	  if(recoMatch!=0 && recoMatch->charge()==0)
	      // 	    {
	      // 	      cout << p.pt() << " " << p.id() << " " << p.qualityFlags() << " " << p.eta() << " | "
	      // 		   << "\t" << recoMatch->pt() << " " << recoMatch->id() << " " << recoMatch->qualityFlags() << endl;
	      // 	    }
	      // 	  else if(recoMatch==0)
	      // 	    {
	      // 	      cout << p.pt() << " " << p.id() << " " << p.qualityFlags() << " " << p.eta() << endl;
	      // 	    }
	      // 	}

	      tue.gen_n++;
	    }


	  //flag if passes selection  
	  tue.gen_passSel=passPresel;
	  tue.gen_nj=jets.size()-bJetsIdx.size();
          tue.gen_nb=bJetsIdx.size();
	  
	  int posLepton( leptons[0].charge()>0 ? 0 : 1 );
          tue.gen_mll    = mll;
          tue.gen_ptpos  = leptons[posLepton].pt();
          tue.gen_phipos = leptons[posLepton].phi();
          tue.gen_ptll   = dil.Pt();
          tue.gen_phill  = dil.Phi();
          tue.gen_sumpt  = leptons[0].pt()+leptons[1].pt();
          tue.gen_dphill = TMath::Abs(leptons[0].p4().DeltaPhi(leptons[1].p4()));

	  //save ttbar and pseudo-ttbar kinematics
	  if(ev.ngtop>0)
	    {
	      TLorentzVector ttbar(0,0,0,0);
	      Int_t nfs(0);
	      for(Int_t it=0; it<ev.ngtop; it++)
		{
		  int absid(abs(ev.gtop_id[it]));
		  if(absid!=5000 && absid!=11000 && absid!=13000 && absid!=0) continue;
		  nfs++;
		  TLorentzVector p4(0,0,0,0);
		  p4.SetPtEtaPhiM(ev.gtop_pt[it],ev.gtop_eta[it],ev.gtop_phi[it],ev.gtop_m[it]);
		  ttbar += p4;
		}
	      if(nfs!=5) ttbar.SetPtEtaPhiM(0,0,0,0);
	      tue.gen_ptttbar=ttbar.Pt();
	      tue.gen_phittbar=ttbar.Phi();
	    }
	}
      
      //check it if it passed at least one of b-tagging selections
      if(tue.passSel==0 && tue.gen_passSel==0) continue;

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
  t->Branch("nj",           tue.nj,          "nj[11]/I");
  t->Branch("nb",           tue.nb,          "nb[11]/I");
  t->Branch("nvtx",        &tue.nvtx,        "nvtx/I");
  t->Branch("N_cut",	 &tue.NJettiness,     "N_cut/I");
  t->Branch("N_slope",		 &tue.NJettSlope, "N_slope/I");

  //event weights
  t->Branch("nw",     &tue.nw, "nw/I");
  t->Branch("weight",  tue.weight, "weight[nw]/F");

  //ptttbar
  t->Branch("gen_ptttbar",         &tue.gen_ptttbar,      "gen_pttbar/F");
  t->Branch("gen_phittbar",        &tue.gen_phittbar,     "gen_phittbar/F");
  t->Branch("ptttbar",              tue.ptttbar,          "ptttbar[11]/F");
  t->Branch("phittbar",             tue.phittbar,         "phittbar[11]/F");

  //leptonic quantities
  t->Branch("ptpos",      tue.ptpos ,       "ptpos[11]/F");
  t->Branch("phipos",     tue.phipos ,      "phipos[11]/F");
  t->Branch("ptll",       tue.ptll ,        "ptll[11]/F");
  t->Branch("phill",      tue.phill ,       "phill[11]/F");
  t->Branch("mll",        tue.mll ,         "mll[11]/F");
  t->Branch("sumpt",      tue.sumpt ,       "sumpt[11]/F");
  t->Branch("dphill",     tue.dphill ,      "dphill[11]/F");
  t->Branch("gen_ptpos",  &tue.gen_ptpos ,  "gen_ptpos/F");
  t->Branch("gen_phipos", &tue.gen_phipos , "gen_phipos/F");
  t->Branch("gen_ptll",   &tue.gen_ptll ,   "gen_ptll/F");
  t->Branch("gen_phill",  &tue.gen_phill ,  "gen_phill/F");
  t->Branch("gen_mll",    &tue.gen_mll ,    "gen_mll/F");
  t->Branch("gen_sumpt",  &tue.gen_sumpt ,  "gen_sumpt/F");
  t->Branch("gen_dphill", &tue.gen_dphill , "gen_dphill/F");

  //charged particles
  t->Branch("n",          &tue.n,            "n/I");
  t->Branch("id",          tue.id ,          "id[n]/I");
  t->Branch("pt",          tue.pt ,          "pt[n]/F");
  t->Branch("eta",         tue.eta ,         "eta[n]/F");
  t->Branch("phi",         tue.phi ,         "phi[n]/F");

  //gen charged particles
  t->Branch("gen_n",   &tue.gen_n,     "gen_n/I");
  t->Branch("gen_id",   tue.gen_id ,   "gen_id[gen_n]/I");
  t->Branch("gen_rec",  tue.gen_rec ,  "gen_rec[gen_n]/I");
  t->Branch("gen_pt",   tue.gen_pt ,   "gen_pt[gen_n]/F");
  t->Branch("gen_eta",  tue.gen_eta ,  "gen_eta[gen_n]/F");
  t->Branch("gen_phi",  tue.gen_phi ,  "gen_phi[gen_n]/F");
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
  tue.gen_ptttbar=0; 
  tue.gen_phittbar=0;
  tue.gen_ptpos=0;
  tue.gen_phipos=0;
  tue.gen_ptll=0;
  tue.gen_mll=0;
  tue.gen_sumpt=0;
  tue.gen_dphill=0;
  tue.gen_nj=0;             
  tue.gen_nb=0;
  
  //reset particle counters
  tue.n=0; 
  tue.gen_n=0; 
}
