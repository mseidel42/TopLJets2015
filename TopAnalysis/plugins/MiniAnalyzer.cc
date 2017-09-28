//
// -*- C++ -*-
//
// Package:    TopLJets2015/TopAnalysis
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TopLJets2015/TopAnalysis/interface/MiniEvent.h"
#include "TopLJets2015/TopAnalysis/interface/MyIPTools.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  virtual void endRun(const edm::Run&,const edm::EventSetup&);  
private:
  virtual void beginJob() override;
  int genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  int recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  float getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			 const reco::Candidate* ptcl,  
                         float r_iso_min, float r_iso_max, float kt_scale,
                         bool charged_only);

  bool isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex);
  bool isMediumMuon2016ReReco(const reco::Muon & recoMu);

  // member data 
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorevtToken_;
  edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
  edm::EDGetTokenT<LHERunInfoProduct> generatorRunInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>  > genLeptonsToken_,   genJetsToken_;
  edm::EDGetTokenT<reco::METCollection> genMetsToken_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticlesToken_;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_,metFilterBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<edm::View<pat::Electron>  >  electronToken_;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_, puppiMetToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  
  //Electron Decisions
  edm::EDGetTokenT<edm::ValueMap<float> > eleMvaIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_,eleLooseIdMapToken_,eleMediumIdMapToken_,eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleVetoIdFullInfoMapToken_,eleLooseIdFullInfoMapToken_,eleMediumIdFullInfoMapToken_,eleTightIdFullInfoMapToken_;
  unsigned int evetoIsoBit_, elooseIsoBit_, emediumIsoBit_, etightIsoBit_;
  edm::EDGetTokenT<bool> BadChCandFilterToken_,BadPFMuonFilterToken_;

  //  edm::EDGetTokenT<edm::ValueMap<float> > petersonFragToken_;

  std::unordered_map<std::string,TH1*> histContainer_;

  PFJetIDSelectionFunctor pfjetIDLoose_;

  std::vector<std::string> triggersToUse_,metFiltersToUse_;

  bool saveTree_,savePF_,runOnGEN_;
  TTree *tree_;
  MiniEvent_t ev_;
  
  bool useRawLeptons_;
  edm::Service<TFileService> fs;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
  generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
  generatorevtToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator",""))),
  generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
  generatorRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>({"externalLHEProducer"})),
  puToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),  
  genLeptonsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:leptons"))),
  genJetsToken_(consumes<std::vector<reco::GenJet> >(edm::InputTag("particleLevel:jets"))),
  genMetsToken_(consumes<reco::METCollection>(edm::InputTag("particleLevel:mets"))),
  genParticlesToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
  prunedGenParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  metFilterBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  jetToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),  
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  puppiMetToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("puppimets"))),
  pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  eleMvaIdMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("eleMvaIdMap"))),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  eleVetoIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdFullInfoMap"))),
  eleLooseIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdFullInfoMap"))),
  eleMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdFullInfoMap"))),
  eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
  evetoIsoBit_(999), elooseIsoBit_(999), emediumIsoBit_(999), etightIsoBit_(999),
  BadChCandFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandFilter"))),
  BadPFMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"))),
  //petersonFragToken_(consumes<edm::ValueMap<float> >(edm::InputTag("bfragWgtProducer:PetersonFrag"))),
  pfjetIDLoose_( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE ),  
  saveTree_( iConfig.getParameter<bool>("saveTree") ),
  savePF_( iConfig.getParameter<bool>("savePF") ),
  runOnGEN_( iConfig.getParameter<bool>("runOnGEN") ),
  useRawLeptons_( iConfig.getParameter<bool>("useRawLeptons") )
{
  //now do what ever initialization is needed
  electronToken_      = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"));
  triggersToUse_      = iConfig.getParameter<std::vector<std::string> >("triggersToUse");
  metFiltersToUse_  = iConfig.getParameter<std::vector<std::string> >("metFiltersToUse");

  histContainer_["triggerList"] = fs->make<TH1F>("triggerList", ";Trigger bits;",triggersToUse_.size(),0,triggersToUse_.size());
  for(size_t i=0; i<triggersToUse_.size(); i++) histContainer_["triggerList"] ->GetXaxis()->SetBinLabel(i+1,triggersToUse_[i].c_str());
  histContainer_["counter"]    = fs->make<TH1F>("counter", ";Counter;Events",2,0,2);
  histContainer_["fidcounter"] = (TH1 *)fs->make<TH2F>("fidcounter",    ";Variation;Events", 1000, 0., 1000.,11,0,11); 
  histContainer_["pu"]         = fs->make<TH1F>("pu",      ";Pileup observed;Events / 1",100,0,100);
  histContainer_["putrue"]     = fs->make<TH1F>("putrue",  ";Pileup true;Events / 0.1",100,0,100);
  for(std::unordered_map<std::string,TH1*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();

  //create a tree for the selected events
  if(saveTree_)
    {
      tree_ = fs->make<TTree>("data","data");
      createMiniEventTree(tree_,ev_);
    }
}


//
MiniAnalyzer::~MiniAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//
int MiniAnalyzer::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //
  // PILEUP
  //
  if(!runOnGEN_) {
    edm::Handle<std::vector <PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(puToken_,PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator ipu;
    for (ipu = PupInfo->begin(); ipu != PupInfo->end(); ++ipu)
      {
        if ( ipu->getBunchCrossing() != 0 ) continue; // storing detailed PU info only for BX=0
        ev_.g_pu=ipu->getPU_NumInteractions();
        ev_.g_putrue=ipu->getTrueNumInteractions();
      }
    histContainer_["pu"]->Fill(ev_.g_pu);
    histContainer_["putrue"]->Fill(ev_.g_putrue);
  }
  
  //
  // GENERATOR WEIGHTS
  //
  ev_.g_nw=0; ev_.g_w[0]=1.0;
  edm::Handle<GenEventInfoProduct> evt;
  iEvent.getByToken( generatorToken_,evt);
  if(evt.isValid())
    {
      ev_.g_w[0] = evt->weight();
      ev_.g_nw++;

      //PDF info
      ev_.g_qscale = evt->pdf()->scalePDF;
      ev_.g_x1     = evt->pdf()->x.first;
      ev_.g_x2     = evt->pdf()->x.second;
      ev_.g_id1    = evt->pdf()->id.first;
      ev_.g_id2    = evt->pdf()->id.second;
    }
  histContainer_["counter"]->Fill(1,ev_.g_w[0]);
  
  //alternative weights for systematics 
  edm::Handle<LHEEventProduct> evet;
  iEvent.getByToken(generatorlheToken_, evet);
  if(evet.isValid())
    {
      double asdd=evet->originalXWGTUP();
      for(unsigned int i=0  ; i<evet->weights().size();i++){
	double asdde=evet->weights()[i].wgt;
	ev_.g_w[ev_.g_nw]=ev_.g_w[0]*asdde/asdd;
	ev_.g_nw++;
      }
    }
     
  //
  // GENERATOR LEVEL EVENT
  //
  ev_.ng=0; 
  edm::Handle<std::vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken_,genJets);  
  std::map<const reco::Candidate *,int> jetConstsMap;
  //edm::Handle<edm::ValueMap<float> > petersonFrag;
  //iEvent.getByToken(petersonFragToken_,petersonFrag);
  int ngjets(0),ngbjets(0);
  for(auto genJet=genJets->begin(); genJet!=genJets->end(); ++genJet)
    {

      edm::Ref<std::vector<reco::GenJet> > genJetRef(genJets,genJet-genJets->begin());

      //map the gen particles which are clustered in this jet
      std::vector< const reco::Candidate * > jconst=genJet->getJetConstituentsQuick();
      for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;

      ev_.g_id[ev_.ng]   = genJet->pdgId();
      ev_.g_pt[ev_.ng]   = genJet->pt();
      ev_.g_eta[ev_.ng]  = genJet->eta();
      ev_.g_phi[ev_.ng]  = genJet->phi();
      ev_.g_m[ev_.ng]    = genJet->mass();       
      ev_.ng++;
      
      //gen level selection
      if(genJet->pt()>25 && fabs(genJet->eta())<2.5)
	{
	  ngjets++;	
	  if(abs(genJet->pdgId())==5) ngbjets++;
	}
    }

  //leptons
  int ngleptons(0);
  edm::Handle<std::vector<reco::GenJet> > dressedLeptons;  
  iEvent.getByToken(genLeptonsToken_,dressedLeptons);
  for(auto genLep = dressedLeptons->begin();  genLep != dressedLeptons->end(); ++genLep)
    {
      //map the gen particles which are clustered in this lepton
      std::vector< const reco::Candidate * > jconst=genLep->getJetConstituentsQuick();
      for(size_t ijc=0; ijc <jconst.size(); ijc++) jetConstsMap[ jconst[ijc] ] = ev_.ng;
      
      ev_.g_pt[ev_.ng]   = genLep->pt();
      ev_.g_id[ev_.ng]   = genLep->pdgId();
      ev_.g_eta[ev_.ng]  = genLep->eta();
      ev_.g_phi[ev_.ng]  = genLep->phi();
      ev_.g_m[ev_.ng]    = genLep->mass();       
      ev_.ng++;

      //gen level selection
      if(genLep->pt()>20 && fabs(genLep->eta())<2.5) ngleptons++;
    }


  //final state particles
  ev_.ngpf=0;
  if(!runOnGEN_) { // MINIAOD
    edm::Handle<pat::PackedGenParticleCollection> genParticles;
    iEvent.getByToken(genParticlesToken_,genParticles);
    for (size_t i = 0; i < genParticles->size(); ++i)
      {
        const pat::PackedGenParticle & genIt = (*genParticles)[i];

        //this shouldn't be needed according to the workbook
        //if(genIt.status()!=1) continue;
        if(genIt.pt()<0.5 || fabs(genIt.eta())>2.5) continue;

        ev_.gpf_id[ev_.ngpf]     = genIt.pdgId();
        ev_.gpf_c[ev_.ngpf]      = genIt.charge();
        ev_.gpf_g[ev_.ngpf]=-1;
        for(std::map<const reco::Candidate *,int>::iterator it=jetConstsMap.begin();
            it!=jetConstsMap.end();
            it++)
          {
            if(it->first->pdgId()!=genIt.pdgId()) continue;
            if(deltaR( *(it->first), genIt)>0.01) continue;
            ev_.gpf_g[ev_.ngpf]=it->second;
            break;
          }
        ev_.gpf_pt[ev_.ngpf]     = genIt.pt();
        ev_.gpf_eta[ev_.ngpf]    = genIt.eta();
        ev_.gpf_phi[ev_.ngpf]    = genIt.phi();
        ev_.gpf_m[ev_.ngpf]      = genIt.mass();
        ev_.ngpf++;
      }
  }
  else { // GEN
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByToken(prunedGenParticlesToken_,genParticles);
    for (size_t i = 0; i < genParticles->size(); ++i)
      {
        const reco::GenParticle & genIt = (*genParticles)[i];

        if(genIt.status()!=1) continue;
        if(genIt.pt()<0.5 || fabs(genIt.eta())>2.5) continue;

        ev_.gpf_id[ev_.ngpf]     = genIt.pdgId();
        ev_.gpf_c[ev_.ngpf]      = genIt.charge();
        ev_.gpf_g[ev_.ngpf]=-1;
        for(std::map<const reco::Candidate *,int>::iterator it=jetConstsMap.begin();
            it!=jetConstsMap.end();
            it++)
          {
            if(it->first->pdgId()!=genIt.pdgId()) continue;
            if(deltaR( *(it->first), genIt)>0.01) continue;
            ev_.gpf_g[ev_.ngpf]=it->second;
            break;
          }
        ev_.gpf_pt[ev_.ngpf]     = genIt.pt();
        ev_.gpf_eta[ev_.ngpf]    = genIt.eta();
        ev_.gpf_phi[ev_.ngpf]    = genIt.phi();
        ev_.gpf_m[ev_.ngpf]      = genIt.mass();
        ev_.ngpf++;
      }
  }


 //Bhadrons and top quarks (lastCopy)
  edm::Handle<reco::GenParticleCollection> prunedGenParticles;
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  ev_.ngtop=0; 
  for (size_t i = 0; i < prunedGenParticles->size(); ++i)
    {
      const reco::GenParticle & genIt = (*prunedGenParticles)[i];
      int absid=abs(genIt.pdgId());
  
      //top quarks
      if(absid==6 && genIt.isLastCopy()) 
	{
	  ev_.gtop_id[ ev_.ngtop ]  = genIt.pdgId();
	  ev_.gtop_pt[ ev_.ngtop ]  = genIt.pt();
	  ev_.gtop_eta[ ev_.ngtop ] = genIt.eta();
	  ev_.gtop_phi[ ev_.ngtop ] = genIt.phi();
	  ev_.gtop_m[ ev_.ngtop ]   = genIt.mass();
	  ev_.ngtop++;
	}
    }
  
  //gen met
  edm::Handle<reco::METCollection> genMet;
  iEvent.getByToken(genMetsToken_,genMet);
  ev_.gtop_id[ ev_.ngtop ]  = 0;
  ev_.gtop_pt[ ev_.ngtop ]  = (*genMet)[0].pt();
  ev_.gtop_eta[ ev_.ngtop ] = 0;
  ev_.gtop_phi[ ev_.ngtop ] = (*genMet)[0].phi();
  ev_.gtop_m[ ev_.ngtop ]   = 0;
  ev_.ngtop++;

  //fiducial counters
  for(Int_t iw=0; iw<ev_.g_nw; iw++)
    {
      Double_t x(iw);
      Double_t wgt(ev_.g_w[iw]);
      TH2F *fidCounter=(TH2F *)histContainer_["fidcounter"];
      fidCounter->Fill(x,0.,wgt);
      if(ngleptons>0)               fidCounter->Fill(x, 1., wgt);
      if(ngleptons>1)               fidCounter->Fill(x, 2., wgt);
      if(ngleptons>0 && ngjets>0)   fidCounter->Fill(x, 3., wgt);
      if(ngleptons>1 && ngjets>0)   fidCounter->Fill(x, 4., wgt);
      if(ngleptons>0 && ngjets>1)   fidCounter->Fill(x, 5., wgt);
      if(ngleptons>1 && ngjets>1)   fidCounter->Fill(x, 6., wgt);
      if(ngleptons>0 && ngjets>2)   fidCounter->Fill(x, 7., wgt);
      if(ngleptons>1 && ngjets>2)   fidCounter->Fill(x, 8., wgt);
      if(ngleptons>0 && ngjets>3)   fidCounter->Fill(x, 9., wgt);
      if(ngleptons>1 && ngjets>3)   fidCounter->Fill(x, 10.,wgt);
    }
  
  //return number of leptons in fiducial range
  return ngleptons;
}


//
int MiniAnalyzer::recAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  return 0;
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Soft_Muon
bool MiniAnalyzer::isSoftMuon(const reco::Muon & recoMu,const reco::Vertex &vertex)
{

  bool isGood(muon::isGoodMuon(recoMu, muon::TMOneStationTight));
  bool passLayersWithMeas(recoMu.innerTrack().isNonnull()
			  && recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
			  && recoMu.innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0 );
  bool matchesVertex(recoMu.innerTrack().isNonnull()
		     && fabs(recoMu.innerTrack()->dxy(vertex.position())) < 0.3 
		     && fabs(recoMu.innerTrack()->dz(vertex.position())) < 20. );
  return (isGood && passLayersWithMeas && matchesVertex);
}

//cf. https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Standard_MediumID_to_be_used_wit
bool MiniAnalyzer::isMediumMuon2016ReReco(const reco::Muon & recoMu) 
{
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.8 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}



// ------------ method called for each event  ------------
void MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  histContainer_["counter"]->Fill(0);

  //analyze the event
  int ngleptons(0),nrecleptons(0);
  if(!iEvent.isRealData()) ngleptons=genAnalysis(iEvent,iSetup);
  if(!runOnGEN_) nrecleptons=recAnalysis(iEvent,iSetup);
  
  //save event if at least one object at gen or reco level
  if((ngleptons==0 && nrecleptons==0) || !saveTree_) return;  
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event(); 
  ev_.isData  = iEvent.isRealData();
  if(!savePF_) { ev_.ngpf=0; ev_.npf=0; }
  tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob(){
}

//
void 
MiniAnalyzer::endRun(edm::Run const& iRun,
		     const EventSetup& iSetup) 
{

}

//-------------
//cf. https://twiki.cern.ch/twiki/bin/view/CMS/MiniIsolationSUSY
float MiniAnalyzer::getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
				     const reco::Candidate* ptcl,  
				     float r_iso_min, float r_iso_max, float kt_scale,
				     bool charged_only) 
{

    if (ptcl->pt()<5.) return 99999.;

    float deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
    if(ptcl->isElectron()) {
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
    } else if(ptcl->isMuon()) {
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  
    } else {
      //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
    }

    float iso_nh(0.), iso_ch(0.), iso_ph(0.), iso_pu(0.);
    float ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    float r_iso = (float)TMath::Max((float)r_iso_min,
				    (float)TMath::Min((float)r_iso_max, (float)(kt_scale/ptcl->pt())));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;
      
      float dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }
    float iso(0.);
    if (charged_only){
      iso = iso_ch;
    } else {
      iso = iso_ph + iso_nh;
      iso -= 0.5*iso_pu;
      if (iso>0) iso += iso_ch;
      else iso = iso_ch;
    }
    iso = iso/ptcl->pt();

    return iso;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
  std::cout << "[MiniAnalyzer::endJob]" << endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
