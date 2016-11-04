#ifndef _topjetshape_h_
#define _topjetshape_h_

#include "TLorentzVector.h"

struct TopJetShapeEvent_t
{
  Int_t nw, nl, ngl, nj, ngj;
  Int_t gen_sel, reco_sel;
  
  Float_t weight[10];
  
  Float_t l_pt[5], l_eta[5], l_phi[5], l_m[5];
  Int_t l_id[5];
  
  Float_t gl_pt[5], gl_eta[5], gl_phi[5], gl_m[5];
  Int_t gl_id[5];
  
  Float_t j_pt[50], j_eta[50], j_phi[50], j_m[50];
  Int_t j_flavor[50], j_overlap[50], j_gj[50];
  Float_t j_ga[50][3][3][3][3];
  
  Float_t gj_pt[50], gj_eta[50], gj_phi[50], gj_m[50];
  Int_t gj_flavor[50], gj_overlap[50], gj_j[50];
  Float_t gj_ga[50][3][3][3][3];
  
  Float_t met_pt,met_phi;
};

struct Particle {
  TLorentzVector p4;
  int charge;
  double puppi;
  
  Particle(TLorentzVector p, int c, double pu)
  : p4(p), charge(c), puppi(pu) {}
  
  double pt()  { return p4.Pt();  }
  double eta() { return p4.Eta(); }
  double phi() { return p4.Phi(); }
  double energy() { return p4.E(); }
  double mass()   { return p4.M(); }
  TLorentzVector momentum() { return p4; }
};

struct Jet {
  TLorentzVector p4;
  std::vector<Particle> particles;
  int flavor;
  int oldidx;
  int overlap;
  
  double ga[3][3][3][3]; // beta/kappa/iptcut/icharge
  
  Jet(TLorentzVector p, int f, int oi)
  : p4(p), flavor(f), oldidx(oi), overlap(0) {}
  
  TLorentzVector momentum() { return p4; }
};

double calcGA(Jet jet, int beta, int kappa, int iptcut, int icharge);
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev);
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev);

double deltaR(TLorentzVector v1, TLorentzVector v2) { return v1.DeltaR(v2); }
double mapAngleMPiToPi(double phi) { return TVector2::Phi_mpi_pi(phi); }
#endif
