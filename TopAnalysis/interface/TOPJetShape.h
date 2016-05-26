#ifndef _topjetshape_h_
#define _topjetshape_h_

#include "TLorentzVector.h"

struct TopJetShapeEvent_t
{
  Int_t cat,nw,nl,nj,nt,ngj,ngl;
  Float_t weight[10];
  Float_t l_pt[5],l_eta[5],l_phi[5],l_m[5];
  Int_t l_id[5];
  Float_t gl_pt[5],gl_eta[5],gl_phi[5],gl_m[5];
  Int_t gl_id[5];
  Float_t j_pt[50],j_eta[50],j_phi[50],j_m[50];
  Float_t t_pt[4],t_eta[4],t_phi[4],t_m[4];
  Int_t t_id[4];
  Float_t met_pt,met_phi;
  
  Int_t gen_sel;
  Float_t gj_pt[50],gj_eta[50],gj_phi[50],gj_m[50];
  Int_t gj_flav[50],gj_hadflav[50];
  Float_t gj_ga[50][3][3][3][3];
  
};

struct Particle {
  TLorentzVector p4;
  int charge;
  double puppi;
  
  Particle(TLorentzVector p, int c, double pu)
  : p4(p), charge(c), puppi(pu) {}
};

struct Jet {
  TLorentzVector p4;
  std::vector<Particle> particles;
  int flav;
  int oldidx;
  
  double ga[3][3][3][3]; // beta/kappa/iptcut/icharge
  
  Jet(TLorentzVector p, int f, int oi)
  : p4(p), flav(f), oldidx(oi) {}
};

enum {charged, all, puppi};
enum {pt0, pt500, pt1000};

double calcGA(Jet jet, int beta, int kappa, int iptcut, int icharge);
void createTopJetShapeEventTree(TTree *t,TopJetShapeEvent_t &tjsev);
void resetTopJetShapeEvent(TopJetShapeEvent_t &tjsev);
#endif
