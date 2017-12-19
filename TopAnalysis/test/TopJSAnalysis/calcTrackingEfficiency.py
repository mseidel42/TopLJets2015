#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *
from math import sqrt
from math import isnan

def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--task',
                            dest='task',   
                            default='fill',
                            help='task: optimize,fill [default: %default]')
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='eos',
                            help='input file, if directory the script will run in batch mode [default: %default]')
    parser.add_option('-r', '--radius', dest='radius' , help='radius for density', default=0.1, type=float)
    (opt, args) = parser.parse_args()

    tree = ROOT.TChain('analysis/data')
    if opt.input == 'eos':
        opt.input = '/eos/cms/store/cmst3/group/top/ReReco2016/b312177/MC13TeV_TTJets/MergedMiniEvents_0_ext0.root'
    if opt.input == 'local':
        opt.input = '/home/mseidel/Projects/MergedMiniEvents_0_ext0.root'
    if opt.input == 'eosdata':
        opt.input = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177/Chunks/Data13TeV_SingleMuon_2016G_27.root'
    if (tree.Add(opt.input) == 0): return
    
    ROOT.TH1.SetDefaultSumw2(True)
    nbins = 50
    h_dr      = ROOT.TH1F('h_dr',      'h_dr',      nbins, 0, 0.01)
    h_ptr     = ROOT.TH1F('h_ptr',     'h_ptr',     100, 0., 2.)
    h_all     = ROOT.TH1F('h_all',     'h_all',     nbins, 0, nbins)
    h_allSS   = ROOT.TH1F('h_allSS',   'h_allSS',   nbins, 0, nbins)
    h_allOS   = ROOT.TH1F('h_allOS',   'h_allOS',   nbins, 0, nbins)
    h_matched = ROOT.TH1F('h_matched', 'h_matched', nbins, 0, nbins)
    h_matched.GetYaxis().SetRangeUser(0., 2.)
    h_matched.GetXaxis().SetTitle('N_{ch} density (R=%s)'%(str(opt.radius)))
    h_matchedSS = ROOT.TH1F('h_matchedSS', 'h_matchedSS', nbins, 0, nbins)
    h_matchedSS.GetYaxis().SetRangeUser(0., 2.)
    h_matchedSS.GetXaxis().SetTitle('N_{ch} density (R=%s)'%(str(opt.radius)))
    h_matchedOS = ROOT.TH1F('h_matchedOS', 'h_matchedOS', nbins, 0, nbins)
    h_matchedOS.GetYaxis().SetRangeUser(0., 2.)
    h_matchedOS.GetXaxis().SetTitle('N_{ch} density (R=%s)'%(str(opt.radius)))
    
    # plot
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
    # loop
    maxevents = 100000
    ievent = 0
    for ev in tree:
        if not ievent%10:
            print('Processed %i events'%(ievent))
        if not ievent%100:
            h_ratio = h_matched.Clone('h_ratio')
            h_ratio.Divide(h_all)
            h_ratio.Draw('e')
            c.Print('trackingEfficiency_R%s.pdf'%(str(opt.radius)))
            c.Print('trackingEfficiency_R%s.png'%(str(opt.radius)))
            c.Print('trackingEfficiency_R%s.root'%(str(opt.radius)))
            h_ratioSS = h_matchedSS.Clone('h_ratioSS')
            h_ratioSS.Divide(h_allSS)
            h_ratioSS.Draw('e')
            c.Print('trackingEfficiencySS_R%s.pdf'%(str(opt.radius)))
            c.Print('trackingEfficiencySS_R%s.png'%(str(opt.radius)))
            c.Print('trackingEfficiencySS_R%s.root'%(str(opt.radius)))
            h_ratioOS = h_matchedOS.Clone('h_ratioOS')
            h_ratioOS.Divide(h_allOS)
            h_ratioOS.Draw('e')
            c.Print('trackingEfficiencyOS_R%s.pdf'%(str(opt.radius)))
            c.Print('trackingEfficiencyOS_R%s.png'%(str(opt.radius)))
            c.Print('trackingEfficiencyOS_R%s.root'%(str(opt.radius)))
            h_dr.Draw('e')
            c.Print('trackingEfficiency_dR_R%s.pdf'%(str(opt.radius)))
            c.Print('trackingEfficiency_dR_R%s.png'%(str(opt.radius)))
            c.Print('trackingEfficiency_dR_R%s.root'%(str(opt.radius)))
            #h_ptr.Draw('e')
            #c.Print('trackingEfficiency_ptratio_R%s.pdf'%(str(opt.radius)))
            #c.Print('trackingEfficiency_ptratio_R%s.png'%(str(opt.radius)))
            #c.Print('trackingEfficiency_ptratio_R%s.root'%(str(opt.radius)))
            
        ievent += 1
        if ievent > maxevents: break
        
        reco_matched = []
        
        # loop over charged gen hadrons
        for g in range(ev.ngpf):
            if (ev.gpf_pt[g] < 1.): continue
            if (abs(ev.gpf_eta[g]) > 2.4): continue
            if (ev.gpf_c[g] == 0): continue
            if (abs(ev.gpf_id[g]) < 100): continue
            
            density = 0
            densitySS = 0
            densityOS = 0
            
            # calculate charged gen hadron density
            for h in range(ev.ngpf):
                if (ev.gpf_pt[h] < 1.): continue
                if (abs(ev.gpf_eta[h]) > 2.4): continue
                if (ev.gpf_c[h] == 0): continue
                if (abs(ev.gpf_id[h]) < 100): continue
                
                dEta = ev.gpf_eta[g] - ev.gpf_eta[h]
                dPhi = ROOT.TVector2.Phi_mpi_pi(ev.gpf_phi[g] - ev.gpf_phi[h])
                dR = sqrt(pow(dEta, 2) + pow(dPhi, 2))
                if (dR < opt.radius):
                    density += 1
                    if ev.gpf_c[g] == ev.gpf_c[h]: densitySS += 1
                    if ev.gpf_c[g] != ev.gpf_c[h]: densityOS += 1
            
            # match to reco tracks
            matchdr = 10.
            kmatch  = -1.
            for k in range(ev.npf):
                if k in reco_matched: continue
                if (ev.pf_pt[k] < 1.0): continue
                if (abs(ev.pf_eta[k]) > 2.4): continue
                if (abs(ev.pf_id[k]) != 211): continue
                dEta = ev.gpf_eta[g] - ev.pf_eta[k]
                dPhi = ROOT.TVector2.Phi_mpi_pi(ev.gpf_phi[g] - ev.pf_phi[k])
                dR   = sqrt(pow(dEta, 2) + pow(dPhi, 2))
                ptratio = ev.pf_pt[k]/ev.gpf_pt[g]
                if (dR < 0.01 and dR < matchdr and 0.8 < ptratio < 1.2 and ev.pf_c[k] == ev.gpf_c[g]):
                    kmatch = k
                    matchdr = dR
            h_dr.Fill(matchdr)
            #h_ptr.Fill(ptratio)
            reco_matched.append(kmatch)

            # fill
            h_all.Fill(density)
            if kmatch > -1:
                h_matched.Fill(density)
            h_allSS.Fill(densitySS)
            if kmatch > -1:
                h_matchedSS.Fill(densitySS)
            h_allOS.Fill(densityOS)
            if kmatch > -1:
                h_matchedOS.Fill(densityOS)
    
    
    

if __name__ == "__main__":
	sys.exit(main())
