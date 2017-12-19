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

def normalize(hist):
    hist.Scale(1./hist.Integral())

def main():
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=35922.,              type=float)
    parser.add_option('--var', dest='var', default='j_mult_charged', help='observable [default: %default]')
    parser.add_option('--sel', dest='sel', default='1', help='selection [default: %default]')
    parser.add_option('--bin', dest='bin', default='20,0,20', help='selection [default: %default]')
    parser.add_option('-n', '--nevents', dest='nevents', help='number of events for draw',              default=ROOT.TVirtualTreePlayer.kMaxEntries, type=long)
    (opt, args) = parser.parse_args()
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    eosdir = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177_new/Chunks/'
    
    t_data_bf = ROOT.TChain('tjsev')
    for trigger in ['SingleElectron', 'SingleMuon']:
        for run in ['B', 'C', 'D', 'E', 'F']:
            if t_data_bf.GetEntries() > opt.nevents: break
            t_data_bf.Add(eosdir+'Data13TeV_'+trigger+'_2016'+run+'_*.root')
    
    t_data_gh = ROOT.TChain('tjsev')
    for trigger in ['SingleElectron', 'SingleMuon']:
        for run in ['G', 'H']:
            if t_data_gh.GetEntries() > opt.nevents: break
            t_data_gh.Add(eosdir+'Data13TeV_'+trigger+'_2016'+run+'_*.root')

    t_mc = ROOT.TChain('tjsev')
    for i in range(137):
        if t_mc.GetEntries() > opt.nevents: break
        t_mc.Add(eosdir+'MC13TeV_TTJets_'+str(i)+'.root')
    print('t_mc.GetEntries()', t_mc.GetEntries())

    t_data_bf.Draw(opt.var+' >> h_data_bf('+opt.bin+')', 'reco_sel == 1 & ' + opt.sel, "", opt.nevents)
    h_data_bf = ROOT.gDirectory.Get("h_data_bf")
    normalize(h_data_bf)
    h_data_bf.SetLineColor(ROOT.kBlack)
    h_data_bf.SetMarkerColor(ROOT.kBlack)
    h_data_bf.SetMarkerStyle(20)
    h_data_bf.SetTitle('')
    h_data_bf.GetXaxis().SetTitle('N_{ch}')
    h_data_bf.GetYaxis().SetRangeUser(0.01, h_data_bf.GetMaximum()*1.6)
    
    t_data_gh.Draw(opt.var+' >> h_data_gh('+opt.bin+')', 'reco_sel == 1 & ' + opt.sel, "", opt.nevents)
    h_data_gh = ROOT.gDirectory.Get("h_data_gh")
    normalize(h_data_gh)
    h_data_gh.SetLineColor(ROOT.kGray)
    h_data_gh.SetMarkerColor(ROOT.kGray)
    h_data_gh.SetMarkerStyle(20)

    t_mc.Draw(opt.var+' >> h_mc_bf('+opt.bin+')', '(reco_sel == 1 & period == 1 & '+opt.sel+')*weight[0]', "", opt.nevents)
    h_mc_bf = ROOT.gDirectory.Get("h_mc_bf")
    normalize(h_mc_bf)
    h_mc_bf.SetLineColor(ROOT.kRed+1)
    
    t_mc.Draw(opt.var+' >> h_mc_gh('+opt.bin+')', '(reco_sel == 1 & period == 2 & '+opt.sel+')*weight[0]', "", opt.nevents)
    h_mc_gh = ROOT.gDirectory.Get("h_mc_gh")
    normalize(h_mc_gh)
    h_mc_gh.SetLineColor(ROOT.kBlue+1)
    h_mc_gh.SetLineStyle(7)
    
    # Plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0)
    c.SetRightMargin(0.00)
    c.cd()
    
    p1=ROOT.TPad('p1','p1',0.0,0.32,1.0,1.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.06)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    
    h_data_bf.GetYaxis().SetLabelSize(0.05)
    h_data_bf.Draw()
    h_data_gh.Draw('same')
    h_mc_bf.Draw('same,hist')
    h_mc_gh.Draw('same,hist')
    
    legend = ROOT.TLegend(0.5,0.6,0.85,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h_data_bf, 'Data B-F, mean = %.2f'%(h_data_bf.GetMean()), "pl")
    legend.AddEntry(h_data_gh, 'Data GH, mean = %.2f'%(h_data_gh.GetMean()), "pl")
    legend.AddEntry(h_mc_bf, 'MC B-F, mean = %.2f'%(h_mc_bf.GetMean()), "l")
    legend.AddEntry(h_mc_gh, 'MC GH, mean = %.2f'%(h_mc_gh.GetMean()), "l")
    legend.Draw()
    
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.2,1.0,0.32)
    p2.Draw()
    p2.SetBottomMargin(0.02)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.00)
    p2.cd()
    
    ratio_bf = h_mc_bf.Clone()
    ratio_bf.Divide(h_data_bf)
    ratio_bf.SetTitle('')
    ratio_bf.SetXTitle(h_data_bf.GetXaxis().GetTitle())
    ratio_bf.SetYTitle('MC/data')
    ratio_bf.SetFillColor(ROOT.kGray)
    ratio_bf.GetXaxis().SetTitleSize(0.2)
    ratio_bf.GetXaxis().SetTitleOffset(0.8)
    ratio_bf.GetXaxis().SetLabelSize(0.18)
    ratio_bf.GetYaxis().SetTitleSize(0.16/0.6)
    ratio_bf.GetYaxis().SetTitleOffset(0.3*0.6)
    ratio_bf.GetYaxis().SetLabelSize(0.16/0.6)
    ratio_bf.GetYaxis().SetRangeUser(0.4,1.6)
    ratio_bf.GetYaxis().SetNdivisions(503)
    ratio_bf.Draw()
    
    ratio_gh = h_mc_gh.Clone()
    ratio_gh.Divide(h_data_gh)
    ratio_gh.Draw('same')
    
    c.cd()
    p3 = ROOT.TPad('p3','p3',0.0,0.0,1.0,0.2)
    p3.Draw()
    p3.SetBottomMargin(0.4)
    p3.SetRightMargin(0.05)
    p3.SetLeftMargin(0.12)
    p3.SetTopMargin(0.01)
    p3.cd()
    
    ratio_periods = ratio_bf.Clone()
    ratio_periods.Divide(ratio_gh)
    ratio_periods.SetYTitle('B-F/GH ')
    ratio_periods.GetYaxis().SetRangeUser(0.88,1.12)
    ratio_periods.GetYaxis().SetTitleSize(0.16)
    ratio_periods.GetYaxis().SetTitleOffset(0.3)
    ratio_periods.GetYaxis().SetLabelSize(0.16)
    ratio_periods.Draw()
    
    filename = filter(str.isalnum, opt.var) + '_' + filter(str.isalnum, opt.sel)
    c.Print('quickplots/%s.pdf'%(filename))
    c.Print('quickplots/%s.png'%(filename))


if __name__ == "__main__":
	sys.exit(main())
