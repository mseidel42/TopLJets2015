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
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=35922.,              type=float)
    parser.add_option('--var', dest='var', default='j_pt[0]', help='observable [default: %default]')
    parser.add_option('--sel', dest='sel', default='1', help='selection [default: %default]')
    (opt, args) = parser.parse_args()
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    eosdir = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177_new/Chunks/'
    
    t_data = ROOT.TChain('tjsev')
    t_data.Add(eosdir+'Data13TeV_*.root')
    print('t_data.GetEntries()', t_data.GetEntries())

    t_mc = ROOT.TChain('tjsev')
    for i in range(137):
        t_mc.Add(eosdir+'MC13TeV_TTJets_'+str(i)+'.root')
    print('t_mc.GetEntries()', t_mc.GetEntries())

    t_data.Draw(opt.var+" >> h_data(50,0,250)", 'reco_sel == 1 & ' + opt.sel)
    h_data = ROOT.gDirectory.Get("h_data")
    h_data.SetLineColor(ROOT.kBlack)

    t_mc.Draw(opt.var+" >> h_mc(50,0,250)", '(reco_sel == 1 & '+opt.sel+')*weight[0]')
    h_mc = ROOT.gDirectory.Get("h_mc")
    h_mc.Scale(opt.lumi*832.)
    h_mc.SetLineColor(ROOT.kRed+1)
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('+.2f') 
    c = ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()
    
    h_data.Draw()
    h_mc.Draw('same')
    
    filename = filter(str.isalnum, opt.var) + '_' + filter(str.isalnum, opt.sel)
    c.Print('quickplots/%s.pdf'%(filename))
    c.Print('quickplots/%s.png'%(filename))


if __name__ == "__main__":
	sys.exit(main())
