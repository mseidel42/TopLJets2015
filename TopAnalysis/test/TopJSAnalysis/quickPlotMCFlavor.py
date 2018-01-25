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
    
    cmsLabel='#bf{CMS} #it{Simulation}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-l', '--lumi', dest='lumi', help='lumi [/pb]', default=35922., type=float)
    parser.add_option('--var', dest='var', default='gj_pt', help='observable [default: %default]')
    parser.add_option('-n', '--nevents', dest='nevents', help='number of events', default=-1, type=int)
    parser.add_option('--sample', dest='sample', default='', help='sample appendix')
    (opt, args) = parser.parse_args()
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    #read lists of samples
    flavors = ['incl', 'bottom', 'light', 'gluon']
    flavormap = {'bottom': 5, 'light': 1, 'gluon': 0}
    # ColorBrewer 3-class Set2 (http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3)
    colors = {'incl': ROOT.kBlack, 'bottom': ROOT.kOrange+8, 'light': ROOT.kBlue-5, 'gluon': ROOT.kTeal-8}
    lightcolors = {'incl': ROOT.kGray, 'bottom': ROOT.kOrange-9, 'light': ROOT.kBlue-10, 'gluon': ROOT.kCyan-10}
    markers = {'incl': 20, 'bottom': 21, 'light': 22, 'gluon': 23}
    fills = {'incl': 1001, 'bottom': 3254, 'light': 3245, 'gluon': 3002}
    hist = {}
    ratio = {}
    
    eosdir = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177_new/Chunks/'

    t_mc = ROOT.TChain('tjsev')
    if opt.sample == '':
        for i in range(137):
            t_mc.Add(eosdir+'MC13TeV_TTJets_'+str(i)+'.root')
    else:
        t_mc.Add(eosdir+'MC13TeV_TTJets_'+opt.sample+'_*.root')
    print('t_mc.GetEntries()', t_mc.GetEntries())
    
    for flavor in flavors:
        hist[flavor] = ROOT.TH1F(flavor, ';Jet p_{T} [GeV];1/N_{jet} dN_{jet} / d p_{T}', 30,0,300)
        if flavor == 'incl':
            t_mc.Draw(opt.var+' >> '+flavor, '(gen_sel == 1 & abs(gj_eta) < 2. & gj_overlap == 0)*weight[0]', '', opt.nevents)
        else:
            t_mc.Draw(opt.var+' >> '+flavor, '(gen_sel == 1 & abs(gj_eta) < 2. & gj_overlap == 0 & gj_flavor == '+str(flavormap[flavor])+')*weight[0]', '', opt.nevents)
        
        hist[flavor].Scale(1./hist[flavor].Integral())
        
        hist[flavor].GetXaxis().SetTitleSize(0.045)
        hist[flavor].GetXaxis().SetLabelSize(0.04)
        hist[flavor].GetYaxis().SetTitleSize(0.05)
        hist[flavor].GetYaxis().SetLabelSize(0.045)
        hist[flavor].GetYaxis().SetTitleOffset(1.2)
        hist[flavor].SetLineColor(colors[flavor])
        hist[flavor].SetMarkerColor(colors[flavor])
        hist[flavor].SetMarkerStyle(markers[flavor])
        print(flavor, hist[flavor].GetEntries())
        
        ratio[flavor] = hist[flavor].Clone(flavor+'_ratio')
        ratio[flavor].Divide(hist['incl'])
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('+.2f') 
    c = ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0)
    c.SetRightMargin(0.00)
    c.cd()
    
    p1=ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.06)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    ROOT.gPad.SetLogy()
    
    hist['incl'].GetYaxis().SetRangeUser(0.0002, hist['incl'].GetMaximum()*25)
    hist['incl'].Draw()
    for flavor in flavors:
        hist[flavor].Draw('same')
    
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
    p2.Draw()
    p2.SetBottomMargin(0.4)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.01)
    p2.cd()
    
    ratio['incl'].SetTitle('')
    ratio['incl'].SetYTitle('#frac{flavor}{incl} ')
    ratio['incl'].SetFillColor(ROOT.kGray)
    ratio['incl'].GetXaxis().SetTitleSize(0.2)
    ratio['incl'].GetXaxis().SetTitleOffset(0.8)
    ratio['incl'].GetXaxis().SetLabelSize(0.18)
    ratio['incl'].GetYaxis().SetTitleSize(0.2)
    ratio['incl'].GetYaxis().SetTitleOffset(0.25)
    ratio['incl'].GetYaxis().SetLabelSize(0.18)
    ratio['incl'].GetYaxis().SetRangeUser(0.4,1.6)
    ratio['incl'].GetYaxis().SetNdivisions(503)
    ratio['incl'].GetYaxis().SetRangeUser(0.,2.25)
    ratio['incl'].Draw()
    for flavor in flavors:
        ratio[flavor].Draw('same')
    
    c.cd()
    
    legend = ROOT.TLegend(0.6,0.7,0.95,0.925)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hist['incl'], 'Inclusive jets',  'lep')
    legend.AddEntry(hist['bottom'], 'Bottom jets',   'lep')
    legend.AddEntry(hist['light'], 'Light-enriched', 'lep')
    legend.AddEntry(hist['gluon'], 'Gluon-enriched', 'lep')
    legend.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.91,cmsLabel)
    txt.DrawLatex(0.15,0.875,'#scale[0.9]{t#bar{t} #rightarrow lepton+jets}')
    
    filename = filter(str.isalnum, opt.var) + opt.sample
    c.Print('quickplots/%s_flavors.pdf'%(filename))
    c.Print('quickplots/%s_flavors.png'%(filename))


if __name__ == "__main__":
	sys.exit(main())
