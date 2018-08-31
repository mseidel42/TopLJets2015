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
    parser.add_option('-n', '--nevents', dest='nevents', help='number of events', default=10000, type=int)
    parser.add_option('--flavor', dest='flavor', default='incl', help='jet flavor')
    (opt, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2(True)

    #read lists of samples
    flavormap = {'bottom': 5, 'light': 1, 'gluon': 0}
    flavor_label = {'incl': 'inclusive', 'bottom': 'bottom', 'light': 'light-enriched', 'gluon': 'gluon-enriched'}

    generators = ['pythia', 'fsrup', 'fsrdn', 'herwig7', 'sherpa', 'dire2002']
    name = {'pythia': 'POWHEG+PYTHIA 8', 'fsrup': 'PP8 FSR up', 'fsrdn': 'PP8 FSR down', 'pythia8_asfsr0.120_meon_crdefault': '#minus FSR tuned', 'herwig7': 'POWHEG+HERWIG 7', 'sherpa': 'SHERPA 2', 'dire2002': 'DIRE 2'}
    color = {'pythia': ROOT.kRed+1, 'fsrup': ROOT.kRed+1, 'fsrdn': ROOT.kRed+1, 'pythia8_asfsr0.120_meon_crdefault': ROOT.kOrange+1, 'herwig7': ROOT.kAzure+2, 'sherpa': ROOT.kGreen+2, 'dire2002': ROOT.kMagenta+3}
    marker = {'pythia': 24, 'fsrup': 26, 'fsrdn': 32, 'pythia8_asfsr0.120_meon_crdefault': 0, 'herwig7': 25, 'sherpa': 27, 'dire2002': 30}
    line = {'pythia': 0, 'fsrup': 0, 'fsrdn': 0, 'pythia8_asfsr0.120_meon_crdefault': 2, 'herwig7': 7, 'sherpa': 5, 'dire2002': 2}
    draw = {'pythia': 'h', 'fsrup': 'same p x0 e1', 'fsrdn': 'same p x0 e1', 'pythia8_asfsr0.120_meon_crdefault': 'same h', 'herwig7': 'same h', 'sherpa': 'same h', 'dire2002': 'same h'}
    legend_opt = {'pythia': 'pl', 'fsrup': 'p', 'fsrdn': 'p', 'pythia8_asfsr0.120_meon_crdefault': 'pl', 'herwig7': 'pl', 'sherpa': 'pl', 'dire2002': 'pl'}

    hist = {}
    ratio = {}

    eosdir = '/eos/cms/store/group/cmst3/user/mseidel/analysis/TopJetShapes/b312177/Chunks/'

    p1correction = 0.8/0.69
    p2correction = 0.2/0.31

    for generator in generators:
        maxweight = 100000
        if generator == 'dire2002': maxweight = 1000
        t_mc = ROOT.TChain('tjsev')
        if generator == 'pythia':
            for i in range(137):
                t_mc.Add(eosdir+'MC13TeV_TTJets_'+str(i)+'.root')
        else:
            t_mc.Add(eosdir+'MC13TeV_TTJets_'+generator+'_*.root')

        hist[generator] = ROOT.TH1F(generator, ';Jet p_{T} [GeV];1/N_{jet} dN_{jet} / d p_{T}', 10,30,230)
        if opt.flavor == 'incl':
            t_mc.Draw(opt.var+' >> '+generator, '(gen_sel == 1 & abs(gj_eta) < 2. & gj_overlap == 0 & abs(weight[0])<%i)*weight[0]'%(maxweight), '', opt.nevents)
        else:
            t_mc.Draw(opt.var+' >> '+generator, '(gen_sel == 1 & abs(gj_eta) < 2. & gj_overlap == 0 & gj_flavor == '+str(flavormap[opt.flavor])+' & abs(weight[0])<%i)*weight[0]'%(maxweight), '', opt.nevents)

        hist[generator].Scale(1./hist[generator].Integral())

        hist[generator].GetXaxis().SetTitleSize(0)
        hist[generator].GetXaxis().SetLabelSize(0)
        hist[generator].GetYaxis().SetTitleSize(0.05*p1correction)
        hist[generator].GetYaxis().SetLabelSize(0.045*p1correction)
        hist[generator].GetYaxis().SetTitleOffset(1.2/p1correction)
        hist[generator].SetLineColor(color[generator])
        hist[generator].SetLineStyle(line[generator])
        hist[generator].SetLineWidth(2)
        hist[generator].SetMarkerColor(color[generator])
        hist[generator].SetMarkerStyle(marker[generator])

        ratio[generator] = hist[generator].Clone(generator+'_ratio')
        ratio[generator].Divide(hist['pythia'])

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('+.2f')
    c = ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0)
    c.SetRightMargin(0.00)
    c.cd()

    p1=ROOT.TPad('p1','p1',0.0,0.31,1.0,1.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.07)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    #ROOT.gPad.SetLogy()

    hist['pythia'].GetYaxis().SetRangeUser(0.00001, hist['pythia'].GetMaximum()*1.75)
    hist['pythia'].Draw()
    for generator in generators:
        hist[generator].Draw(draw[generator])

    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.31)
    p2.Draw()
    p2.SetBottomMargin(0.29)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.01)
    p2.cd()

    ratio['pythia'].SetTitle('')
    ratio['pythia'].SetYTitle('MC/PP8   ')
    ratio['pythia'].GetXaxis().SetTitleSize(0.2*p2correction)
    ratio['pythia'].GetXaxis().SetTitleOffset(0.95)
    ratio['pythia'].GetXaxis().SetLabelSize(0.18*p2correction)
    ratio['pythia'].GetYaxis().SetTitleSize(0.2*p2correction)
    ratio['pythia'].GetYaxis().SetTitleOffset(0.26/p2correction)
    ratio['pythia'].GetYaxis().SetLabelSize(0.18*p2correction)
    ratio['pythia'].GetYaxis().SetNdivisions(503)
    ratio['pythia'].GetYaxis().SetRangeUser(0.35,1.65)
    ratio['pythia'].Draw(draw['pythia'])
    for generator in generators:
        ratio[generator].Draw(draw[generator])

    c.cd()

    legend = ROOT.TLegend(0.55,0.6,0.95,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    for generator in generators:
        legend.AddEntry(hist[generator], name[generator], legend_opt[generator])
    legend.Draw()

    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.15,0.9,'#scale[1.2]{%s}'%cmsLabel)
    txt.DrawLatex(0.15,0.85,'#scale[1.0]{t#bar{t} #rightarrow lepton+jets}')
    txt.DrawLatex(0.15,0.8, '#scale[1.0]{'+flavor_label[opt.flavor]+' jets}')

    filename = filter(str.isalnum, opt.var)
    c.Print('quickplots/%s_%s.pdf'%(filename, opt.flavor))
    c.Print('quickplots/%s_%s.png'%(filename, opt.flavor))


if __name__ == "__main__":
	sys.exit(main())
