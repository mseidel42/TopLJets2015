#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
import itertools
import copy
from array import *

"""
Get correlation factor
"""
def getCorrelation(points, obs1, obs2):
    h = ROOT.TH2D('', '', 100, 0, 1, 100, 0, 1)
    for i in range(len(points[obs1])):
        val1 = points[obs1][i]
        val2 = points[obs2][i]
        
        if any(x in obs1 for x in ['mult', 'nsd']): val1 /= 100.
        if any(x in obs2 for x in ['mult', 'nsd']): val2 /= 100.
        if any(x in obs1 for x in ['n3_b1', 'n3_b2']): val1 /= 5.
        if any(x in obs2 for x in ['n3_b1', 'n3_b2']): val2 /= 5.
        
        if (val1 < 0.) or (val2 < 0.): continue
        
        h.Fill(val1, val2)
    
    return h.GetCorrelationFactor()

"""
Get data points from tree
"""
def getPointsFromTree(tree, observables, recos):
    points = {}
    for obs in observables:
        for reco in recos:
            points[obs+'_'+reco] = []
    count = 0
    for event in tree:
        count += 1
        if count > 10000: break
        if event.gen_sel != 1: continue
        for j in range(event.ngj):
            if event.gj_overlap == 1: continue
            for obs in observables:
                for reco in recos:
                    val = eval('event.gj_'+obs+'_'+reco)[j]
                    points[obs+'_'+reco].append(val)
    
    return points
    
"""
steer
"""
def main():

    cmsLabel='#bf{CMS} #it{Simulation} #it{Prel.}'
    
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file [default: %default]')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='',
                            help='Output directory [default: %default]')
    parser.add_option('--ro', '--rootoutput',
                            dest='rootoutput',
                            default='correlations.root',
                            help='output root file [default: %default]')
    parser.add_option('-r', '--reco',
                            dest='reco',
                            default='charged',
                            help='Use charged/puppi/all particles [default: %default]')
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=35922.,              type=float)
    parser.add_option(     '--keep',          dest='keep'  ,        help='observables to keep [default: %default]', default='ga_width_charged',       type='string')
    parser.add_option('-N', '--nobs',        dest='nobs' ,       help='Target number of observables',              default=4,              type=int)
    (opt, args) = parser.parse_args()
    
    tree = ROOT.TChain('tjsev')
    tree.Add(opt.input)

    rootoutfile = ROOT.TFile(opt.rootoutput, "RECREATE");
    
    recos = []
    for reco in opt.reco.split(','):
        recos.append(reco)
    
    #observables = ["mult", "width", "ptd", "ptds", "ecc", "zg", "zgxdr", "zgdr", "nsd", "ga_width", "ga_lha", "ga_thrust", "tau21", "tau32", "tau43", "c1_00", "c1_02", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "m2_b2", "n2_b1", "n2_b2", "n3_b1", "n3_b2"]
    observables = ["mult", "ptds", "ga_width", "ga_lha", "ga_thrust", "ecc", "zg", "zgdr", "nsd", "tau21", "tau32", "tau43", "c1_00", "c1_02", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "m2_b2", "n2_b1", "n2_b2", "n3_b1", "n3_b2"]
    
    obsreco = []
    for reco in recos:
        for obs in observables:
            obsreco.append(obs+'_'+reco)
    
    nice_observables_root = {"mult": "#lambda_{0}^{0} (N)", "width": "width", "ptd": "p_{T}D", "ptds": "#lambda_{0}^{2}* (p_{T}D*)", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "#DeltaR_{g}", "ga_width": "#lambda_{1}^{1} (width)", "ga_lha": "#lambda_{0.5}^{1} (LHA)", "ga_thrust": "#lambda_{2}^{1} (thrust)", "c1_00": "C_{1}^{(0.0)}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1.0)}", "c1_20": "C_{1}^{(2.0)}", "c2_00": "C_{2}^{(0.0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1.0)}", "c2_20":  "C_{2}^{(2.0)}", "c3_00": "C_{3}^{(0.0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1.0)}", "c3_20": "C_{3}^{(2.0)}", "m2_b1": "M_{ 2}^{ (1)}", "n2_b1": "N_{ 2}^{ (1)}", "n3_b1": "N_{ 3}^{ (1)}", "m2_b2": "M_{ 2}^{ (2)}", "n2_b2": "N_{ 2}^{ (2)}", "n3_b2": "N_{ 3}^{ (2)}", "nsd": "n_{SD}"}
    
    nice_observables_root_short = {"mult": "#lambda_{0}^{0}", "width": "#lambda_{1}^{1}", "ptd": "#lambda_{0}^{2}", "ptds": "#lambda_{0}^{2}*", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "#DeltaR_{g}", "ga_width": "#lambda_{1}^{1}", "ga_lha": "#lambda_{0.5}^{1}", "ga_thrust": "#lambda_{2}^{1}", "c1_00": "C_{1}^{(0.0)}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1.0)}", "c1_20": "C_{1}^{(2.0)}", "c2_00": "C_{2}^{(0.0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1.0)}", "c2_20":  "C_{2}^{(2.0)}", "c3_00": "C_{3}^{(0.0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1.0)}", "c3_20": "C_{3}^{(2.0)}", "m2_b1": "M_{ 2}^{ (1)}", "n2_b1": "N_{ 2}^{ (1)}", "n3_b1": "N_{ 3}^{ (1)}", "m2_b2": "M_{ 2}^{ (2)}", "n2_b2": "N_{ 2}^{ (2)}", "n3_b2": "N_{ 3}^{ (2)}", "nsd": "n_{SD}"}
    
    points = getPointsFromTree(tree, observables, recos)
    
    # 11-class RdBu http://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11
    stops = array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    red   = array('d')
    green = array('d')
    blue  = array('d')
    colors = [[103,0,31],
              [178,24,43],
              [214,96,77],
              [244,165,130],
              [253,219,199],
              [247,247,247],
              [209,229,240],
              [146,197,222],
              [67,147,195],
              [33,102,172],
              [5,48,97]]
    for color in colors:
        red.append(color[0]/255.)
        green.append(color[1]/255.)
        blue.append(color[2]/255.)
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 30)
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('+2.0f') 
    c = ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()
    
    nbins = len(obsreco)
    h_correlations = ROOT.TH2D('h_correlations', '', nbins, 0, nbins, nbins, 0, nbins)
    
    for i in range(len(obsreco)):
        h_correlations.GetXaxis().SetBinLabel(i+1, nice_observables_root[obsreco[i].rsplit('_', 1)[0]])
        h_correlations.GetYaxis().SetBinLabel(i+1, nice_observables_root[obsreco[i].rsplit('_', 1)[0]])
        for j in range(i+1):
            correlation = 100.*getCorrelation(points, obsreco[i], obsreco[j])
            h_correlations.SetBinContent(i+1, j+1, correlation)
            h_correlations.SetBinContent(j+1, i+1, correlation)
    
    h_correlations.GetXaxis().LabelsOption('v')
    h_correlations.GetXaxis().SetLabelSize(0.02)
    h_correlations.GetYaxis().SetLabelSize(0.02)
    h_correlations.GetZaxis().SetRangeUser(-100., 100.)
    h_correlations.GetZaxis().SetTitle('Correlation [%]     ')
    h_correlations.GetZaxis().SetTitleSize(0.04)
    h_correlations.SetMarkerSize(0.5)
    h_correlations.Draw('colz,text')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.14,0.97, '#scale[1.1]{%s}'%cmsLabel)
    txt.DrawLatex(0.59,0.97, '#scale[1.0]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.Print('correlations.pdf')
    c.Print('correlations.png')
    
    # actually try to remove some that cannot be measured too well by subjective criteria...
    blacklist = ['n3_b1_charged', 'n3_b1_all']
    if opt.keep == 'none':
        keeplist = []
    else:
        keeplist  = opt.keep.split(',')
    
    # algorithm deleting observables with highest correlation
    observables_low = copy.copy(obsreco)
    for obs in blacklist:
        try:
            observables_low.remove(obs)
        except:
            print('not in observable list', obs)
    while (len(observables_low) > 15):
        maxCorrelation = 0.
        maxCorrelationPair = None
        for pair in itertools.combinations(range(len(obsreco)),2):
            if not obsreco[pair[0]] in observables_low: continue
            if not obsreco[pair[1]] in observables_low: continue
            pairCorrelation = abs(h_correlations.GetBinContent(pair[0]+1, pair[1]+1))
            if (pairCorrelation > maxCorrelation):
                maxCorrelation = pairCorrelation
                maxCorrelationPair = pair
        sumcorr = {}
        for i in maxCorrelationPair:
            sumcorr[i] = 0.
            for j in range(len(obsreco)):
                if not obsreco[j] in observables_low: continue
                sumcorr[i] += abs(h_correlations.GetBinContent(i+1, j+1))
        maxCorrelationObs = max(sumcorr.iterkeys(), key=(lambda key: sumcorr[key]))
        minCorrelationObs = min(sumcorr.iterkeys(), key=(lambda key: sumcorr[key]))
        print('remove', maxCorrelationObs, obsreco[maxCorrelationObs], maxCorrelation, minCorrelationObs, obsreco[minCorrelationObs])
        if obsreco[maxCorrelationObs] in keeplist:
            observables_low.remove(obsreco[minCorrelationObs])
        else:
            observables_low.remove(obsreco[maxCorrelationObs])
    print(observables_low)
    
    # finding set of observables with smallest global correlation now
    bestCombination = None
    bestSumCorr = 9999999.
    for combination in itertools.combinations(observables_low, opt.nobs):
        if not set(keeplist).issubset(combination): continue
        sumcorr = 0.
        for si in combination:
            for sj in combination:
                corr = abs(h_correlations.GetBinContent(obsreco.index(si)+1, obsreco.index(sj)+1))
                if corr > 30.: sumcorr += 1000. # punish correlations >30%
                else:           sumcorr += corr
        if (sumcorr < bestSumCorr):
            bestSumCorr = sumcorr
            bestCombination = combination
    observables_low = list(bestCombination)
    #observables_low = ["ptds", "tau43", "zg", "zgdr", "n3_b2", "n2_b2"]
    print(observables_low)
    
    h_correlations_low = ROOT.TH2D('h_correlations_low', '', len(observables_low), 0, len(observables_low), len(observables_low), 0, len(observables_low))
    
    for i in range(len(observables_low)):
        h_correlations_low.GetXaxis().SetBinLabel(i+1, nice_observables_root_short[observables_low[i].rsplit('_', 1)[0]])
        h_correlations_low.GetYaxis().SetBinLabel(i+1, nice_observables_root_short[observables_low[i].rsplit('_', 1)[0]])
        for j in range(len(observables_low)):
            correlation = 100.*getCorrelation(points, observables_low[i], observables_low[j])
            h_correlations_low.SetBinContent(i+1, j+1, correlation)
    
    h_correlations_low.GetXaxis().LabelsOption('h')
    h_correlations_low.GetXaxis().SetLabelSize(0.07)
    h_correlations_low.GetYaxis().SetLabelSize(0.07)
    h_correlations_low.GetZaxis().SetRangeUser(-100., 100.)
    h_correlations_low.GetZaxis().SetTitle('Correlation [%]     ')
    h_correlations_low.GetZaxis().SetTitleSize(0.04)
    h_correlations_low.SetMarkerSize(2.)
    h_correlations_low.Draw('colz,text')
    
    ROOT.gPad.Update()
    tl1 = ROOT.TLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl1.Draw()
    tl2 = ROOT.TLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                     ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax())
    tl2.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.14,0.97, '#scale[1.1]{%s}'%cmsLabel)
    txt.DrawLatex(0.59,0.97, '#scale[1.0]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    
    c.Print('correlations_low_%s_%s.pdf'%(filter(str.isalnum, opt.keep), opt.nobs))
    c.Print('correlations_low_%s_%s.png'%(filter(str.isalnum, opt.keep), opt.nobs))

if __name__ == "__main__":
	sys.exit(main())
