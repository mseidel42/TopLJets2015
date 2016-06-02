#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *

"""
Analysis loop
"""
def optimize(fileName):
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.05);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetTitleXOffset(1.5);
    ROOT.gStyle.SetTitleYOffset(1.75);
    
    # Viridis palette
    stops = array('d', [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000])
    red   = array('d', [26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.])
    green = array('d', [9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.])
    blue  = array('d', [30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.])
    ROOT.TColor.CreateGradientColorTable(9, stops, red, green, blue, 255)
    #ROOT.gStyle.SetPalette(53)
    
    tree = ROOT.TChain('tjsev')
    tree.AddFile(fileName)
    totalEntries = tree.GetEntries()
    
    #initial histogram
    nbins   = 20
    lowbin  = 30
    highbin = 230
    print("Starting with bin width " + str((highbin-lowbin)/nbins))
    
    labels = ";generated;reconstructed"
    
    h = ROOT.TH2F("", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    
    fillHist(h, tree)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    h.SetMinimum(-1e-10)
    h.Draw("colz")
    c.Print("h.eps")
    
    indices  = []
    diagonal = []
    gensums  = []
    recosums = []
    
    for g in range(1, h.GetNbinsX()+1):
        gensum = 0
        indices.append(g)
        diagonal.append(h.GetBinContent(g, g))
        for r in range(1, h.GetNbinsY()+1):
            gensum += h.GetBinContent(g, r)
        gensums.append(gensum)
    
    for r in range(1, h.GetNbinsY()+1):
        recosum = 0
        for g in range(0, h.GetNbinsX()+1):
            recosum += h.GetBinContent(g, r)
        recosums.append(recosum)

    #print diagonal
    #print gensums
    #print recosums
    
    purities    = []
    stabilities = []
    
    for i in range(len(gensums)):
        purity    = -1
        stability = -1
        if gensums[i] > 0:
            purity = diagonal[i]/gensums[i]
        purities.append(purity)
        if recosums[i] > 0:
            stability = diagonal[i]/recosums[i]
        stabilities.append(stability)
        
    #print purities
    #print stabilities
    
    #splitForEqualPurity(h, gensums, recosums, indices)
    bins = splitForMinPurity(h, gensums, recosums, indices)
    
    print(bins)
    
    hnorm = ROOT.TH2F("", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    for g in range(1, h.GetNbinsX()+1):
        for r in range(1, h.GetNbinsY()+1):
            if gensums[g-1] > 0:
                hnorm.SetBinContent(g, r, h.GetBinContent(g, r)/gensums[g-1])
            else: hnorm.SetBinContent(g, r, 0)
    
    hnorm.SetMinimum(-1e-10)
    hnorm.Draw("colz")
    c.Print("hnorm.eps")
    
    bins2 = []
    for i in range(len(bins)-1):
        bins2.append(bins[i])
        bins2.append((bins[i]+bins[i+1])/2.)
    bins2.append(bins[-1])
    
    print(bins2)
    
    binArray = array('d', bins)
    bin2Array = array('d', bins) # put bins2 for reco bin split
    
    hopt = ROOT.TH2F("", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    fillHist(hopt, tree)
    
    hopt.SetMinimum(-1e-10)
    hopt.Draw("colz")
    c.Print("hopt.eps")
    
    optgensums = []
    for g in range(1, hopt.GetNbinsX()+1):
        gensum = 0
        for r in range(1, hopt.GetNbinsY()+1):
            gensum += hopt.GetBinContent(g, r)
        optgensums.append(gensum)

    hoptnorm = ROOT.TH2F("", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    for g in range(1, hopt.GetNbinsX()+1):
        for r in range(1, hopt.GetNbinsY()+1):
            if optgensums[g-1] > 0:
                hoptnorm.SetBinContent(g, r, hopt.GetBinContent(g, r)/optgensums[g-1])
            else: hoptnorm.SetBinContent(g, r, 0)
    
    hoptnorm.SetMinimum(-1e-10)
    hoptnorm.Draw("colz")
    c.Print("hoptnorm.eps")



def fillHist(h, tree):
    for event in tree:
        if event.gen_sel*event.reco_sel == -1: continue
        for j in range(event.nj):
            if event.j_gj[j] >= 0:
                h.Fill(event.gj_pt[event.j_gj[j]], event.j_pt[j])



def splitForMinPurity(h, gensums, recosums, indices):
    for i in range(len(gensums)):
        if len(gensums) <= 1: break
        if sum(gensums[:i]) == 0: continue
        if sum(gensums[i:]) == 0: continue
        if len(recosums) <= 1: break
        if sum(recosums[:i]) == 0: continue
        if sum(recosums[i:]) == 0: continue
        sumdiagonal1 = 0
        for j in range(i):
            for k in range(i):
                sumdiagonal1 += h.GetBinContent(indices[j], indices[k])
        sumpurity1 = sumdiagonal1/sum(gensums[:i])
        sumstability1 = sumdiagonal1/sum(recosums[:i])
        
        threshold = 0.7
        
        if sumpurity1 > threshold and sumstability1 > threshold:
            
            print(h.GetXaxis().GetBinLowEdge(indices[0]), h.GetXaxis().GetBinUpEdge(indices[i-1]), sumpurity1, sumstability1)
            
            return [h.GetXaxis().GetBinLowEdge(indices[0])] + splitForMinPurity(h, gensums[i:], recosums[i:], indices[i:])
            
    return [h.GetXaxis().GetBinUpEdge(indices[-1])]
            
            #return True



def splitForEqualPurity(h, gensums, recosums, indices):
    for i in range(len(gensums)):
        if len(gensums) <= 1: break
        if sum(gensums[:i]) == 0: continue
        if sum(gensums[i:]) == 0: continue
        if len(recosums) <= 1: break
        if sum(recosums[:i]) == 0: continue
        if sum(recosums[i:]) == 0: continue
        sumdiagonal1 = 0
        for j in range(i):
            for k in range(i):
                sumdiagonal1 += h.GetBinContent(indices[j], indices[k])
        sumpurity1 = sumdiagonal1/sum(gensums[:i])
        sumstability1 = sumdiagonal1/sum(recosums[:i])
        sumdiagonal2 = 0
        for j in range(i,len(gensums)):
            for k in range(i,len(gensums)):
                sumdiagonal2 += h.GetBinContent(indices[j], indices[k])
        sumpurity2 = sumdiagonal2/sum(gensums[i:])
        sumstability2 = sumdiagonal2/sum(recosums[i:])
        #print sumpurity1
        #print sumpurity2
        
        threshold = 0.4
        
        if sumpurity1 > sumpurity2:
            if sumpurity1 < threshold or sumpurity2 < threshold or sumstability1 < threshold or sumstability2 < threshold: return False
            
            success1 = splitForPurity(h, gensums[:i], recosums[:i], indices[:i]) 
            success2 = splitForPurity(h, gensums[i:], recosums[i:], indices[i:])
            
            if not success1:
                print(h.GetXaxis().GetBinLowEdge(indices[0]), h.GetXaxis().GetBinUpEdge(indices[i-1]), sumpurity1, sumstability1)
            if not success2:
                print(h.GetXaxis().GetBinLowEdge(indices[i]), h.GetXaxis().GetBinUpEdge(indices[-1]), sumpurity2, sumstability2)
            
            return True
            
            #return indices[i], sumpurity1, sumstability1, sumpurity2, sumstability2
        
    
"""
steer
"""
def main():
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--input',
                            dest='input',   
                            default='analysis.root',
                            help='input file [default: %default]')
    parser.add_option('--jobs',
                            dest='jobs', 
                            default=1,
                            type=int,
                            help='# of jobs to process in parallel the trees [default: %default]')
    parser.add_option('--only',
                            dest='only', 
                            default='',
                            type='string',
                            help='csv list of tags to process')
    parser.add_option('-o', '--output',
                            dest='output', 
                            default='analysis',
                            help='Output directory [default: %default]')
    parser.add_option('-q', '--queue',
                            dest='queue',
                            default='local',
                            help='Batch queue to use [default: %default]')
    (opt, args) = parser.parse_args()

          #ROOT.FWLiteEnabler.enable() 
    os.system('mkdir -p %s' % opt.output)

    optimize(opt.input)
        

if __name__ == "__main__":
	sys.exit(main())
