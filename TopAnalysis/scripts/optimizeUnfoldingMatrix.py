#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import time
from array import *
from math import sqrt

"""
Analysis loop
"""
def optimize(fileName, obs):
    
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.05);
    ROOT.gStyle.SetPadBottomMargin(0.13);
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.15)
    ROOT.gStyle.SetTitleXOffset(1.5);
    ROOT.gStyle.SetTitleYOffset(1.75);
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    
    # Viridis palette reversed + white
    stops = array('d', [0.0, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0000])
    red   = array('d', [26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255., 1., 1.])
    green = array('d', [9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255., 1., 1.])
    blue  = array('d', [30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255., 1., 1.])
    ROOT.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 255)
    #ROOT.gStyle.SetPalette(53)
    
    tree = ROOT.TChain('tjsev')
    tree.AddFile(fileName)
    totalEntries = tree.GetEntries()
    
    print("Observable: " + obs)
    
    #initial histogram
    # jet pt
    #nbins   = 20
    #lowbin  = 30
    #highbin = 230
    
    if (obs == "mult"):
        nbins   = 30
        lowbin  = 0
        highbin = 30
        label = "N_{ch}"
    if (obs == "width"):
        nbins   = 40
        lowbin  = 0
        highbin = 0.4
        label = "width"
    if (obs == "ptd"):
        nbins   = 50
        lowbin  = 0
        highbin = 1
        label = "p_{T}D"
    
    label1 = ";" + label
    labels = ";generated " + label + ";reconstructed " + label
    
    print("Starting with bin width " + str((highbin-lowbin)/float(nbins)))
    
    h = ROOT.TH2F("", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    
    fillHist(h, tree, obs)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    h.SetMinimum(-1e-10)
    h.Draw("colz")
    #h.Fit("pol1")
    c.Print("unfolding/" + obs + "_h.eps")
    
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
        if recosums[i] > 0:
            purity = diagonal[i]/recosums[i]
        purities.append(purity)
        if gensums[i] > 0:
            stability = diagonal[i]/gensums[i]
        stabilities.append(stability)
        
    #print purities
    #print stabilities
    
    #splitForEqualPurity(h, gensums, recosums, indices)
    #bins = splitForMinPurity(h, gensums, recosums, indices)
    bins = splitForMinSigma(h, 1.)
    
    c = ROOT.TCanvas('c', 'c', 500, 450)
    c.cd()
    
    print(bins)
    
    hnorm = ROOT.TH2F("", labels, nbins, lowbin, highbin, nbins, lowbin, highbin)
    for g in range(1, h.GetNbinsX()+1):
        for r in range(1, h.GetNbinsY()+1):
            if recosums[r-1] > 0:
                hnorm.SetBinContent(g, r, h.GetBinContent(g, r)/recosums[r-1])
            else: hnorm.SetBinContent(g, r, 0)
    
    hnorm.SetMinimum(-1e-10)
    #hnorm.SetMaximum(1)
    hnorm.Draw("colz")
    #hnorm.Fit("pol1")
    c.Print("unfolding/" + obs + "_hnorm.eps")
    
    # reco bin splitting
    divisor = 1
    bins2 = []
    for i in range(len(bins)-1):
        for j in range(divisor):
            bins2.append(bins[i] + abs(bins[i]-bins[i+1])/divisor*j)
    bins2.append(bins[-1])
    
    #print(bins2)
    
    # original reco bins
    bins3 = []
    for i in range(nbins+1):
        bins3.append(lowbin + (highbin-lowbin)/nbins*i)
    
    #print bins3
    
    binArray = array('d', bins)
    bin2Array = array('d', bins2) # put bins2 for reco bin split
    
    hopt = ROOT.TH2F("", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    fillHist(hopt, tree, obs)
    
    hopt.SetMinimum(-1e-10)
    hopt.Draw("colz")
    c.Print("unfolding/" + obs + "_hopt.eps")
    
    optgensums = []
    for g in range(1, hopt.GetNbinsX()+1):
        gensum = 0
        for r in range(1, hopt.GetNbinsY()+1):
            gensum += hopt.GetBinContent(g, r)
        optgensums.append(gensum)
    
    optrecosums = []
    for r in range(1, hopt.GetNbinsY()+1):
        recosum = 0
        for g in range(1, hopt.GetNbinsX()+1):
            recosum += hopt.GetBinContent(g, r)
        optrecosums.append(recosum)
    
    hoptpur  = ROOT.TH1F("", label1, len(binArray)-1, binArray)
    hoptsta  = ROOT.TH1F("", label1, len(binArray)-1, binArray)
    hoptnorm = ROOT.TH2F("", labels, len(binArray)-1, binArray, len(bin2Array)-1, bin2Array)
    for g in range(1, hopt.GetNbinsX()+1):
        for r in range(1, hopt.GetNbinsY()+1):
            purity    = 0
            if optrecosums[r-1] > 0:
              purity    = hopt.GetBinContent(g, r)/optrecosums[r-1]
            hoptnorm.SetBinContent(g, r, purity)
            if g == r:
              hoptpur.SetBinContent(g, purity)
              stability = 0
              if optgensums[g-1] > 0:
                stability = hopt.GetBinContent(g, r)/optgensums[g-1]
              hoptsta.SetBinContent(g, stability)
    
    hoptnorm.SetMinimum(-1e-10)
    #hoptnorm.SetMaximum(1)
    hoptnorm.SetMarkerColor(ROOT.kWhite)
    hoptnorm.SetMarkerSize(1.5)
    hoptnorm.Draw("colz,text")
    c.Print("unfolding/" + obs + "_hoptnorm.eps")
    
    hoptpur.GetYaxis().SetRangeUser(0., 1.)
    hoptpur.SetLineColor(ROOT.kRed+1)
    hoptpur.Draw()
    hoptsta.SetLineStyle(7)
    hoptsta.Draw("same")
    c.Print("unfolding/" + obs + "_hoptpursta.eps")


def fillHist(h, tree, obs):
    for event in tree:
        if event.gen_sel*event.reco_sel == -1: continue
        for j in range(event.nj):
            if event.j_gj[j] >= 0:
                #j_p4 = ROOT.TLorentzVector()
                #j_p4.SetPtEtaPhiM(event.j_pt[j], event.j_eta[j], event.j_phi[j], event.j_m[j])
                #
                #g = event.j_gj[j]
                #gj_p4 = ROOT.TLorentzVector()
                #gj_p4.SetPtEtaPhiM(event.gj_pt[g], event.gj_eta[g], event.gj_phi[g], event.gj_m[g])
                #
                ##h.Fill(event.gj_pt[event.j_gj[j]], event.j_pt[j]) # pt
                ##h.Fill(event.gj_m[event.j_gj[j]], event.j_m[j]) # m
                ##h.Fill(gj_p4.M(), j_p4.M()) # m
                #h.Fill(gj_p4.M()/gj_p4.E(), j_p4.M()/j_p4.E()) # m/E
                
                #ga = []
                #for i in range(len(event.gj_ga)): ga.append(event.gj_ga[i])
                #print ga
                
                c = 0 #0=charged, 1=charged+neutral, 2=puppi
                
                i = event.j_gj[j]
                if (obs == "mult"):  h.Fill(event.gj_ga[i*81 + 0+c],  event.j_ga[j*81 + 0+c]) # charged mult
                #if (obs == "mult2"):  h.Fill(event.gj_ga[i*81 + 0+c],  (event.j_ga[j*81 + 0+c] - 1.01747)/1.14334) # charged mult
                if (obs == "width"): h.Fill(event.gj_ga[i*81 + 36+c] * 0.4, event.j_ga[j*81 + 36+c] * 0.4) # charged width
                if (obs == "ptd"):   h.Fill(sqrt(event.gj_ga[i*81 + 18+c]), sqrt(event.j_ga[j*81 + 18+c])) # charged ptd



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
        sumpurity1 = sumdiagonal1/sum(recosums[:i])
        sumstability1 = sumdiagonal1/sum(gensums[:i])
        
        threshold = 0.5
        
        if sumpurity1 > threshold and sumstability1 > threshold:
            
            print(h.GetXaxis().GetBinLowEdge(indices[0]), h.GetXaxis().GetBinUpEdge(indices[i-1]), sumpurity1, sumstability1)
            
            return [h.GetXaxis().GetBinLowEdge(indices[0])] + splitForMinPurity(h, gensums[i:], recosums[i:], indices[i:])
            
    return [h.GetXaxis().GetBinUpEdge(indices[-1])]
            
            #return True

def splitForMinSigma(h, factor = 1):   
    slices = ROOT.TObjArray()
    h.FitSlicesX(0, 0, -1, 0, "QNR", slices)
    
    c2 = ROOT.TCanvas('c', 'c', 500, 450)
    c2.cd()
    
    slices[0].Draw()
    c2.Print("unfolding/slices0.eps")
    slices[1].GetYaxis().SetRangeUser(0, 30)
    slices[1].Draw()
    c2.Print("unfolding/slices1.eps")
    slices[2].GetYaxis().SetRangeUser(0, 10)
    slices[2].Draw()
    c2.Print("unfolding/slices2.eps")
    
    bins  = [0]
    exact = [0]
    
    for i in range(1, h.GetNbinsX()+1):
      mean  = slices[1].GetBinContent(i)
      sigma = slices[2].GetBinContent(i) * factor
      #print(bins[:-1])
      #print(mean)
      #print(sigma)
      if (mean - sigma) > exact[-1]:
        exact.append(mean+sigma)
        bins.append(h.GetXaxis().GetBinUpEdge((h.GetXaxis().FindBin(mean+sigma))))
    
    bins[-1] = h.GetXaxis().GetXmax()
    
    return bins


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
    parser.add_option('--obs',
                            dest='obs',   
                            default='mult',
                            help='observable [default: %default]')
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

    optimize(opt.input, opt.obs)
        

if __name__ == "__main__":
	sys.exit(main())
