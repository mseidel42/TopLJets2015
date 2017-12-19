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
    ROOT.TH1.SetDefaultSumw2(True)
    
    treeBF = ROOT.TChain('tjsev')
    treeBF.Add('analysisBCDEF.root')

    treeGH = ROOT.TChain('tjsev')
    treeGH.Add('analysisGH.root')

    treeBF.Draw("j_mult_charged:gj_mult_charged >> h_bf(20,0,20,20,0,20)", "", "colz,text")
    h_bf = ROOT.gDirectory.Get("h_bf")

    treeGH.Draw("j_mult_charged:gj_mult_charged >> h_gh(20,0,20,20,0,20)", "", "colz,text")
    h_gh = ROOT.gDirectory.Get("h_gh")

    h_gh.Divide(h_bf)
    
    for i in range(h_gh.GetNbinsX()+2):
        for j in range(h_gh.GetNbinsY()+2):
            content = h_gh.GetBinContent(i, j)
            error = h_gh.GetBinError(i, j)
            if content != 0. and error != 0:
                h_gh.SetBinContent(i, j, (content - 1.)/error)

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
    ROOT.gStyle.SetPaintTextFormat('+.2f') 
    c = ROOT.TCanvas('c','c',500,500)
    c.SetRightMargin(0.15)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    c.cd()
    
    h_gh.SetMarkerSize(0.5)
    h_gh.GetZaxis().SetRangeUser(-1., 1.)
    h_gh.Draw('colz,text')
        
    c.Print('compareRunsMigration.pdf')    


if __name__ == "__main__":
	sys.exit(main())
