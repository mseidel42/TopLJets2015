import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt
from array import *
import random
import numpy
import copy
import pickle
import numpy as np

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS} #it{Preliminary}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('', '--inDirToys',       dest='inDirToys' ,      help='input toy directory',                default='unfolding/toys',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/covariance',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=35922.,              type=float)
    parser.add_option('--obs', dest='obs',  default='lowcornew', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='incl', help='flavor [default: %default]')
    parser.add_option('-r', '--reco', dest='reco', default='charged', help='Use charged/puppi/all particles [default: %default]')
    parser.add_option('-g', '--generator', dest='generator', default='pythia8', help='Parton shower generator [default: %default]')
    (opt, args) = parser.parse_args()
    
    observables = ["mult", "width", "ptd", "ptds", "ecc", "tau21", "tau32", "tau43", "zg", "zgxdr", "zgdr", "ga_width", "ga_lha", "ga_thrust", "c1_02", "c1_00", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "n2_b1", "n3_b1", "m2_b2", "n2_b2", "n3_b2", "nsd"]
    
    nice_observables_tex = {"mult": "N", "width": "$\\lambda_{1}^{1}$ (width)", "ptd": "$p_{T}D$", "ptds": "$p_{T}D^{s}$", "ecc": "$\\varepsilon$", "tau21": "$\\tau_{21}$", "tau32": "$\\tau_{32}$", "tau43": "$\\tau_{43}$", "zg": "$z_{g}$", "zgxdr": "$z_{g} \\times \\Delta R$", "zgdr": "$z_{g} \\Delta R$", "ga_width": "$\\lambda_{1}^{1}$ (width)", "ga_lha": "$\\lambda_{0.5}^{1}$ (LHA)", "ga_thrust": "$\\lambda_{2}^{1}$ (thrust)", "c1_02": "$C_{1}^{(0.2)}$", "c1_05": "$C_{1}^{(0.5)}$", "c1_10": "$C_{1}^{(1.0)}$", "c1_20": "$C_{1}^{(2.0)}$", "c2_02": "$C_{2}^{(0.2)}$", "c2_05": "$C_{2}^{(0.5)}$", "c2_10": "$C_{2}^{(1.0)}$", "c2_20":  "$C_{2}^{(2.0)}$", "c3_02": "$C_{3}^{(0.2)}$", "c3_05": "$C_{3}^{(0.5)}$", "c3_10": "$C_{3}^{(1.0)}$", "c3_20": "$C_{3}^{(2.0)}$", "m2_b1": "$M_{2}^{(1)}$", "n2_b1": "$N_{2}^{(1)}$", "n3_b1": "$N_{3}^{(1)}$", "m2_b2": "$M_{2}^{(2)}$", "n2_b2": "$N_{2}^{(2)}$", "n3_b2": "$N_{3}^{(2)}$"}
    
    nice_observables_root = {"mult": "#lambda_{0}^{0} (N)", "width": "#lambda_{1}^{1} (width)", "ptd": "#lambda_{0}^{2} (p_{T}^{d})", "ptds": "#lambda_{0}^{2}* (p_{T}^{d,}*)", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "#DeltaR_{g}", "ga_width": "#lambda_{1}^{1} (width)", "ga_lha": "#lambda_{0.5}^{1} (LHA)", "ga_thrust": "#lambda_{2}^{1} (thrust)", "c1_00": "C_{1}^{(0.0)}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1.0)}", "c1_20": "C_{1}^{(2.0)}", "c2_00": "C_{2}^{(0.0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1.0)}", "c2_20":  "C_{2}^{(2.0)}", "c3_00": "C_{3}^{(0.0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1.0)}", "c3_20": "C_{3}^{(2.0)}", "m2_b1": "M_{ 2}^{ (1)}", "n2_b1": "N_{ 2}^{ (1)}", "n3_b1": "N_{ 3}^{ (1)}", "m2_b2": "M_{ 2}^{ (2)}", "n2_b2": "N_{ 2}^{ (2)}", "n3_b2": "N_{ 3}^{ (2)}", "nsd": "n_{SD}"}
    
    generator_labels = {'pythia8': "POWHEG+PYTHIA 8", 'herwigpp': 'POWHEG+HERWIG++'}
    
    properSet = True
    if opt.obs in observables:
        observables = [opt.obs]
    elif opt.obs == 'lowcor':
        observables = ["ptds", "ecc", "zg", "zgdr", "tau43"]
        properSet = False
    elif opt.obs == 'c1all':
        observables = ["c1_00", "c1_02", "c1_05", "c1_10", "c1_20"]
        properSet = False
    elif opt.obs == 'c1pert':
        observables = ["c1_05", "c1_10", "c1_20"]
        properSet = False
    elif opt.obs == 'lowcornew':
        observables = ["ga_width", "ecc", "zg", "tau43"]
        properSet = False
    elif opt.obs == 'lambda':
        observables = ['mult', 'ptds', 'ga_lha', 'ga_width', 'ga_thrust']
        properSet = False
    
    styles = {'ptds': [ROOT.kGreen+1, 3],
              'ecc': [ROOT.kGreen+1, 3],
              'tau43': [ROOT.kAzure+1, 4],
              'zg': [ROOT.kViolet+1, 8],
              'zgdr': [ROOT.kMagenta+1, 6],
              'c1_00': [ROOT.kCyan+1, 2],
              'c1_02': [ROOT.kGreen+1, 3],
              'c1_05': [ROOT.kAzure+1, 4],
              'c1_10': [ROOT.kMagenta+1, 6],
              'c1_20': [ROOT.kViolet+1, 8],
              'ga_lha': [ROOT.kAzure+1, 4],
              'ga_width': [ROOT.kMagenta+1, 6],
              'ga_thrust': [ROOT.kViolet+1, 8],
              'mult': [ROOT.kCyan+1, 2],
             }
    defaultStyle = [ROOT.kBlack, 1]

    modelsToTest = []
    if (opt.generator == 'pythia8'):        
        modelsToTest.append('pythia8_asfsr0.070_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.080_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.090_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.100_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.105_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.110_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.115_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.120_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.125_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.130_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.135_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.140_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.150_meon_crdefault')
        modelsToTest.append('pythia8_asfsr0.160_meon_crdefault')
    elif (opt.generator == 'herwigpp'):
        modelsToTest.append('herwigpp_asfsr0.100_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.110_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.115_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.120_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.125_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.130_meon_crdefault')
    else:
        return
    
    x0 = []
    for model in modelsToTest:
        x0.append(float('0.'+model.split('.')[1][0:3]))
    
    fIn = ROOT.TFile.Open('unfolding/result/%s_%s_%s_result.root'%(opt.obs, opt.reco, opt.flavor))
    href = fIn.Get('FSRDownGen')
    
    y0 = []
    ymax = 0.
    for model in modelsToTest:
        htest = fIn.Get('MC13TeV_TTJets_%s_gen'%model)
        y0.append(href.Chi2Test(htest, 'NORM,CHI2'))
        if (y0[-1] > ymax):
            ymax = y0[-1]
    
    gr_sum = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0))
    #gr_sum.SetBit(ROOT.TGraph.kIsSortedX)
    gr_sum.SetLineColor(ROOT.kWhite)
    gr_sum.SetTitle()
    gr_sum.GetXaxis().SetTitle('#alpha_{s}^{FSR}(m(Z))')
    gr_sum.GetXaxis().SetTitleOffset(1.)
    gr_sum.GetXaxis().SetTitleSize(0.045)
    gr_sum.GetXaxis().SetLabelSize(0.04)
    gr_sum.GetYaxis().SetTitle('#chi^{2}')
    gr_sum.GetYaxis().SetTitleOffset(1.3)
    gr_sum.GetYaxis().SetTitleSize(0.045)
    gr_sum.GetYaxis().SetLabelSize(0.04)
    gr_sum.GetXaxis().SetRangeUser(x0[0], x0[-1])
    gr_sum.GetYaxis().SetRangeUser(0., 100)
    if opt.obs == 'ga_width':
        gr_sum.GetXaxis().SetRangeUser(0.115, 0.13)
        gr_sum.GetYaxis().SetRangeUser(0., 40)
    
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
    inix = 0.15
    iniy = 0.9 - 0.05*len(observables)
    if properSet: iniy -= 2*0.05
    legend = ROOT.TLegend(inix,iniy,inix+0.45,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    
    gr_sum.Draw('AC')
    #gr.Fit('pol2', 'em')
    
    x = []
    y = []
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        x.append(xval)
        y.append(yval)
    grfine = ROOT.TGraph(len(x), array('d',x) ,array('d',y))
    if properSet: grfine.Draw('same')
    grfine.SetLineColor(ROOT.kGray)
    grfine.SetLineWidth(4)
    if opt.obs in styles:
        grfine.SetLineColor(styles[opt.obs][0])
        grfine.SetLineStyle(styles[opt.obs][1])
    
    ymin = min(y)
    
    x2x = []
    y2x = []
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        if (yval > ymin * 2.0): continue
        x2x.append(xval)
        y2x.append(yval)
    gr_unc_2x = ROOT.TGraph(len(x2x), array('d',x2x) ,array('d',y2x))
    gr_unc_2x.SetLineColor(ROOT.kOrange+1)
    gr_unc_2x.SetLineWidth(4)
    gr_unc_2x.SetFillColor(ROOT.kOrange-9)
    gr_unc_2x.SetFillStyle(3854)
    #if properSet:
    #    gr_unc_2x.Draw('B')
    #    gr_unc_2x.Draw('C')
    
    x1 = []
    y1 = []
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        if (yval > ymin + 1): continue
        x1.append(xval)
        y1.append(yval)
    gr_unc_1 = ROOT.TGraph(len(x1), array('d',x1) ,array('d',y1))
    gr_unc_1.SetLineColor(ROOT.kRed+1)
    gr_unc_1.SetLineWidth(8)
    gr_unc_1.SetFillColor(ROOT.kRed-9)
    gr_unc_1.SetFillStyle(3254)
    if properSet:
        gr_unc_1.Draw('B')
        gr_unc_1.Draw('C')
    
    asfsr = x[y.index(min(y))]
    
    dummy = ROOT.TH1F("", "", 0, 0, 0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    if properSet:
        #legend.AddEntry(grfine, "Sum", "l")
        legend.AddEntry(gr_unc_1, "#delta #chi^{2} = 1", "lf")
        #legend.AddEntry(gr_unc_2x, "#delta #chi^{2} =  #chi^{2}_{min}", "lf")
        legend.AddEntry(dummy, "", "")
    
    
    gr = {}
    for obs in observables:
        #print(obs, 'chi2 values', y0[obs])
        gr[obs] = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0))
        #gr[obs].SetBit(ROOT.TGraph.kIsSortedX)
        gr[obs].SetFillColor(ROOT.kWhite)
        gr[obs].SetMarkerColor(styles.get(obs, defaultStyle)[0])
        gr[obs].SetMarkerStyle(5)
        gr[obs].SetMarkerSize(1.5)
        gr[obs].SetLineColor(styles.get(obs, defaultStyle)[0])
        gr[obs].SetLineStyle(styles.get(obs, defaultStyle)[1])
        gr[obs].SetLineWidth(4)
        gr[obs].Draw('p,same')
        legend.AddEntry(gr[obs], nice_observables_root[obs], "pl")
    
    legend.Draw()
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.041)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.66,0.97,'#scale[1.0]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    txt.SetTextAlign(32)
    txt.DrawLatex(0.90,0.91,'#scale[1.2]{%s}'%cmsLabel)
    txt.DrawLatex(0.90,0.85, '#scale[1.0]{'+generator_labels[opt.generator]+'}')
    txt.DrawLatex(0.90,0.8, '#scale[1.0]{'+opt.flavor+' jets}')
    txt.DrawLatex(0.90,0.75, '#scale[1.0]{'+opt.reco+' particles}')
    if properSet:
        txt.SetTextColor(ROOT.kRed+1)
        txt.DrawLatex(0.90,0.675, "#scale[1.0]{#alpha_{s}^{FSR}(m_{Z}) = %.4f_{%+.4f}^{%+.4f}}"%(asfsr, x1[0] - asfsr, x1[-1] - asfsr))
    
    ROOT.gPad.RedrawAxis()
    
    c.Print('fit/fitGenAlphaS_'+opt.obs+'_'+opt.generator+'_'+opt.flavor+'_'+opt.reco+'.pdf')
    c.Print('fit/fitGenAlphaS_'+opt.obs+'_'+opt.generator+'_'+opt.flavor+'_'+opt.reco+'.png')
    
    print('Best fit asfsr = %.4f'%(asfsr))
    print('dChi2 = 1 \n  envelope %.4f, %.4f \n  errors %.4f, %.4f'%(x1[0], x1[-1], x1[0] - asfsr, x1[-1] - asfsr))
    print('dChi2 = min(Chi2) = %.4f \n  envelope %.4f, %.4f \n  errors %.4f, %.4f'%(ymin, x2x[0], x2x[-1], x2x[0] - asfsr, x2x[-1] - asfsr))

def normalizeAndDivideByBinWidth(hist):
    hist.Scale(1./hist.Integral())
    for i in range(1, hist.GetNbinsX()+1):
        hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
        hist.SetBinError  (i, hist.GetBinError(i)  /hist.GetBinWidth(i))
    return hist
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

