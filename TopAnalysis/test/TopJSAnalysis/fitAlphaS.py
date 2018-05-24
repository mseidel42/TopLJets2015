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
    
    cmsLabel='#bf{CMS}'
    
    #configuration
    usage = 'usage: %prog [opt]'
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
    parser.add_option('', '--cmw', dest='cmw', default=False, action='store_true', help='Use CMW 2-loop (Pythia 8) [default: %default]')
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
    
    unsummedChi2 = pickle.load(open("unsummedChi2_"+opt.reco+".pkl", "rb"))

    modelsToTest = []
    if (opt.generator == 'pythia8'):
        #modelsToTest.append('pythia8_asfsr0.070_meon_crdefault')
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
        #modelsToTest.append('pythia8_asfsr0.150_meon_crdefault')
        #modelsToTest.append('pythia8_asfsr0.160_meon_crdefault')
    elif (opt.generator == 'herwigpp'):
        modelsToTest.append('herwigpp_asfsr0.100_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.110_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.115_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.120_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.125_meon_crdefault')
        modelsToTest.append('herwigpp_asfsr0.130_meon_crdefault')
    else:
        return

    cmw_models = {}
    if opt.cmw:
        cmw_as = ['0.080', '0.090', '0.100', '0.105', '0.110', '0.115', '0.120', '0.125', '0.130', '0.135', '0.140']
        cmw_scale = ['0.5', '1.0', '2.0']
        for scale in cmw_scale:
            cmw_models[scale] = []
            for alphas in cmw_as:
                cmw_models[scale].append('pythia8_asfsr%s_scale%s_CMW_2loop'%(alphas, scale))
    else:
        cmw_models['1.0'] = modelsToTest
    
    varModel = [['nominalGen'], #bias test
                ['evtgen'],
                ['m171v5', 'm173v5'],
                ['isrup', 'isrdn'],
                ['hdampup', 'hdampdn'],
                ['ueup', 'uedn'],
                ['erdON'],
                ['qcdBased', 'gluonMove'],
                ['pythia8_asfsr0.1365_meon_croff_flavrope'],
                ['bfrag_up', 'bfrag_down'], # b frag Bowler-Lund up/down
                ['bfrag_peterson'], # b frag Peterson
                ['slbr_up', 'slbr_down'], # B hadron semilep BR
                ['top_pt'], # top pt reweighting
                ['id1005muR2muF2hdampmt272.7225', 'id1009muR0.5muF0.5hdampmt272.7225'], # muF+muR
                ['id3001PDFset13100', 'id4001PDFset25200'], # PDF
               ]
    if not opt.cmw:
        varModel.append(['fsrup', 'fsrdn'])
    
    # START
    x0 = []
    for model in modelsToTest:
        x0.append(float('0.'+model.split('.')[1][0:3]))
    
    y0 = {'sum': []}
    ymax = 0.
    for obs in observables:
        y0[obs] = []
    for model in cmw_models['1.0']:
        y0['sum'].append(0.)
        for obs in observables:
            chi2ndf = unsummedChi2[obs]['data'][model][opt.flavor]
            y0['sum'][-1] += chi2ndf
            y0[obs].append(chi2ndf)
        if (y0['sum'][-1] > ymax):
            ymax = y0['sum'][-1]
    
    gr_sum = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0['sum']))
    
    # GET BEST ASFSR
    x = []
    y = []
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        x.append(xval)
        y.append(yval)
    ymin = max(0, min(y))
    asfsr = x[y.index(min(y))]
    print('ymin,asfsr', ymin, asfsr)
    
    if opt.cmw:
        # Get FSR scale uncertainties
        asfsr_scale = {}
        for scale in cmw_scale:
            y0_syst = {'sum': []}
            for model in cmw_models[scale]:
                y0_syst['sum'].append(0.)
                for obs in observables:
                    chi2ndf = unsummedChi2[obs]['data'][model][opt.flavor]
                    y0_syst['sum'][-1] += chi2ndf
            
            gr_syst = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0_syst['sum']))
            
            x_syst = []
            y_syst = []
            for xval in np.arange(x0[0], x0[-1], 0.000001):
                yval = gr_syst.Eval(xval, 0, 'S')
                x_syst.append(xval)
                y_syst.append(yval)
            ymin_syst = max(0, min(y_syst))
            asfsr_scale[scale] = x[y_syst.index(min(y_syst))]
        print(asfsr_scale)
        for a in asfsr_scale.iteritems():
            print(a[1] - asfsr_scale['1.0'])
    
    # Let's get full theory uncertainties
    ref = 0.1365
    up2 = 0.
    dn2 = 0.
    for var in varModel:
        for vardir in var:
            y0_syst = {'sum': []}
            for model in modelsToTest:
                y0_syst['sum'].append(0.)
                for obs in observables:
                    chi2ndf = unsummedChi2[obs][vardir][model][opt.flavor]
                    y0_syst['sum'][-1] += chi2ndf
            
            gr_syst = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0_syst['sum']))
            
            x_syst = []
            y_syst = []
            for xval in np.arange(x0[0], x0[-1], 0.000001):
                yval = gr_syst.Eval(xval, 0, 'S')
                x_syst.append(xval)
                y_syst.append(yval)
            ymin_syst = max(0, min(y_syst))
            asfsr_syst = x[y_syst.index(min(y_syst))]
            
            if asfsr_syst > ref:
                up2 += (asfsr_syst-ref)**2
            if asfsr_syst < ref:
                dn2 += (asfsr_syst-ref)**2
            print(vardir, '%.3f'%(asfsr_syst-ref))
    
    up2_total = up2
    dn2_total = dn2
    if opt.cmw:
        up2_total = up2 + (asfsr_scale['2.0'] - asfsr)**2
        dn2_total = dn2 + (asfsr_scale['0.5'] - asfsr)**2
    
    print('syst up   = %.3f'%sqrt(up2_total))
    print('syst down = %.3f'%sqrt(dn2_total))
    print('--> asfsr = %.3f + %.3f (stat) +/- %.3f (syst)')
    
    # DRAWING
    #gr_sum.SetBit(ROOT.TGraph.kIsSortedX)
    gr_sum.SetLineColor(ROOT.kWhite)
    gr_sum.SetTitle()
    gr_sum.GetXaxis().SetTitle('#alpha_{S}^{FSR}(m(Z))')
    if opt.cmw:
        gr_sum.GetXaxis().SetTitle('#alpha_{S,CMW}^{FSR}(m(Z))')
    gr_sum.GetXaxis().SetTitleOffset(1.)
    gr_sum.GetXaxis().SetTitleSize(0.045)
    gr_sum.GetXaxis().SetLabelSize(0.04)
    gr_sum.GetYaxis().SetTitle('#chi^{2}')
    gr_sum.GetYaxis().SetTitleOffset(1.3)
    gr_sum.GetYaxis().SetTitleSize(0.045)
    gr_sum.GetYaxis().SetLabelSize(0.04)
    gr_sum.GetXaxis().SetRangeUser(x0[0], x0[-1])
    gr_sum.GetYaxis().SetRangeUser(0, 750)
    #if opt.obs in ['ga_width', 'zgdr']:
    #    gr_sum.GetXaxis().SetRangeUser(0.08, 0.16)
    #    gr_sum.GetYaxis().SetRangeUser(0, 400)
    
    c = ROOT.TCanvas('c','c',500,500)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.1)
    c.SetTopMargin(0.05)
    
    inix = 0.15
    iniy = 0.9 - 0.03*(len(observables)+6)
    if properSet: iniy -= 2*0.05
    legend = ROOT.TLegend(inix,iniy,inix+0.45,0.9)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    
    gr_sum.Draw('AC')
    #gr.Fit('pol2', 'em')
    
    
    grfine = ROOT.TGraph(len(x), array('d',x) ,array('d',y))
    if properSet: grfine.Draw('same')
    grfine.SetLineColor(ROOT.kGray)
    grfine.SetLineWidth(4)
    if opt.obs in styles:
        grfine.SetLineColor(ROOT.kGray+3)
        grfine.SetLineStyle(styles[opt.obs][1])
    
    # get exp uncertainty
    x1 = []
    y1 = []
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        if (yval > ymin + 1): continue
        x1.append(xval)
        y1.append(yval)
    up2_total += (x1[-1] - asfsr)**2
    dn2_total += (x1[0] - asfsr)**2
    
    # draw total uncertainty
    x2x = []
    y2x = []
    for xval in np.arange(asfsr-sqrt(dn2_total), asfsr+sqrt(up2_total), 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        x2x.append(xval)
        y2x.append(yval)
    gr_unc_2x = ROOT.TGraph(len(x2x), array('d',x2x) ,array('d',y2x))
    gr_unc_2x.SetLineWidth(0)
    gr_unc_2x.SetFillColor(ROOT.kGreen+2)
    gr_unc_2x.SetFillStyle(1001)
    if properSet:
        gr_unc_2x.Draw('B')
    
    # draw fsr uncertainty (CMW)
    x2f = []
    y2f = []
    if opt.cmw:
        for xval in np.arange(asfsr_scale['0.5'], asfsr_scale['2.0'], 0.000001):
            yval = gr_sum.Eval(xval, 0, 'S')
            x2f.append(xval)
            y2f.append(yval)
        gr_unc_2f = ROOT.TGraph(len(x2f), array('d',x2f) ,array('d',y2f))
        gr_unc_2f.SetLineWidth(0)
        gr_unc_2f.SetFillColor(ROOT.kGreen-10)
        gr_unc_2f.SetFillStyle(1001)
        if properSet:
            gr_unc_2f.Draw('B')
    
    # draw exp uncertainty
    gr_unc_1 = ROOT.TGraph(len(x1), array('d',x1) ,array('d',y1))
    gr_unc_1.SetLineWidth(0)
    gr_unc_1.SetFillColor(ROOT.kAzure+1)
    gr_unc_1.SetFillStyle(1001)
    if properSet:
        gr_unc_1.Draw('B')
        grfine.Draw('same')
    
    # draw model uncertainty
    x2m = []
    y2m = []
    for xval in np.arange(asfsr-sqrt(dn2), asfsr+sqrt(up2), 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        x2m.append(xval)
        y2m.append(yval)
    gr_unc_2m = ROOT.TGraph(len(x2m), array('d',x2m) ,array('d',y2m))
    gr_unc_2m.SetLineWidth(0)
    gr_unc_2m.SetFillColor(ROOT.kBlack)
    gr_unc_2m.SetFillStyle(3254)
    if properSet:
        gr_unc_2m.Draw('B')
    
    gr = {}
    for obs in observables:
        #print(obs, 'chi2 values', y0[obs])
        gr[obs] = ROOT.TGraph(len(x0), array('d',x0) ,array('d',y0[obs]))
        #gr[obs].SetBit(ROOT.TGraph.kIsSortedX)
        gr[obs].SetFillColor(ROOT.kWhite)
        gr[obs].SetMarkerColor(ROOT.kGray+3)
        gr[obs].SetMarkerStyle(5)
        gr[obs].SetMarkerSize(1.5)
        gr[obs].SetLineColor(ROOT.kGray+3)
        gr[obs].SetLineStyle(styles.get(obs, defaultStyle)[1])
        gr[obs].SetLineWidth(4)
        gr[obs].Draw('p,same')
        legend.AddEntry(gr[obs], nice_observables_root[obs], "pl")

    dummy = ROOT.TH1F("", "", 0, 0, 0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    if properSet:
        #legend.AddEntry(grfine, "Sum", "l")
        legend.AddEntry(dummy, "", "")
        legend.AddEntry(dummy, "Uncertainties", "")
        legend.AddEntry(gr_unc_1, "Experimental", "lf")
        legend.AddEntry(gr_unc_2m, "Model", "lf")
        if opt.cmw: legend.AddEntry(gr_unc_2f, "FSR scale", "lf")
        legend.AddEntry(gr_unc_2x, "Total", "lf")
    
    
    legend.Draw()
    
    box = ROOT.TLegend(0.5, 0.625, 0.9, 0.875)
    box.SetFillStyle(1001)
    box.SetFillColor(ROOT.kWhite)
    box.SetLineColor(ROOT.kWhite)
    box.Draw()

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
        #txt.SetTextColor(ROOT.kRed+1)
        if opt.cmw:
            txt.DrawLatex(0.90,0.675, "#scale[1.0]{#alpha_{S,CMW}^{FSR}(m_{Z}) = %.3f_{%+.3f}^{ %+.3f}}"%(asfsr, x2x[0] - asfsr, x2x[-1] - asfsr))
        else:
            txt.DrawLatex(0.90,0.675, "#scale[1.0]{#alpha_{S}^{FSR}(m_{Z}) = %.3f_{%+.3f}^{ %+.3f}}"%(asfsr, x2x[0] - asfsr, x2x[-1] - asfsr))
    if opt.cmw:
        txt.DrawLatex(0.90,0.6, 'LO, CMW, 2-loop')
    
    ROOT.gPad.RedrawAxis()
    
    if opt.cmw:
        c.Print('fit/fitAlphaS_'+opt.obs+'_'+opt.generator+'_cmw_2-loop_'+opt.flavor+'_'+opt.reco+'.pdf')
        c.Print('fit/fitAlphaS_'+opt.obs+'_'+opt.generator+'_cmw_2-loop_'+opt.flavor+'_'+opt.reco+'.png')
    else:
        c.Print('fit/fitAlphaS_'+opt.obs+'_'+opt.generator+'_'+opt.flavor+'_'+opt.reco+'.pdf')
        c.Print('fit/fitAlphaS_'+opt.obs+'_'+opt.generator+'_'+opt.flavor+'_'+opt.reco+'.png')
    
    print('Best fit asfsr = %.4f'%(asfsr))
    print('dChi2 = 1 \n  envelope %.4f, %.4f \n  errors %.4f, %.4f'%(x1[0], x1[-1], x1[0] - asfsr, x1[-1] - asfsr))
    print('exp+model \n  envelope %.4f, %.4f \n  errors %.4f, %.4f'%(x2x[0], x2x[-1], x2x[0] - asfsr, x2x[-1] - asfsr))
    
    # What alpha_s do fsr scale up/down correspond to?
    chi2fsrdn = unsummedChi2[obs]['data']['fsrdn'][opt.flavor]
    chi2fsrup = unsummedChi2[obs]['data']['fsrup'][opt.flavor]
    found_down = False
    for xval in np.arange(x0[0], x0[-1], 0.000001):
        yval = gr_sum.Eval(xval, 0, 'S')
        if not found_down and yval < chi2fsrdn:
            found_down = True
            print('FSR down ~= %.4f'%xval)
        if found_down and yval > chi2fsrup:
            print('FSR up   ~= %.4f'%xval)
            break

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

