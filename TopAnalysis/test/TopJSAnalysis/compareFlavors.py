import ROOT
ROOT.gROOT.SetBatch(True)
import optparse
import os,sys
import json
import re
from collections import OrderedDict
from math import sqrt
from array import *

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding/result',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=35922.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--reco', dest='reco',  default='charged', help='Jet shapes from charged or all particles [default: %default]')
    (opt, args) = parser.parse_args()

    #read lists of samples
    flavors = ['incl', 'bottom', 'light', 'gluon']
    # ColorBrewer 3-class Set2 (http://colorbrewer2.org/#type=qualitative&scheme=Set2&n=3)
    colors = {'incl': ROOT.kBlack, 'bottom': ROOT.kOrange+8, 'light': ROOT.kBlue-5, 'gluon': ROOT.kTeal-8}
    lightcolors = {'incl': ROOT.kGray, 'bottom': ROOT.kOrange-9, 'light': ROOT.kBlue-10, 'gluon': ROOT.kCyan-10}
    markers = {'incl': 20, 'bottom': 21, 'light': 22, 'gluon': 23}
    fills = {'incl': 1001, 'bottom': 3254, 'light': 3245, 'gluon': 3002}
    offsets = {'incl': 0, 'bottom': -1, 'light': 0, 'gluon': 1}
    infiles = {}
    hists = {}
    unchists = {}
    ratios = {}
    uncratios = {}
    pythiauncratios = {}
    pythiaratios = {}
    pythiafsrupratios = {}
    pythiafsrdownratios = {}
    #pythiafsrtuneratios = {}
    herwigratios = {}
    sherparatios = {}
    direratios = {}
    line = {}
    
    nice_observables_root = {"mult": "#lambda_{0}^{0} (N)", "width": "#lambda_{1}^{1} (width)", "ptd": "#lambda_{0}^{2} (p_{T}^{d})", "ptds": "#lambda_{0}^{2}* (p_{T}^{d,}*)  ", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "#DeltaR_{g}", "ga_width": "#lambda_{1}^{1} (width)", "ga_lha": "#lambda_{0.5}^{1} (LHA)", "ga_thrust": "#lambda_{2}^{1} (thrust)", "c1_00": "C_{1}^{(0)}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1)}", "c1_20": "C_{1}^{(2)}", "c2_00": "C_{2}^{(0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1)}", "c2_20":  "C_{2}^{(2)}", "c3_00": "C_{3}^{(0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1)}     ", "c3_20": "C_{3}^{(2)}", "m2_b1": "M_{ 2}^{ (1)}", "n2_b1": "N_{ 2}^{ (1)}", "n3_b1": "N_{ 3}^{ (1)}", "m2_b2": "M_{ 2}^{ (2)}", "n2_b2": "N_{ 2}^{ (2)}", "n3_b2": "N_{ 3}^{ (2)}", "nsd": "n_{SD}"}
    
    nice_observables_root_short = {"mult": "#lambda_{0}^{0}", "width": "#lambda_{1}^{1}", "ptd": "#lambda_{0}^{2}", "ptds": "#lambda_{0}^{2}*", "ecc": "#varepsilon", "tau21": "#tau_{21}", "tau32": "#tau_{32}", "tau43": "#tau_{43}", "zg": "z_{g}", "zgxdr": "z_{g} #times #DeltaR", "zgdr": "#DeltaR_{g}", "ga_width": "#lambda_{1}^{1}", "ga_lha": "#lambda_{0.5}^{1}", "ga_thrust": "#lambda_{2}^{1}", "c1_00": "C_{1}^{(0)}", "c1_02": "C_{1}^{(0.2)}", "c1_05": "C_{1}^{(0.5)}", "c1_10": "C_{1}^{(1)}", "c1_20": "C_{1}^{(2)}", "c2_00": "C_{2}^{(0)}", "c2_02": "C_{2}^{(0.2)}", "c2_05": "C_{2}^{(0.5)}", "c2_10": "C_{2}^{(1)}", "c2_20":  "C_{2}^{(2)}", "c3_00": "C_{3}^{(0)}", "c3_02": "C_{3}^{(0.2)}", "c3_05": "C_{3}^{(0.5)}", "c3_10": "C_{3}^{(1)}", "c3_20": "C_{3}^{(2)}", "m2_b1": "M_{ 2}^{ (1)}", "n2_b1": "N_{ 2}^{ (1)}", "n3_b1": "N_{ 3}^{ (1)}", "m2_b2": "M_{ 2}^{ (2)}", "n2_b2": "N_{ 2}^{ (2)}", "n3_b2": "N_{ 3}^{ (2)}", "nsd": "n_{SD}"}
    
    for flavor in flavors:
        infiles[flavor] = ROOT.TFile.Open('%s/%s_%s_%s_result.root'%(opt.inDir, opt.obs, opt.reco, flavor))
        
        hists[flavor] = infiles[flavor].Get('MC13TeV_TTJets_Unfolded').Clone()
        hists[flavor].SetTitle('')
        hists[flavor].SetMarkerColor(colors[flavor])
        hists[flavor].SetMarkerStyle(markers[flavor])
        hists[flavor].SetLineColor(colors[flavor])
        hists[flavor].SetFillColor(lightcolors[flavor])
        #hists[flavor].SetFillStyle(fills[flavor])
        
        unchists[flavor] = infiles[flavor].Get('dataUnfoldedSys').Clone()
        unchists[flavor].SetMarkerStyle(0)
        unchists[flavor].SetFillColor(colors[flavor])
        unchists[flavor].SetFillStyle(fills[flavor])
        unchists[flavor].SetLineColor(colors[flavor])
        
        pythiaratios[flavor] = infiles[flavor].Get('nominalGenRatio')
        pythiafsrupratios[flavor] = infiles[flavor].Get('FSRUpGenRatio')
        pythiafsrdownratios[flavor] = infiles[flavor].Get('FSRDownGenRatio')
        #pythiafsrtuneratios[flavor] = infiles[flavor].Get('tunedGenRatio')
        herwigratios[flavor] = infiles[flavor].Get('herwigGenRatio')
        sherparatios[flavor] = infiles[flavor].Get('sherpaGenRatio')
        direratios[flavor] = infiles[flavor].Get('direGenRatio')
        
        pythiauncratios[flavor] = infiles[flavor].Get('dataUnfoldedSysRatio').Clone()
        pythiauncratios[flavor].SetMarkerColor(colors[flavor])
        pythiauncratios[flavor].SetMarkerStyle(0)
        pythiauncratios[flavor].SetLineColor(colors[flavor])
        pythiauncratios[flavor].SetFillColor(lightcolors[flavor])
        #pythiauncratios[flavor].SetFillStyle(fills[flavor])
        
        line[flavor] = infiles[flavor].Get('line')
        line[flavor].SetLineColor(colors[flavor])
        line[flavor].SetMarkerColor(colors[flavor])
        line[flavor].SetMarkerStyle(markers[flavor])
        #for xbin in range(line[flavor].GetNbinsX()+1):
        #    line[flavor].SetBinError(xbin, 0.0001)
        
        if (colors[flavor] == ROOT.kBlack):
            hists[flavor].SetFillColor(ROOT.kGray)
            unchists[flavor].SetFillColor(ROOT.kGray)
            pythiauncratios[flavor].SetFillColor(ROOT.kGray)
    
    for flavor in flavors:
        ratios[flavor] = hists[flavor].Clone()
        ratios[flavor].Divide(hists['incl'])
        uncratios[flavor] = unchists[flavor].Clone()
        uncratios[flavor].Divide(hists['incl'])
        yrange = [-0.25,2.25]
        if opt.obs == 'zg': yrange = [0.75,1.25]
        if opt.obs == 'zgdr': yrange = [0.7,1.3]
        limitToRange(uncratios[flavor], yrange)
        limitToRange(pythiauncratios[flavor], yrange)
        if flavor == 'incl': continue
        shift(hists[flavor], offsets[flavor])
        shift(unchists[flavor], offsets[flavor])
        shift(ratios[flavor], offsets[flavor])
        shift(uncratios[flavor], offsets[flavor])
    
    #plot
    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas('c','c',500,500)
    c.SetBottomMargin(0.0)
    c.SetLeftMargin(0.0)
    c.SetTopMargin(0)
    c.SetRightMargin(0.00)
    c.cd()
    
    p1=ROOT.TPad('p1','p1',0.0,0.60,1.0,1.0)
    p1.SetRightMargin(0.05)
    p1.SetLeftMargin(0.12)
    p1.SetTopMargin(0.06*0.68/0.40)
    p1.SetBottomMargin(0.01)
    p1.Draw()
    p1.cd()
    
    unchists['incl'].SetTitle('')
    unchists['incl'].SetYTitle('1/N_{jet} dN_{jet} / d '+nice_observables_root_short[opt.obs])
    unchists['incl'].GetXaxis().SetTitleOffset(1.)
    unchists['incl'].GetXaxis().SetLabelOffset(1.)
    unchists['incl'].GetYaxis().SetTitleSize(0.0375*pow(0.68/0.40,2))
    unchists['incl'].GetYaxis().SetTitleOffset(1.4/pow(0.68/0.40,2))
    unchists['incl'].GetYaxis().SetLabelSize(0.08)
    unchists['incl'].GetYaxis().SetNdivisions(505)
    unchists['incl'].GetYaxis().SetRangeUser(0.0001, unchists['incl'].GetMaximum()*1.75)
    m = re.search('(.*) (\(.+\))', hists['incl'].GetXaxis().GetTitle())
    unchists['incl'].GetXaxis().SetTitle('')
    unchists['incl'].Draw('e2')
    unchists['bottom'].Draw('x0 e1,same')
    unchists['light' ].Draw('x0 e1,same')
    unchists['gluon' ].Draw('x0 e1,same')
    hists['incl'].Draw  ('e, same')
    hists['bottom'].Draw('x0 e, same')
    hists['light' ].Draw('x0 e, same')
    hists['gluon' ].Draw('x0 e, same')
    
    fliplegend = False
    flip_threshold = 0.4
    if (unchists['incl'].GetXaxis().GetBinCenter(unchists['incl'].GetMaximumBin()) > (unchists['incl'].GetXaxis().GetXmax() - unchists['incl'].GetXaxis().GetXmin())*flip_threshold):
        fliplegend = True
    inix = 0.6 if not fliplegend else 0.15
    
    legend = ROOT.TLegend(inix,0.25,inix+0.3,0.85)
    legend.SetLineWidth(0)
    legend.SetFillStyle(0)
    legend.AddEntry(hists['incl'], 'Inclusive jets',  'fepl')
    legend.AddEntry(hists['bottom'], 'Bottom jets',   'fepl')
    legend.AddEntry(hists['light'], 'Light-enriched', 'fepl')
    legend.AddEntry(hists['gluon'], 'Gluon-enriched', 'fepl')
    dummy = ROOT.TH1F("", "", 1, 0, 1)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    legend.AddEntry(dummy, '(p_{T} > 30 GeV)', '')
    legend.Draw()
    
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.10)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.66,0.95,'#scale[1.0]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
    inix = 0.15 if not fliplegend else 0.92
    if fliplegend: txt.SetTextAlign(32)
    txt.DrawLatex(inix,0.80,'#scale[1.2]{%s}'%cmsLabel)
    txt.DrawLatex(inix,0.71,'#scale[1.0]{t#bar{t} #rightarrow lepton+jets}')
    
    c.cd()
    p2 = ROOT.TPad('p2','p2',0.0,0.48,1.0,0.60)
    p2.Draw()
    p2.SetBottomMargin(0.02)
    p2.SetRightMargin(0.05)
    p2.SetLeftMargin(0.12)
    p2.SetTopMargin(0.00)
    p2.cd()
    
    uncratios['incl'].SetTitle('')
    uncratios['incl'].SetYTitle('#frac{flavor}{incl} ')
    uncratios['incl'].SetFillColor(ROOT.kGray)
    uncratios['incl'].GetYaxis().SetTitleSize(0.2*(0.32/0.2))
    uncratios['incl'].GetYaxis().SetTitleOffset(0.15)
    uncratios['incl'].GetYaxis().SetLabelSize(0.18*(0.32/0.2))
    uncratios['incl'].GetYaxis().SetRangeUser(yrange[0], yrange[1])
    uncratios['incl'].GetYaxis().SetNdivisions(503)
    
    uncratios['incl'].Draw  ('e2')
    uncratios['bottom'].Draw('x0 e1,same')
    uncratios['light' ].Draw('x0 e1,same')
    uncratios['gluon' ].Draw('x0 e1,same')
    ratios['incl'].Draw  ('e same')
    ratios['gluon' ].Draw('x0 e1 same')
    ratios['bottom'].Draw('x0 e1 same')
    ratios['light' ].Draw('x0 e1 same')
    
    c.cd()
    p3 = ROOT.TPad('p3','p3',0.0,0.32,1.0,0.44)
    p3.Draw()
    p3.SetBottomMargin(0.02)
    p3.SetRightMargin(0.05)
    p3.SetLeftMargin(0.12)
    p3.SetTopMargin(0.00)
    p3.cd()
    
    #pythiauncratios['bottom'].GetYaxis().SetTitleColor(colors['bottom'])
    pythiauncratios['bottom'].GetYaxis().SetTitle('#splitline{MC/data}{(bottom)}  ')
    pythiauncratios['bottom'].GetYaxis().SetTitleSize(0.25)
    pythiauncratios['bottom'].GetYaxis().SetTitleOffset(0.2)
    pythiauncratios['bottom'].GetYaxis().SetLabelSize(0.18*(0.32/0.2))
    pythiauncratios['bottom'].Draw('e2')
    line['bottom'].Draw('x0 e1 same')
    line['bottom'].Draw('same')
    pythiaratios['bottom'].Draw('same h')
    pythiafsrupratios['bottom'].Draw('SAME P X0 E1')
    pythiafsrdownratios['bottom'].Draw('SAME P X0 E1')
    #pythiafsrtuneratios['bottom'].Draw('same h')
    herwigratios['bottom'].Draw('same h')
    sherparatios['bottom'].Draw('same h')
    direratios['bottom'].Draw('same h')
    
    c.cd()
    plegend = ROOT.TPad('plegend','plegend',0.0,0.4415,1.0,0.48)
    plegend.Draw()
    plegend.SetBottomMargin(0.0)
    plegend.SetRightMargin(0.05)
    plegend.SetLeftMargin(0.12)
    plegend.SetTopMargin(0.00)
    plegend.cd()

    mclegend = ROOT.TLegend(0.12,0,0.95,1)
    mclegend.SetLineWidth(0)
    mclegend.SetFillStyle(0)
    mclegend.SetNColumns(6)
    mclegend.AddEntry(pythiaratios['bottom'], "POWHEG+PYTHIA 8 ", "pl")
    mclegend.AddEntry(pythiafsrupratios['bottom'], "FSR up", "p")
    mclegend.AddEntry(pythiafsrdownratios['bottom'], "down  ", "p")
    #mclegend.AddEntry(pythiafsrtuneratios['bottom'], "tuned", "l")
    mclegend.AddEntry(herwigratios['bottom'], "POWHEG+HERWIG 7  ", "pl")
    mclegend.AddEntry(sherparatios['bottom'], "SHERPA 2 ", "pl")
    mclegend.AddEntry(direratios['bottom'], "DIRE 2", "pl")
    mclegend.Draw()
    
    c.cd()
    p4 = ROOT.TPad('p4','p4',0.0,0.20,1.0,0.32)
    p4.Draw()
    p4.SetBottomMargin(0.02)
    p4.SetRightMargin(0.05)
    p4.SetLeftMargin(0.12)
    p4.SetTopMargin(0.00)
    p4.cd()
    
    #pythiauncratios['light'].GetYaxis().SetTitleColor(colors['light'])
    pythiauncratios['light'].GetYaxis().SetTitle('#splitline{MC/data}{  (light)}  ')
    pythiauncratios['light'].GetYaxis().SetTitleSize(0.25)
    pythiauncratios['light'].GetYaxis().SetTitleOffset(0.2)
    pythiauncratios['light'].GetYaxis().SetLabelSize(0.18*(0.32/0.2))
    pythiauncratios['light'].Draw('e2')
    line['light'].Draw('x0 e1 same')
    line['light'].Draw('same')
    pythiaratios['light'].Draw('same h')
    pythiafsrupratios['light'].Draw('SAME P X0 E1')
    pythiafsrdownratios['light'].Draw('SAME P X0 E1')
    #pythiafsrtuneratios['light'].Draw('same h')
    herwigratios['light'].Draw('same h')
    sherparatios['light'].Draw('same h')
    direratios['light'].Draw('same h')
    
    c.cd()
    p5 = ROOT.TPad('p5','p5',0.0,0.0,1.0,0.20)
    p5.Draw()
    p5.SetBottomMargin(0.4)
    p5.SetRightMargin(0.05)
    p5.SetLeftMargin(0.12)
    p5.SetTopMargin(0.01)
    p5.cd()
    
    pythiauncratios['gluon'].GetXaxis().SetTitle(nice_observables_root[opt.obs])
    pythiauncratios['gluon'].GetXaxis().SetTitleSize(0.18)
    pythiauncratios['gluon'].GetXaxis().SetLabelSize(0.18)
    pythiauncratios['gluon'].GetXaxis().SetTitleOffset(0.9)
    #pythiauncratios['gluon'].GetYaxis().SetTitleColor(colors['gluon'])
    pythiauncratios['gluon'].GetYaxis().SetTitle('#splitline{MC/data}{ (gluon)}  ')
    pythiauncratios['gluon'].GetYaxis().SetTitleSize(0.16)
    pythiauncratios['gluon'].GetYaxis().SetLabelSize(0.17)
    pythiauncratios['gluon'].GetYaxis().SetTitleOffset(0.305)
    pythiauncratios['gluon'].Draw('e2')
    line['gluon'].Draw('x0 e1 same')
    line['gluon'].Draw('same')
    pythiaratios['gluon'].Draw('same h')
    pythiafsrupratios['gluon'].Draw('SAME P X0 E1')
    pythiafsrdownratios['gluon'].Draw('SAME P X0 E1')
    #pythiafsrtuneratios['gluon'].Draw('same h')
    herwigratios['gluon'].Draw('same h')
    sherparatios['gluon'].Draw('same h')
    direratios['gluon'].Draw('same h')
    
    c.Print(opt.outDir+'/'+opt.obs+'_'+opt.reco+'_flavors.pdf')
    c.Print(opt.outDir+'/'+opt.obs+'_'+opt.reco+'_flavors.png')

"""
Adapt plots to limited range
(Useful for ratio plots! Otherwise, they are not drawn when the central point is outside the range.)
"""
def limitToRange(h, ratiorange):
    h.GetYaxis().SetRangeUser(ratiorange[0], ratiorange[1])
    for i in xrange(1,h.GetNbinsX()+1):
        up = h.GetBinContent(i) + h.GetBinError(i)
        if (up > ratiorange[1]):
            up = ratiorange[1]
        dn = h.GetBinContent(i) - h.GetBinError(i)
        if (dn < ratiorange[0]):
            dn = ratiorange[0]
        h.SetBinContent(i, (up + dn)/2.)
        h.SetBinError(i, (up - dn)/2.)

def shift(h, step):
   axis = h.GetXaxis()
   width = axis.GetXmax() - axis.GetXmin()
   shift = width / 200 * step
   newbins = [x+shift for x in axis.GetXbins()]
   axis.Set(len(newbins)-1, array('d', newbins))

"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

