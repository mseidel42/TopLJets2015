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

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS} #it{Simulation}'
    
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
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='incl', help='flavor [default: %default]')
    parser.add_option('-r', '--reco', dest='reco', default='charged', help='Use charged/puppi/all particles [default: %default]')
    (opt, args) = parser.parse_args()
    
    observables = ["mult", "ptds", "ga_lha", "ga_width", "ga_thrust", "ecc", "zg", "zgdr", "nsd", "tau21", "tau32", "tau43", "c1_00", "c1_02", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "n2_b1", "n3_b1", "m2_b2", "n2_b2", "n3_b2"]
    
    nice_observables_tex = {"mult": "$\\lambda_{0}^{0}$ (N)", "ptds": "$\\lambda_{0}^{2}$ ($p_{T}D^{*})$", "ecc": "$\\varepsilon$", "tau21": "$\\tau_{21}$", "tau32": "$\\tau_{32}$", "tau43": "$\\tau_{43}$", "zg": "$z_{g}$", "zgdr": "$\\Delta R_{g}$", "ga_width": "$\\lambda_{1}^{1}$ (width)", "ga_lha": "$\\lambda_{0.5}^{1}$ (LHA)", "ga_thrust": "$\\lambda_{2}^{1}$ (thrust)", "c1_00": "$C_{1}^{(0.0)}$", "c1_02": "$C_{1}^{(0.2)}$", "c1_05": "$C_{1}^{(0.5)}$", "c1_10": "$C_{1}^{(1.0)}$", "c1_20": "$C_{1}^{(2.0)}$", "c2_00": "$C_{2}^{(0.0)}$", "c2_02": "$C_{2}^{(0.2)}$", "c2_05": "$C_{2}^{(0.5)}$", "c2_10": "$C_{2}^{(1.0)}$", "c2_20":  "$C_{2}^{(2.0)}$", "c3_00": "$C_{3}^{(0.0)}$", "c3_02": "$C_{3}^{(0.2)}$", "c3_05": "$C_{3}^{(0.5)}$", "c3_10": "$C_{3}^{(1.0)}$", "c3_20": "$C_{3}^{(2.0)}$", "m2_b1": "$M_{2}^{(1)}$", "n2_b1": "$N_{2}^{(1)}$", "n3_b1": "$N_{3}^{(1)}$", "m2_b2": "$M_{2}^{(2)}$", "n2_b2": "$N_{2}^{(2)}$", "n3_b2": "$N_{3}^{(2)}$", "nsd": "$n_{SD}$"}
    
    observables_low = ["ptds", "ecc", "tau43", "zg", "zgdr"]
    observables_low = ["ga_width", "ecc", "zg", "tau43"]
    
    flavors = ['incl', 'bottom', 'light', 'gluon']

    # Read lists of syst samples
    varList = []
    varExp = ['jec_CorrelationGroupMPFInSitu',
              'jec_RelativeFSR',
              'jec_CorrelationGroupUncorrelated',
              'jec_FlavorPureGluon',
              'jec_FlavorPureQuark',
              'jec_FlavorPureCharm',
              'jec_FlavorPureBottom',
              'jer',
              'btag_heavy',
              'btag_light',
              'csv_heavy',
              'csv_light',
              'tracking',
              'singletop',
              'wjets'
             ]
    for var in varExp:
        varList.append([var+'_up', var+'_down'])
    varModel = [['evtgen'],
                ['m171v5', 'm173v5'],
                ['herwig'],
                ['isrup', 'isrdn'],
                ['fsrup', 'fsrdn'],
                ['hdampup', 'hdampdn'],
                ['ueup', 'uedn'],
                ['erdON'],
                ['qcdBased'],
                ['gluonMove'],
                ['wgt7', 'wgt8'], # b frag Bowler-Lund up/down
                ['wgt9'], # b frag Peterson
                ['wgt10', 'wgt11'], # B hadron semilep BR
                ['wgt12'], # top pt reweighting
                ['wgt13', 'wgt14'], # muF
                ['wgt15', 'wgt18'], # muR
                ['wgt16', 'wgt20'], # muF+muR
               ]
    varList += varModel
    
    varExpWgt = [['wgt1', 'wgt2'], # PU
                 ['wgt3', 'wgt4'], # lepton trigger
                 ['wgt5', 'wgt6'], # lepton selection
                ]
    varList += varExpWgt

    modelsToTest = varModel + [['cflip'], ['nominalGen']]
    modelsToTest.append(['pythia8_asfsr0.1365_meoff_crdefault'])
    modelsToTest.append(['herwig7'])
    modelsToTest.append(['sherpa'])
    #FSR scan
    modelsToTest.append(['herwigpp_asfsr0.100_meon_crdefault'])
    modelsToTest.append(['herwigpp_asfsr0.110_meon_crdefault'])
    modelsToTest.append(['herwigpp_asfsr0.115_meon_crdefault'])
    modelsToTest.append(['herwigpp_asfsr0.120_meon_crdefault'])
    modelsToTest.append(['herwigpp_asfsr0.125_meon_crdefault'])
    modelsToTest.append(['herwigpp_asfsr0.130_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.070_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.080_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.090_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.100_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.105_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.110_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.115_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.120_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.125_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.130_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.135_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.140_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.150_meon_crdefault'])
    modelsToTest.append(['pythia8_asfsr0.160_meon_crdefault'])
    
    modelsToTex = ['fsrdn', 'nominalGen', 'fsrup', 'herwig7', 'sherpa']
    
    obsgroups = ['all', 'low']
    
    unsummedChi2 = OrderedDict()
    for obs in observables + obsgroups:
        unsummedChi2[obs] = {}
        for var in modelsToTest:
            for vardir in var:
                unsummedChi2[obs][vardir] = {}
                for flavor in flavors:
                    unsummedChi2[obs][vardir][flavor] = 0.
    
    varModelDict = {'evtgen': 'EvtGen',
                    'm171v5': 'mt down',
                    'm173v5': 'mt up',
                    'herwig': 'Herwig++',
                    'herwig7': 'Herwig7',
                    'sherpa': 'Sherpa',
                    'isrup': 'ISR up',
                    'isrdn': 'ISR down',
                    'fsrup': 'FSR up',
                    'fsrdn': 'FSR down',
                    'hdampup': 'hdamp up',
                    'hdampdn': 'hdamp down',
                    'ueup': 'UE up',
                    'uedn': 'UE down',
                    'erdON': 'CR: erd on',
                    'qcdBased': 'CR: QCD-inspired',
                    'gluonMove': 'CR: gluon-move',
                    'wgt7': 'b frag up',
                    'wgt8': 'b frag down', # b frag Bowler-Lund up/down
                    'wgt9': 'b frag Peterson', # b frag Peterson
                    'wgt10': 'B semilep BR up',
                    'wgt11': 'B semilep BR down', # B hadron semilep BR
                    'wgt12': 'top pt', # top pt reweighting
                    'wgt13': 'muF up',
                    'wgt14': 'muF down', # muF
                    'wgt15': 'muR up',
                    'wgt18': 'muR down', # muR
                    'wgt16': 'muF+muR up',
                    'wgt20': 'muF+muR down', # muF+muR
                    'cflip': 'Color octet W',
                    'nominalGen': 'nominal sample',
                    'herwigpp_asfsr0.100_meon_crdefault' : 'H++ asfsr=0.100',
                    'herwigpp_asfsr0.110_meon_crdefault' : 'H++ asfsr=0.110',
                    'herwigpp_asfsr0.115_meon_crdefault' : 'H++ asfsr=0.115',
                    'herwigpp_asfsr0.120_meon_crdefault' : 'H++ asfsr=0.120',
                    'herwigpp_asfsr0.125_meon_crdefault' : 'H++ asfsr=0.125',
                    'herwigpp_asfsr0.130_meon_crdefault' : 'H++ asfsr=0.130',
                    'pythia8_asfsr0.070_meon_crdefault' : 'P8 asfsr=0.070',
                    'pythia8_asfsr0.080_meon_crdefault' : 'P8 asfsr=0.080',
                    'pythia8_asfsr0.090_meon_crdefault' : 'P8 asfsr=0.090',
                    'pythia8_asfsr0.100_meon_crdefault' : 'P8 asfsr=0.100',
                    'pythia8_asfsr0.105_meon_crdefault' : 'P8 asfsr=0.105',
                    'pythia8_asfsr0.110_meon_crdefault' : 'P8 asfsr=0.110',
                    'pythia8_asfsr0.115_meon_crdefault' : 'P8 asfsr=0.115',
                    'pythia8_asfsr0.120_meon_crdefault' : 'P8 asfsr=0.120',
                    'pythia8_asfsr0.125_meon_crdefault' : 'P8 asfsr=0.125',
                    'pythia8_asfsr0.130_meon_crdefault' : 'P8 asfsr=0.130',
                    'pythia8_asfsr0.135_meon_crdefault' : 'P8 asfsr=0.135',
                    'pythia8_asfsr0.140_meon_crdefault' : 'P8 asfsr=0.140',
                    'pythia8_asfsr0.150_meon_crdefault' : 'P8 asfsr=0.150',
                    'pythia8_asfsr0.160_meon_crdefault' : 'P8 asfsr=0.160',
                    'pythia8_asfsr0.1365_meoff_crdefault' : 'ME corr. off',
                    }
    
    sumNominal = 0.
    sumFSRUp = 0.
    sumFSRDown = 0.
    sumHerwig = 0.
    sumSherpa = 0.
    
    sumLowNominal = {}
    sumLowFSRUp   = {}
    sumLowFSRDown = {}
    sumLowHerwig  = {}
    sumLowSherpa  = {}
    sumLowNominal['total'] = 0.
    sumLowFSRUp  ['total'] = 0.
    sumLowFSRDown['total'] = 0.
    sumLowHerwig ['total'] = 0.
    sumLowSherpa ['total'] = 0.
    for flavor in flavors:
        sumLowNominal[flavor] = 0.
        sumLowFSRUp  [flavor] = 0.
        sumLowFSRDown[flavor] = 0.
        sumLowHerwig [flavor] = 0.
        sumLowSherpa [flavor] = 0.
    
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

    with open('unfolding/chi2table_%s.tex'%(opt.reco), 'w') as tex:
        for obs in observables:
            for flavor in flavors:
        
                # statistical covariance
                
                toyfile = '%s/%s_%s_%s_toys.root'%(opt.inDirToys, obs, opt.reco, flavor)
                fInToy = ROOT.TFile.Open(toyfile)
                if not fInToy: continue
                
                pseudoresults = []
                counter = 0
                while True:
                    pseudoresult = []
                    h = fInToy.Get('Unfolded_' + str(counter))
                    if not h: break
                    for i in range(1, h.GetNbinsX()+1):
                        pseudoresult.append(h.GetBinContent(i)/h.GetBinWidth(i))
                    integral = sum(pseudoresult)
                    for i in range(len(pseudoresult)):
                        pseudoresult[i] = pseudoresult[i]/integral
                    pseudoresults.append(pseudoresult)
                    counter += 1
                
                print(obs, flavor, 'imported', counter, 'toy experiments')
                
                x = numpy.array(pseudoresults).T
                #print(x)
                statcov = numpy.cov(x)
                #print(statcov)
                #statcov_reduced = numpy.delete(statcov, 0, 0)
                #statcov_reduced = numpy.delete(statcov_reduced, 0, 1)
                #print(statcov_reduced)
                #print(numpy.linalg.det(statcov_reduced))
                #print(numpy.linalg.inv(statcov_reduced))
                
                #rootoutfile = ROOT.TFile.Open(opt.outDir+'/'+obs+'_charged_'+flavor+'_cov.root', 'RECREATE')
                #rootoutfile.cd()
                
                c = ROOT.TCanvas('c','c',500,500)
                c.SetRightMargin(0.15)
                c.SetLeftMargin(0.12)
                c.SetTopMargin(0.1)
                c.SetTopMargin(0.05)
                c.cd()
                
                h = fInToy.Get('Unfolded_0')
                dataStatCovNorm = ROOT.TH2D('dataStatCovNorm', '', h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray(), h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray())
                axistitle = h.GetXaxis().GetTitle().replace('generated ', '')
                dataStatCovNorm.SetXTitle(axistitle + ' (%s)'%(flavor))
                dataStatCovNorm.SetYTitle(axistitle + ' (%s)'%(flavor))
                
                for i in range(1, dataStatCovNorm.GetNbinsX()+1):
                    for j in range(1, dataStatCovNorm.GetNbinsY()+1):
                        dataStatCovNorm.SetBinContent(i, j, statcov[i-1][j-1])
                
                dataStatCovNorm.GetZaxis().SetRangeUser(-max(abs(statcov.min()), abs(statcov.max())),
                                                         max(abs(statcov.min()), abs(statcov.max())))
                dataStatCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Statistical covariance')
                
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov_stat.pdf')
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov_stat.png')
                
                resultfile = '%s/%s_%s_%s_result.root'%(opt.inDir, obs, opt.reco, flavor)
                fIn=ROOT.TFile.Open(resultfile)
                
                # reference
                hdata = fIn.Get('MC13TeV_TTJets_Unfolded')
                data  = []
                for i in range(1, hdata.GetNbinsX()+1):
                    data.append(hdata.GetBinContent(i))
                
                systcov = copy.copy(statcov)
                systcov -= systcov
                
                for var in varList:
                    if len(var) == 1:
                        hsyst = fIn.Get('MC13TeV_TTJets_'+var[0]+'_Unfolded')
                        if not hsyst:
                            print(var, 'not loaded')  
                            continue
                        
                        syst = []
                        for i in range(1, hsyst.GetNbinsX()+1):
                            syst.append(hsyst.GetBinContent(i))
                        
                        x = numpy.array([data, syst]).T
                        cov = numpy.cov(x)
                        systcov += cov
                    if len(var) == 2:
                        hsyst_up = fIn.Get('MC13TeV_TTJets_'+var[0]+'_Unfolded')
                        hsyst_dn = fIn.Get('MC13TeV_TTJets_'+var[1]+'_Unfolded')
                        if not hsyst_up or not hsyst_dn:
                            print(var, 'not loaded')  
                            continue
                        
                        up = []
                        dn = []
                        maxdelta = []
                        for i in range(1, hsyst_up.GetNbinsX()+1):
                            up.append(hsyst_up.GetBinContent(i))
                            dn.append(hsyst_dn.GetBinContent(i))
                            maxdelta.append(max(abs(hsyst_up.GetBinContent(i) - data[i-1]), abs(hsyst_dn.GetBinContent(i) - data[i-1])))
                        
                        x = numpy.array([up, dn]).T
                        cov = numpy.cov(x)
                        
                        for i in range(len(maxdelta)):
                            for j in range(len(maxdelta)):
                                cov[i][j] = maxdelta[i] * maxdelta[j] * numpy.sign(cov[i][j])
                        
                        systcov += cov
                
                dataSystCovNorm = dataStatCovNorm.Clone('dataSystCovNorm')
                dataSystCovNorm.Reset()
                
                for i in range(1, dataSystCovNorm.GetNbinsX()+1):
                    for j in range(1, dataSystCovNorm.GetNbinsY()+1):
                        dataSystCovNorm.SetBinContent(i, j, systcov[i-1][j-1])
                
                dataSystCovNorm.GetZaxis().SetRangeUser(-max(abs(systcov.min()), abs(systcov.max())),
                                                         max(abs(systcov.min()), abs(systcov.max())))
                dataSystCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Systematic covariance')
                
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov_syst.pdf')
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov_syst.png')
                
                dataCovNorm = dataStatCovNorm.Clone('dataCovNorm')
                dataCovNorm.Reset()
                
                cov = statcov + systcov
                
                for i in range(1, dataCovNorm.GetNbinsX()+1):
                    for j in range(1, dataCovNorm.GetNbinsY()+1):
                        dataCovNorm.SetBinContent(i, j, cov[i-1][j-1])
                
                dataCovNorm.GetZaxis().SetRangeUser(-max(abs(cov.min()), abs(cov.max())),
                                                     max(abs(cov.min()), abs(cov.max())))
                dataCovNorm.Draw('colz')
                
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
                txt.DrawLatex(0.16,0.91, cmsLabel)
                txt.DrawLatex(0.63,0.97, '#scale[0.8]{%3.1f fb^{-1} (%s)}' % (opt.lumi/1000.,opt.com) )
                txt.DrawLatex(0.16,0.85, 'Total covariance')
                
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov.pdf')
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov.png')
                c.Print(opt.outDir+'/'+obs+'_'+opt.reco+'_'+flavor+'_cov.root')
                
                #print(numpy.linalg.det(cov))
                cov_reduced = numpy.delete(cov, 0, 0)
                cov_reduced = numpy.delete(cov_reduced, 0, 1)
                #print(cov_reduced)
                #print(numpy.linalg.det(cov_reduced))
                #print(numpy.linalg.inv(cov_reduced))
                
                for var in modelsToTest:
                    for vardir in var:
                        if vardir in ['nominalGen']:
                            prediction = vardir
                        else:
                            prediction = 'MC13TeV_TTJets_'+vardir+'_gen'
                        unsummedChi2[obs][vardir][flavor] = returnChi2(fIn, cov_reduced, data, prediction)
                        unsummedChi2['all'][vardir][flavor] += unsummedChi2[obs][vardir][flavor]
                        if obs in observables_low:
                            unsummedChi2['low'][vardir][flavor] += unsummedChi2[obs][vardir][flavor]
                
        for obs in observables:
            tex.write('\\hline\n%s\n'%(nice_observables_tex[obs]))
            for flavor in flavors:
                tex.write(' & %s'%(flavor))
                for model in modelsToTex:
                    tex.write(' & %.1f'%(unsummedChi2[obs][model][flavor]))
                tex.write('\\\\\n')
        
        nice_groups  = {'all': 'All observables', 'low': 'Low correlation'}
        for group in obsgroups:
            tex.write('\\hline\n\\hline\n%s'%(nice_groups[group]))
            for flavor in flavors:
                tex.write(' & %s'%(flavor))
                for model in modelsToTex:
                    tex.write(' & %.1f'%(unsummedChi2[group][model][flavor]))
                tex.write('\\\\\n')
                
        #tex.write('\\hline\nTotal &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(sumNominal, sumFSRUp, sumFSRDown, sumHerwig))
        #tex.write('\\hline\nTotal low corr. &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(sumLowNominal['total'], sumLowFSRUp['total'], sumLowFSRDown['total'], sumLowHerwig['total']))
        #for flavor in flavors:
        #    tex.write('-- %s &  & %.1f & %.1f & %.1f & %.1f \\\\\n'%(flavor, sumLowNominal[flavor], sumLowFSRUp[flavor], sumLowFSRDown[flavor], sumLowHerwig[flavor]))
        #    probLowNominal = ROOT.TMath.Prob(sumLowNominal[flavor], len(observables_low))
        #    probLowFSRUp   = ROOT.TMath.Prob(sumLowFSRUp[flavor], len(observables_low))
        #    probLowFSRDown = ROOT.TMath.Prob(sumLowFSRDown[flavor], len(observables_low))
        #    probLowHerwig  = ROOT.TMath.Prob(sumLowHerwig[flavor], len(observables_low))
        #    tex.write('$P(\\chi^{2})$ &  & %.3f & %.3f & %.3f & %.3f \\\\\n'%(probLowNominal, probLowFSRUp, probLowFSRDown, probLowHerwig))
        #    
        
        tex.write('\n'*4)
        
        for var in modelsToTest:
            for vardir in var:
                tex.write('%s & $\chi^{2}$/ndf'%(varModelDict[vardir]))
                for flavor in flavors:
                    tex.write(' & %.1f'%(unsummedChi2['low'][vardir][flavor]))
                tex.write('\\\\\n')
        
        pickle.dump(unsummedChi2, open("unsummedChi2_"+opt.reco+".pkl", "wb"))
    
def returnChi2(fIn, cov_reduced, data, prediction):
    hpred = fIn.Get(prediction)
    pred  = []
    for i in range(1, hpred.GetNbinsX()+1):
        pred.append(hpred.GetBinContent(i))
    diff = []
    for i in range(len(data)):
        diff.append(pred[i] - data[i])
    
    chi2 = numpy.array(diff[1:]).T.dot(numpy.linalg.inv(cov_reduced).dot(numpy.array(diff[1:])))
    ndf  = hpred.GetNbinsX()-1
    #prob = ROOT.TMath.Prob(chi2, ndf)
    return chi2/ndf

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

