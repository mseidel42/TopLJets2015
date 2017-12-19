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

debug = True

"""
steer the script
"""
def main():
    
    cmsLabel='#bf{CMS} #it{preliminary}'
    
    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option(     '--mcUnc',        dest='mcUnc'  ,      help='common MC related uncertainty (e.g. lumi)',        default=0,              type=float)
    parser.add_option(     '--com',          dest='com'  ,        help='center of mass energy',                            default='13 TeV',       type='string')
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',        default='data/era2016/samples.json',              type='string')
    parser.add_option( '--systJson', dest='systJson', help='json with list of systematics', default='data/era2016/syst_samples.json', type='string')
    parser.add_option('-i', '--inDir',       dest='inDir' ,      help='input directory',                default='unfolding/result',              type='string')
    parser.add_option('', '--inDirToys',       dest='inDirToys' ,      help='input toy directory',                default='unfolding/toys',              type='string')
    parser.add_option('-O', '--outDir',      dest='outDir' ,     help='output directory',                default='unfolding',              type='string')
    parser.add_option('-o', '--outName',     dest='outName' ,    help='name of the output file',        default='plotter.root',    type='string')
    parser.add_option(      '--silent',      dest='silent' ,     help='only dump to ROOT file',         default=False,             action='store_true')
    parser.add_option(      '--saveTeX',     dest='saveTeX' ,    help='save as tex file as well',       default=False,             action='store_true')
    parser.add_option('-l', '--lumi',        dest='lumi' ,       help='lumi [/pb]',              default=16551.,              type=float)
    parser.add_option('--obs', dest='obs',  default='mult', help='observable [default: %default]')
    parser.add_option('--flavor', dest='flavor',  default='all', help='flavor [default: %default]')
    parser.add_option('-r', '--reco', dest='reco', default='charged', help='Use charged/puppi/all particles [default: %default]')
    (opt, args) = parser.parse_args()
    
    observables = ["mult", "ptds", "ga_lha", "ga_width", "ga_thrust", "ecc", "zg", "zgdr", "nsd", "tau21", "tau32", "tau43", "c1_00", "c1_02", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "n2_b1", "n3_b1", "m2_b2", "n2_b2", "n3_b2"]
    
    nice_observables_tex = {"mult": "$\\lambda_{0}^{0}$ (N)", "ptds": "$\\lambda_{0}^{2}$ ($p_{T}D^{*})$", "ecc": "$\\varepsilon$", "tau21": "$\\tau_{21}$", "tau32": "$\\tau_{32}$", "tau43": "$\\tau_{43}$", "zg": "$z_{g}$", "zgdr": "$\\Delta R_{g}$", "ga_width": "$\\lambda_{1}^{1}$ (width)", "ga_lha": "$\\lambda_{0.5}^{1}$ (LHA)", "ga_thrust": "$\\lambda_{2}^{1}$ (thrust)", "c1_00": "$C_{1}^{(0.0)}$", "c1_02": "$C_{1}^{(0.2)}$", "c1_05": "$C_{1}^{(0.5)}$", "c1_10": "$C_{1}^{(1.0)}$", "c1_20": "$C_{1}^{(2.0)}$", "c2_00": "$C_{2}^{(0.0)}$", "c2_02": "$C_{2}^{(0.2)}$", "c2_05": "$C_{2}^{(0.5)}$", "c2_10": "$C_{2}^{(1.0)}$", "c2_20":  "$C_{2}^{(2.0)}$", "c3_00": "$C_{3}^{(0.0)}$", "c3_02": "$C_{3}^{(0.2)}$", "c3_05": "$C_{3}^{(0.5)}$", "c3_10": "$C_{3}^{(1.0)}$", "c3_20": "$C_{3}^{(2.0)}$", "m2_b1": "$M_{2}^{(1)}$", "n2_b1": "$N_{2}^{(1)}$", "n3_b1": "$N_{3}^{(1)}$", "m2_b2": "$M_{2}^{(2)}$", "n2_b2": "$N_{2}^{(2)}$", "n3_b2": "$N_{3}^{(2)}$", "nsd": "$n_{SD}$"}
    
    observables_low = ["ga_width", "ecc", "zg", "tau43"]
    
    flavors = ['incl', 'bottom', 'light', 'gluon']

    # Read lists of syst samples
    varList = []
    varExp = [['JEC/JER', ['jec_CorrelationGroupMPFInSitu',
              'jec_RelativeFSR',
              'jec_CorrelationGroupUncorrelated',
              'jec_FlavorPureGluon',
              'jec_FlavorPureQuark',
              'jec_FlavorPureCharm',
              'jec_FlavorPureBottom',
              'jer']],
              ['BTAG', ['btag_heavy',
              'btag_light',
              'csv_heavy',
              'csv_light']],
              ['TRACK', ['tracking']],
              ['BKG', ['singletop',
              'wjets']]
             ]
    for vargroup in varExp:
        grouplist = []
        for var in vargroup[1]:
            grouplist.append([var, [var+'_up', var+'_down']])
        varList.append([vargroup[0], grouplist])
    varModel = [['MT', [['mt', ['m171v5', 'm173v5']]]],
                ['HERWIG', [['herwig', ['herwig']]]],
                ['FSR', [['fsr', ['fsrup', 'fsrdn']]]],
                ['UE/CR', [['ue', ['ueup', 'uedn']],
                 ['erdON', ['erdON']],
                 ['qcdBased', ['qcdBased']],
                 #['gluonMove', ['gluonMove']], # TODO: needs reprocessing with new normCache
                ]],
                ['FRAG', [['evtgen', ['evtgen']],
                 ['bfragbl', ['wgt7', 'wgt8']], # b frag Bowler-Lund up/down
                 ['bfragpeterson', ['wgt9']], # b frag Peterson
                 ['semilepbr', ['wgt10', 'wgt11']], # B hadron semilep BR
                ]],
                ['ME/ISR', [['toppt', ['wgt12']], # top pt reweighting
                 ['muf', ['wgt13', 'wgt14']], # muF
                 ['mur', ['wgt15', 'wgt18']], # muR
                 ['mufmur', ['wgt16', 'wgt20']], # muF+muR
                 ['hdamp', ['hdampup', 'hdampdn']],
                 ['isr', ['isrup', 'isrdn']],
                ]]
               ]
    varList += varModel
    
    varExpWgt = [['PU', [['pu', ['wgt1', 'wgt2']]]], # PU
                 #['lepton', ['leptrigger', ['wgt3', 'wgt4']], # lepton trigger
                 # ['lepsel', ['wgt5', 'wgt6']], # lepton selection
                 #]
                 ['STAT', [['stat']]],
                ]
    varList += varExpWgt
    
    varModelDict = {'evtgen': 'EvtGen',
                    'm171v5': 'mt down',
                    'm173v5': 'mt up',
                    'herwig': 'Herwig++',
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
                    #'gluonMove': 'CR: gluon-move',
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
                    'nominalGen': 'nominal sample'
                    }

    with open('%s/syst_%s.tex'%(opt.outDir, opt.reco), 'w') as tex:
        tex.write('Observable & Flavor')
        for vargroup in varList:
            for var in vargroup[1]:
                tex.write(' & %s'%(var[0]))
        tex.write(' \\\\\n \\hline \n')
        
        for obs in observables:
            for flavor in flavors:                
                resultfile = '%s/%s_%s_%s_result.root'%(opt.inDir, obs, opt.reco, flavor)
                fIn=ROOT.TFile.Open(resultfile)
                
                # reference
                hdata = fIn.Get('MC13TeV_TTJets_Unfolded')
                
                tex.write('%s & %s'%(nice_observables_tex[obs], flavor))
                
                for vargroup in varList:
                    for var in vargroup[1]:
                        syst_min = 999.
                        syst_max = 0.
                        if var[0] == 'stat':
                            for i in range(1, hdata.GetNbinsX()+1):
                                stat = hdata.GetBinError(i)/hdata.GetBinContent(i)*100.
                                if stat < syst_min: syst_min = stat
                                if stat > syst_max: syst_max = stat
                        else:
                            for vardir in var[1]:
                                hsyst = fIn.Get('MC13TeV_TTJets_'+vardir+'_Unfolded')
                                if not hsyst:
                                    print(vardir, 'not loaded')  
                                    continue
                                
                                shifts = []
                                for i in range(1, hdata.GetNbinsX()+1):
                                    shifts.append(abs(hsyst.GetBinContent(i) - hdata.GetBinContent(i))/hdata.GetBinContent(i)*100.)
                                
                                if min(shifts) < syst_min: syst_min = min(shifts)
                                if max(shifts) > syst_max: syst_max = max(shifts)
                                
                        tex.write(' & %.1f'%(syst_max))
                tex.write(' \\\\\n')

    with open('%s/syst_groups_%s.tex'%(opt.outDir, opt.reco), 'w') as tex:
        tex.write('Observable & Flavor')
        for vargroup in varList:
            tex.write(' & %s'%(vargroup[0]))
        tex.write(' \\\\\n \\hline \n \\hline \n')
        
        for obs in observables:
            tex.write('%s\n'%(nice_observables_tex[obs]))
            for flavor in flavors:
                resultfile = '%s/%s_%s_%s_result.root'%(opt.inDir, obs, opt.reco, flavor)
                fIn=ROOT.TFile.Open(resultfile)
                
                # reference
                hdata = fIn.Get('MC13TeV_TTJets_Unfolded')
                
                tex.write(' & %s'%(flavor))
                
                for vargroup in varList:
                    groupbins = hdata.GetNbinsX()*[0.]
                    for var in vargroup[1]:
                        if var[0] == 'stat':
                            for i in range(1, hdata.GetNbinsX()+1):
                                groupbins[i-1] = hdata.GetBinError(i)/hdata.GetBinContent(i)*100.
                            continue
                        varbins = hdata.GetNbinsX()*[0.]
                        for vardir in var[1]:
                            hsyst = fIn.Get('MC13TeV_TTJets_'+vardir+'_Unfolded')
                            if not hsyst:
                                print(vardir, 'not loaded')  
                                continue
                            
                            for i in range(1, hdata.GetNbinsX()+1):
                                shift = abs(hsyst.GetBinContent(i) - hdata.GetBinContent(i))/hdata.GetBinContent(i)*100.
                                if shift > varbins[i-1]:
                                    varbins[i-1] = shift
                            for i in range(len(groupbins)):
                                groupbins[i] = sqrt(groupbins[i]**2 + varbins[i]**2)
                            
                    tex.write(' & %.1f -- %.1f'%(min(groupbins), max(groupbins)))
                tex.write(' \\\\\n')
            tex.write('\\hline \n')

    # FIXME: this part does not work yet
    with open('%s/syst_groups_low_incl.tex'%(opt.outDir), 'w') as tex:
        for obs in observables:
            if not obs in observables_low: continue
            tex.write(' & %s'%(nice_observables_tex[obs]))
        tex.write(' \\\\\n \\hline \n')
        
        for vargroup in varList:
            tex.write('%s'%(vargroup[0]))
            for obs in observables:
                if not obs in observables_low: continue
                for flavor in flavors:
                    if not flavor == 'all': continue
                    resultfile = '%s/%s_charged_%s_result.root'%(opt.inDir, obs, flavor)
                    fIn=ROOT.TFile.Open(resultfile)
                    
                    # reference
                    hdata = fIn.Get('MC13TeV_TTJets_Unfolded')
                    
                    
                    groupbins = hdata.GetNbinsX()*[0.]
                    for var in vargroup[1]:
                        if var[0] == 'stat':
                            for i in range(1, hdata.GetNbinsX()+1):
                                groupbins[i-1] = hdata.GetBinError(i)/hdata.GetBinContent(i)*100.
                            continue
                        varbins = hdata.GetNbinsX()*[0.]
                        for vardir in var[1]:
                            hsyst = fIn.Get('MC13TeV_TTJets_'+vardir+'_Unfolded')
                            if not hsyst:
                                print(vardir, 'not loaded')  
                                continue
                            
                            for i in range(1, hdata.GetNbinsX()+1):
                                shift = abs(hsyst.GetBinContent(i) - hdata.GetBinContent(i))/hdata.GetBinContent(i)*100.
                                if shift > varbins[i-1]:
                                    varbins[i-1] = shift
                            for i in range(len(groupbins)):
                                groupbins[i] = sqrt(groupbins[i]**2 + varbins[i]**2)
                            
                    tex.write(' & %.1f -- %.1f'%(min(groupbins), max(groupbins)))
            tex.write(' \\\\\n')
        
"""
for execution from another script
"""
if __name__ == "__main__":
    main()
    #sys.exit(main())

