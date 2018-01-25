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
    parser.add_option('--var', dest='var', default='abs(gj_partonflavor)', help='observable [default: %default]')
    parser.add_option('-n', '--nevents', dest='nevents', help='number of events', default=10000, type=int)
    parser.add_option('--sample', dest='sample', default='', help='sample appendix')
    (opt, args) = parser.parse_args()
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    #read lists of samples
    flavors = ['bottom', 'light', 'gluon']
    flavormap = {'bottom': 5, 'light': 1, 'gluon': 0}
    partonflavors = {'bottom': [5], 'light': range(1,5), 'gluon': [0,21]}
    hist = {}
    fraction = {}
    
    eosdir = '/eos/user/m/mseidel/analysis/TopJetShapes/b312177_new/Chunks/'

    t_mc = ROOT.TChain('tjsev')
    if opt.sample == '':
        for i in range(137):
            t_mc.Add(eosdir+'MC13TeV_TTJets_'+str(i)+'.root')
    else:
        t_mc.Add(eosdir+'MC13TeV_TTJets_'+opt.sample+'_*.root')
    print('t_mc.GetEntries()', t_mc.GetEntries())
    
    for flavor in flavors:
        hist[flavor] = ROOT.TH1F(flavor, '', 21, 1, 22)
        t_mc.Draw(opt.var+' >> '+flavor, '(gen_sel == 1 & abs(gj_eta) < 2. & gj_overlap == 0 & gj_flavor == '+str(flavormap[flavor])+')*weight[0]', '', opt.nevents)
        
        hist[flavor].Scale(1./hist[flavor].Integral())
        
        for partonflavor in flavors:
            fraction[partonflavor] = 0.
            for xbin in partonflavors[partonflavor]:
                fraction[partonflavor] += hist[flavor].GetBinContent(xbin)
        
        print('= = = = = = = = = =')
        print('Particle level: %s'%flavor)
        print('= = = = = = = = = =')
        
        print('Parton level')
        print('bottom: %2.2f %%'%(100.*fraction['bottom']))
        print('light:  %2.2f %%'%(100.*fraction['light']))
        print('gluon:  %2.2f %%'%(100.*fraction['gluon']))


if __name__ == "__main__":
	sys.exit(main())
