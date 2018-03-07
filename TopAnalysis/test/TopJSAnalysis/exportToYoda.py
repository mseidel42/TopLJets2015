import os
import pickle
import ROOT
import sys

recos = ['charged', 'all']

observables = ["mult", "ptds", "ga_lha", "ga_width", "ga_thrust", "ecc", "zg", "zgdr", "nsd", "tau21", "tau32", "tau43", "c1_00", "c1_02", "c1_05", "c1_10", "c1_20", "c2_00", "c2_02", "c2_05", "c2_10", "c2_20", "c3_00", "c3_02", "c3_05", "c3_10", "c3_20", "m2_b1", "n2_b1", "n3_b1", "m2_b2", "n2_b2", "n3_b2"]
    
nice_observables_tex = {"mult": "$\\lambda_{0}^{0}$ (N)", "ptds": "$\\lambda_{0}^{2}$ ($p_{T}^{d,*})$", "ecc": "$\\varepsilon$", "tau21": "$\\tau_{21}$", "tau32": "$\\tau_{32}$", "tau43": "$\\tau_{43}$", "zg": "$z_{g}$", "zgdr": "$\\Delta R_{g}$", "ga_width": "$\\lambda_{1}^{1}$ (width)", "ga_lha": "$\\lambda_{0.5}^{1}$ (LHA)", "ga_thrust": "$\\lambda_{2}^{1}$ (thrust)", "c1_00": "$C_{1}^{(0.0)}$", "c1_02": "$C_{1}^{(0.2)}$", "c1_05": "$C_{1}^{(0.5)}$", "c1_10": "$C_{1}^{(1.0)}$", "c1_20": "$C_{1}^{(2.0)}$", "c2_00": "$C_{2}^{(0.0)}$", "c2_02": "$C_{2}^{(0.2)}$", "c2_05": "$C_{2}^{(0.5)}$", "c2_10": "$C_{2}^{(1.0)}$", "c2_20":  "$C_{2}^{(2.0)}$", "c3_00": "$C_{3}^{(0.0)}$", "c3_02": "$C_{3}^{(0.2)}$", "c3_05": "$C_{3}^{(0.5)}$", "c3_10": "$C_{3}^{(1.0)}$", "c3_20": "$C_{3}^{(2.0)}$", "m2_b1": "$M_{2}^{(1)}$", "n2_b1": "$N_{2}^{(1)}$", "n3_b1": "$N_{3}^{(1)}$", "m2_b2": "$M_{2}^{(2)}$", "n2_b2": "$N_{2}^{(2)}$", "n3_b2": "$N_{3}^{(2)}$", "nsd": "$n_{SD}$"}

flavors = ['incl', 'bottom', 'light', 'gluon']

baseDir='./'
outDir='./'
if len(sys.argv)>1:  baseDir=sys.argv[1]
if len(sys.argv)>2: outDir=sys.argv[2]

name = 'CMS_2017_PAS_TOP_17_013'

fOut={'data' : open(os.path.join(outDir,name+'.yoda'),'w'),
      'mc'   : open(os.path.join(outDir,name+'_MC.yoda'),'w'),
      'plot' : open(os.path.join(outDir,name+'.plot'),'w')}


fOut['plot'].write('''# BEGIN PLOT /CMS_2017_PAS_TOP_17_013/*
LogY=0
LogX=0
XTwosidedTicks=1
YTwosidedTicks=1
LegendXPos=0.55
# END PLOT
''')

ireco = 0
for reco in recos:
    ireco += 1
    iobs = 0
    for obs in observables:
        iobs += 1
        iflavor = 0
        for flavor in flavors:
            iflavor += 1
            pname = 'd{0:0>2}-x{1:0>2}-y{2:0>2}'.format(ireco,iobs,iflavor)
            
            #readout data
            gr = {}
            resultfile = 'unfolding/result/%s_%s_%s_result.root'%(obs, reco, flavor)
            fIn=ROOT.TFile.Open(resultfile)
            gr['data'] = fIn.Get('MC13TeV_TTJets_Unfolded')
            gr['syst'] = fIn.Get('dataUnfoldedSys')
            gr['mc']   = fIn.Get('nominalGen')
            
            fOut['plot'].write('''
# BEGIN PLOT /CMS_2017_PAS_TOP_17_013/%s
Title=CMS, 13 TeV, $t\\bar{t}$ lepton+jets, %s jets, %s particles
XLabel=%s
YLabel=1/N dN/d%s
'''%(pname, flavor, reco, nice_observables_tex[obs], nice_observables_tex[obs]))
            if (gr['mc'].GetXaxis().GetBinCenter(gr['mc'].GetMaximumBin()) > (gr['mc'].GetXaxis().GetXmax() - gr['mc'].GetXaxis().GetXmin())*0.5):
                fOut['plot'].write('LegendXPos=0.05\n')
            fOut['plot'].write('# END PLOT\n')
            
            #dump plot in yoda format
            for t in ['data','mc']:

                fOut[t].write('BEGIN YODA_SCATTER2D /REF/%s/%s\n'%(name,pname))
                if t == 'data': fOut[t].write('Path=/REF/%s/%s\n'%(name,pname))
                else: fOut[t].write('Path=/%s/%s\n'%(name,pname))
                fOut[t].write('Type=Scatter2D\n')
                fOut[t].write('# HISTID: %s %s %s\n'%(reco,obs,flavor))
                fOut[t].write('# xval xerr- xerr+ yval yerr- yerr+\n')
                eylo,eyhi = 0., 0.
                for i in range(1,gr[t].GetNbinsX()+1):
                    xref = gr[t].GetBinCenter(i)
                    exlo = gr[t].GetBinWidth(i)/2.
                    exhi = gr[t].GetBinWidth(i)/2.
                    
                    y = gr[t].GetBinContent(i)
                    if t == 'data':
                        eylo = y - (gr['syst'].GetBinContent(i) - gr['syst'].GetBinError(i))
                        eyhi = gr['syst'].GetBinContent(i) + gr['syst'].GetBinError(i) - y
                    
                    fOut[t].write('%6.6g %6.6g %6.6g %.6g %.6g %.6g\n'%(float(xref),exlo,exhi,float(y),eylo,eyhi))
                fOut[t].write('END YODA_SCATTER2D\n')
                fOut[t].write('\n')

#close open files
for t in ['data','mc']: 
    fOut[t].write('\n\n# this file was generated from %s\n'%t)
    fOut[t].close()
