import ROOT
import math
import os,sys
from collections import OrderedDict

from TopLJets2015.TopAnalysis.rounding import *

"""
increments the first and the last bin to show the under- and over-flows
"""
def fixExtremities(h,addOverflow=True,addUnderflow=True):
    if addUnderflow :
        fbin  = h.GetBinContent(0) + h.GetBinContent(1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(0)*h.GetBinError(0) + h.GetBinError(1)*h.GetBinError(1))
	h.SetBinContent(1,fbin)
	h.SetBinError(1,fbine)
	h.SetBinContent(0,0)
	h.SetBinError(0,0)
    if addOverflow:
        nbins = h.GetNbinsX();
	fbin  = h.GetBinContent(nbins) + h.GetBinContent(nbins+1)
	fbine = ROOT.TMath.Sqrt(h.GetBinError(nbins)*h.GetBinError(nbins)  + h.GetBinError(nbins+1)*h.GetBinError(nbins+1))
	h.SetBinContent(nbins,fbin)
	h.SetBinError(nbins,fbine)
	h.SetBinContent(nbins+1,0)
	h.SetBinError(nbins+1,0)

"""
Adapt plots to limited range
(Useful for ratio plots! Otherwise, they are not drawn when the central point is outside the range.)
"""
def limitToRange(h, ratiorange):
    for i in xrange(1,h.GetNbinsX()+1):
        up = h.GetBinContent(i) + h.GetBinError(i)
        if (up > ratiorange[1]):
            up = ratiorange[1]
        dn = h.GetBinContent(i) - h.GetBinError(i)
        if (dn < ratiorange[0]):
            dn = ratiorange[0]
        h.SetBinContent(i, (up + dn)/2.)
        h.SetBinError(i, (up - dn)/2.)
        

"""
A wrapper to store data and MC histograms for comparison
"""
class Plot(object):

    def __init__(self,name,com='13 TeV'):
        self.name = name
        self.cmsLabel='#bf{CMS}'
        self.com=com
        self.wideCanvas = True if 'ratevsrun' in self.name else False
        self.doPoissonErrorBars=True
        self.mc = OrderedDict()
        self.mcsyst = {}
        self.totalMCUnc = None
        self.spimpose={}
        self.dataH = None
        self.data = None
        self._garbageList = []
        self.plotformats = ['pdf','png']
        self.savelog = False
        self.doChi2 = False
        self.ratiorange = (0.4,1.6)
        self.frameMin=0.01
        self.frameMax=2.
        self.mcUnc=0

    def add(self, h, title, color, isData, spImpose, isSyst):

        if 'ratevsrun' in self.name and not isData: return

        h.SetTitle(title)
        
        #x-title replacements
        if h.GetXaxis().GetTitle() == 'N (charged)':
            h.GetXaxis().SetTitle('#lambda_{0}^{0} (N) (charged)')

        #check if color is given in hexadec format
        try:
            if '#' in color : color=ROOT.TColor.GetColor(color)
        except:
            pass

        if isData:
            try:
                self.dataH.Add(h)
            except:
                self.dataH=h
                self.dataH.SetDirectory(0)
                self.dataH.SetMarkerStyle(20)
                self.dataH.SetMarkerSize(1.)
                self.dataH.SetMarkerColor(color)
                self.dataH.SetLineColor(ROOT.kBlack)
                self.dataH.SetLineWidth(2)
                self.dataH.SetFillColor(ROOT.kBlack)
                self.dataH.SetFillStyle(3001)
                self._garbageList.append(h)
        elif isSyst:
            try:
                self.mcsyst[title].Add(h)
            except:
                self.mcsyst[title]=h
                self.mcsyst[title].SetName('%s_%s' % (h.GetName(), title ) )
                self.mcsyst[title].SetDirectory(0)
                self.mcsyst[title].SetMarkerStyle(1)
                self.mcsyst[title].SetMarkerColor(color)
                self.mcsyst[title].SetLineColor(ROOT.kBlack)
                self.mcsyst[title].SetLineWidth(1)
                self.mcsyst[title].SetFillColor(color)
                self.mcsyst[title].SetFillStyle(1001)
                self._garbageList.append(h)
        else:
            try:
                if spImpose : self.spimpose[title].Add(h)
                else        : self.mc[title].Add(h)
            except:
                h.SetName('%s_%s' % (h.GetName(), title ) )
                h.SetDirectory(0)
                h.SetMarkerStyle(1)
                h.SetMarkerColor(color)
                if spImpose :
                    self.spimpose[title]=h
                    h.SetFillStyle(0)
                    h.SetLineColor(color)
                    h.SetLineWidth(2)
                else :
                    if color == ROOT.kWhite:
                        h.SetLineColor(ROOT.kGray+3)
                    else:
                        h.SetLineColor(ROOT.TColor.GetColorDark(color))
                    h.SetLineWidth(1)
                    h.SetFillColor(color)
                    h.SetFillStyle(1001)
                    self.mc[title]=h
                self._garbageList.append(h)

    def finalize(self):
        if self.doPoissonErrorBars:
            self.data = convertToPoissonErrorGr(self.dataH)
        else:
            self.data=ROOT.TGraphErrors(self.dataH)

    def appendTo(self,outUrl):
        outF = ROOT.TFile.Open(outUrl,'UPDATE')
        if not outF.cd(self.name):
            outDir = outF.mkdir(self.name)
            outDir.cd()
        for m in self.mc :
            self.mc[m].Write(self.mc[m].GetName(), ROOT.TObject.kOverwrite)
        for m in self.mcsyst :
            self.mcsyst[m].Write(self.mcsyst[m].GetName(), ROOT.TObject.kOverwrite)
        if self.totalMCUnc:
            self.totalMCUnc.Write(self.totalMCUnc.GetName(), ROOT.TObject.kOverwrite)
        for m in self.spimpose:
            self.spimpose[m].Write(self.spimpose[m].GetName(), ROOT.TObject.kOverwrite)
        if self.dataH :
            self.dataH.Write(self.dataH.GetName(), ROOT.TObject.kOverwrite)
        if self.data :
            self.data.Write(self.data.GetName(), ROOT.TObject.kOverwrite)
        outF.Close()

    def reset(self):
        for o in self._garbageList:
            try:
                o.Delete()
            except:
                pass

    def show(self, outDir,lumi,noStack=False,saveTeX=False,extraText=None,noRatio=False):

        if len(self.mc)<1 and self.dataH is None:
            print '%s has 0 or 1 MC!' % self.name
            return

        if len(self.mc)>0 and self.mc.values()[0].InheritsFrom('TH2') :
            print 'Skipping TH2'
            return

        cwid=1000 if self.wideCanvas else 500
        c = ROOT.TCanvas('c','c',cwid,500)
        c.SetBottomMargin(0.0)
        c.SetLeftMargin(0.0)
        c.SetTopMargin(0)
        c.SetRightMargin(0.00)

        #holds the main plot
        c.cd()
        p1 = None
        if self.dataH and not noRatio:
            p1=ROOT.TPad('p1','p1',0.0,0.2,1.0,1.0) if cwid!=1000 else ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.06)
            p1.SetBottomMargin(0.01)
            if self.wideCanvas and len(self.mc)==0 : p1.SetBottomMargin(0.12)
        else:
            p1=ROOT.TPad('p1','p1',0.0,0.0,1.0,1.0)
            p1.SetRightMargin(0.05)
            p1.SetLeftMargin(0.12)
            p1.SetTopMargin(0.08)
            if noRatio:
                p1.SetTopMargin(0.05)
            p1.SetBottomMargin(0.12)
        p1.Draw()

        p1.SetGridx(False)
        p1.SetGridy(False) #True)
        self._garbageList.append(p1)
        p1.cd()

        # legend
        iniy=0.9
        dy=0.0775
        ndy= max(len(self.mc),1)
        inix,dx =0.55,0.4
        if noRatio: inix=0.6
        if noStack:
            inix,dx=0.6,0.35
            iniy,dy,ndy=0.85,0.03,len(self.mc)

        leg = ROOT.TLegend(inix, iniy-dy*ndy, inix+dx, iniy)

        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        #leg.SetTextSize(0.045 if self.wideCanvas else 0.04)
        if noRatio : leg.SetTextSize(0.035)
        nlegCols = 0

        if self.dataH is not None:
            if self.data is None: self.finalize()
            leg.AddEntry( self.data, self.data.GetTitle(),'ep')
            nlegCols += 1
        for h in self.mc:

            #compare
            if noStack:
                refH=self.mc.values()[0]
                if refH!=self.mc[h] and self.doChi2:
                    chi2=refH.Chi2Test( self.mc[h], 'WW CHI2')
                    pval=refH.Chi2Test( self.mc[h], 'WW')
                    self.mc[h].SetTitle('#splitline{%s}{#chi^{2}=%3.1f (p-val: %3.3f)}'%(self.mc[h].GetTitle(),chi2,pval))
                else:
                    refH.SetLineWidth(2)

            leg.AddEntry(self.mc[h], self.mc[h].GetTitle(), 'f')
            nlegCols += 1
        if nlegCols ==0 :
            print '%s is empty'%self.name
            return

        #if not noStack:
        #    leg.SetNColumns(ROOT.TMath.Min(nlegCols/2,3))

        # Build the stack to plot from all backgrounds
        totalMC = None
        stack = ROOT.THStack('mc','mc')
        for h in reversed(self.mc):

            if noStack:
                self.mc[h].SetFillStyle(0)
                self.mc[h].SetLineColor(self.mc[h].GetFillColor())

            stack.Add(self.mc[h],'hist')

            try:
                totalMC.Add(self.mc[h])
            except:
                totalMC = self.mc[h].Clone('totalmc')
                self._garbageList.append(totalMC)
                totalMC.SetDirectory(0)
            nominalTTbar=None
            if h=='t#bar{t}':
                nominalTTbar = self.mc[h].Clone('nomttbar')
                self._garbageList.append(nominalTTbar)
                nominalTTbar.SetDirectory(0)

        #systematics
        if (totalMC and nominalTTbar and len(self.mcsyst)>0):
            # complete
            systUp=[0.]
            systDown=[0.]
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                systUp.append(0.)
                systDown.append(0.)
                for hname in self.mcsyst:
                    if any(x in hname for x in ['fsr', 'Herwig']):
                        continue # these samples should have extra b-tag efficiencies, don't use them to calculate the normalization uncertainty
                    diff = self.mcsyst[hname].GetBinContent(xbin) - nominalTTbar.GetBinContent(xbin)
                    if (diff > 0):
                        systUp[xbin] = math.sqrt(systUp[xbin]**2 + diff**2)
                    else:
                        systDown[xbin] = math.sqrt(systDown[xbin]**2 + diff**2)
            totalMCUnc = totalMC.Clone('totalmcunc')
            self._garbageList.append(totalMCUnc)
            totalMCUnc.SetDirectory(0)
            totalMCUnc.SetFillColor(ROOT.TColor.GetColor('#333333'))
            ROOT.gStyle.SetHatchesLineWidth(1)
            totalMCUnc.SetFillStyle(3254)
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                totalMCUnc.SetBinContent(xbin, totalMCUnc.GetBinContent(xbin) + (systUp[xbin]-systDown[xbin])/2.)
                totalMCUnc.SetBinError(xbin, math.sqrt(totalMCUnc.GetBinError(xbin)**2 + ((systUp[xbin]+systDown[xbin])/2.)**2))
            # shape
            nominalIntegral = nominalTTbar.Integral()
            for hname,h in self.mcsyst.iteritems():
                print(hname,h.Integral())
                if (h.Integral()>0.): h.Scale(nominalIntegral/h.Integral())
            systUpShape=[0.]
            systDownShape=[0.]
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                systUpShape.append(0.)
                systDownShape.append(0.)
                for hname in self.mcsyst:
                    diff = self.mcsyst[hname].GetBinContent(xbin) - nominalTTbar.GetBinContent(xbin)
                    if (diff > 0):
                        systUpShape[xbin] = math.sqrt(systUpShape[xbin]**2 + diff**2)
                    else:
                        systDownShape[xbin] = math.sqrt(systDownShape[xbin]**2 + diff**2)
            totalMCUncShape = totalMC.Clone('totalmcuncshape')
            self._garbageList.append(totalMCUncShape)
            totalMCUncShape.SetDirectory(0)
            totalMCUncShape.SetFillColor(ROOT.TColor.GetColor('#d73027'))
            totalMCUncShape.SetFillStyle(3254)
            for xbin in xrange(1,nominalTTbar.GetNbinsX()+1):
                totalMCUncShape.SetBinContent(xbin, totalMCUncShape.GetBinContent(xbin) + (systUpShape[xbin]-systDownShape[xbin])/2.)
                totalMCUncShape.SetBinError(xbin, math.sqrt(totalMCUncShape.GetBinError(xbin)**2 + ((systUpShape[xbin]+systDownShape[xbin])/2.)**2))
            self.totalMCUnc = totalMCUnc

        #test for null plots
        if totalMC :
            if totalMC.Integral()==0:
                if self.dataH is None : return
                if self.dataH.Integral()==0: return
        elif self.dataH is None : return
        elif self.dataH.Integral()==0 : return


        frame = totalMC.Clone('frame') if totalMC is not None else self.dataH.Clone('frame')
        frame.Reset('ICE')
        if noStack:
            try:
                maxY=self.dataH.GetMaximum()*1.25 
            except:
                maxY=stack.GetStack().At(0).GetMaximum()/1.25
        elif totalMC:
            maxY = totalMC.GetMaximum()
            if self.dataH:
                if maxY<self.dataH.GetMaximum():
                    maxY=self.dataH.GetMaximum()
        else:
            maxY=self.dataH.GetMaximum()

        frame.GetYaxis().SetRangeUser(self.frameMin,self.frameMax*maxY)

        frame.SetDirectory(0)
        frame.Reset('ICE')
        self._garbageList.append(frame)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.045)
        ROOT.TGaxis.SetMaxDigits(4)
        frame.GetYaxis().SetTitleOffset(1.2)
        if noRatio:
            frame.GetYaxis().SetTitleOffset(1.4)
        if self.dataH:
            frame.GetXaxis().SetTitleSize(0.0)
            frame.GetXaxis().SetLabelSize(0.0)
        else :
            frame.GetXaxis().SetTitleSize(0.045)
            frame.GetXaxis().SetLabelSize(0.04)
        if noRatio:
            frame.GetXaxis().SetTitleSize(0.045)
            frame.GetXaxis().SetLabelSize(0.04)
        if self.wideCanvas and totalMC is None :
            frame.GetXaxis().SetLabelSize(0.03)
            frame.GetXaxis().SetTitleSize(0.035)
        frame.Draw()
        if totalMC is not None   :
            if noStack: stack.Draw('nostack same')
            else:
                stack.Draw('hist same')
                if (len(self.mcsyst)>0):
                    totalMCUnc.Draw('e2 same')
                    totalMCUncShape.Draw('e2 same')
                    leg.AddEntry(totalMCUnc, "MC total unc.", 'f')
                    leg.AddEntry(totalMCUncShape, "MC shape unc.", 'f')
        for m in self.spimpose:
            self.spimpose[m].Draw('histsame')
            leg.AddEntry(self.spimpose[m],self.spimpose[m].GetTitle(),'l')
        if self.data is not None : self.data.Draw('P X0 E1')


        if (totalMC is not None and totalMC.GetMaximumBin() > totalMC.GetNbinsX()/2.):
            inix = 0.15
            leg.SetX1(inix)
            leg.SetX2(inix+dx)
        
        leg.Draw()
        
        txt=ROOT.TLatex()
        txt.SetNDC(True)
        txt.SetTextFont(42)
        txt.SetTextSize(0.05)
        txt.SetTextAlign(12)
        iniy=0.88 if self.wideCanvas else 0.88
        inix=0.15
        if (totalMC is not None and totalMC.GetMaximumBin() > totalMC.GetNbinsX()/2.):
            inix = 0.92
        inixlumi=0.66
        if not self.dataH or noRatio:
            inix=0.56
            inixlumi=0.83

        txt.DrawLatex(inix,iniy,'#scale[1.2]{%s}'%self.cmsLabel)
        if lumi<1:
            txt.DrawLatex(inixlumi,0.97,'#scale[1.0]{%3.1f nb^{-1} (%s)}' % (lumi*1000.,self.com) )
        elif lumi<100:
            txt.DrawLatex(inixlumi,0.97,'#scale[1.0]{%3.1f pb^{-1} (%s)}' % (lumi,self.com) )
        else:
            txt.DrawLatex(inixlumi,0.97,'#scale[1.0]{%3.1f fb^{-1} (%s)}' % (lumi/1000.,self.com) )
        txt.DrawLatex(inix,0.815,'#scale[1.0]{t#bar{t} #rightarrow lepton+jets}')
        txt.DrawLatex(inix,0.75,'#scale[1.0]{inclusive jets}')
        txt.DrawLatex(inix,0.68,'#scale[1.0]{p_{T} > 30 GeV}')
        try:
            extraCtr=1
            for extra in extraText.split('\\'):
                txt.DrawLatex(inix,iniy-0.05*extraCtr,'#scale[0.8]{#it{%s}}'%extra)
                extraCtr+=1
        except:
            pass


        #holds the ratio
        c.cd()
        if not noRatio and self.dataH and len(self.mc)>0 :
            p2 = ROOT.TPad('p2','p2',0.0,0.0,1.0,0.2)
            p2.Draw()
            p2.SetBottomMargin(0.4)
            p2.SetRightMargin(0.05)
            p2.SetLeftMargin(0.12)
            p2.SetTopMargin(0.01)
            p2.SetGridx(False)
            p2.SetGridy(True)
            self._garbageList.append(p2)
            p2.cd()
            ratioframe=frame.Clone('ratioframe')
            ratioframe.GetYaxis().SetTitle('MC/data')
            ratioframe.GetYaxis().SetRangeUser(self.ratiorange[0], self.ratiorange[1])
            self._garbageList.append(ratioframe)
            ratioframe.GetYaxis().SetNdivisions(503)
            ratioframe.GetYaxis().SetLabelSize(0.18)
            ratioframe.GetYaxis().SetTitleSize(0.2)
            ratioframe.GetYaxis().SetTitleOffset(0.3)
            ratioframe.GetXaxis().SetLabelSize(0.18)
            ratioframe.GetXaxis().SetTitleSize(0.2)
            ratioframe.GetXaxis().SetTitleOffset(0.8)
            ratioframe.SetFillStyle(3254)
            ratioframe.SetFillColor(ROOT.TColor.GetColor('#333333'))


            #in case we didn't stack compare each distribution
            ratioGrs=[]
            if noStack:

                ratioframe.Draw('e2')

                for title in self.mc:
                    ratio=self.dataH.Clone('ratio')
                    ratio.SetDirectory(0)
                    self._garbageList.append(ratio)
                    ratio.Divide(self.mc[title])
                    ratioGrs.append( ROOT.TGraphAsymmErrors(ratio) )
                    ratioGrs[-1].SetMarkerStyle(self.data.GetMarkerStyle())
                    ratioGrs[-1].SetMarkerSize(self.data.GetMarkerSize())
                    ci=self.mc[title].GetLineColor()
                    ratioGrs[-1].SetMarkerColor(ci)
                    ratioGrs[-1].SetLineColor(ci)
                    ratioGrs[-1].SetLineWidth(self.data.GetLineWidth())
                    ratioGrs[-1].Draw('p')
                    
            #in case we stacked compare only the total
            else:
                if (len(self.mcsyst)>0):
                    ratioframeshape=ratioframe.Clone('ratioframeshape')
                    self._garbageList.append(ratioframeshape)
                    ratioframeshape.SetFillColor(ROOT.TColor.GetColor('#d73027'))            
                
                totalMCnoUnc=totalMC.Clone('totalMCnounc')
                self._garbageList.append(totalMCnoUnc)
                for xbin in xrange(1,totalMC.GetNbinsX()+1):
                    ratioframe.SetBinContent(xbin,1)
                    val=self.dataH.GetBinContent(xbin)
                    totalMCnoUnc.SetBinError(xbin,0.)
                    if val==0 : continue
                    if (len(self.mcsyst)>0):
                        totalUnc=ROOT.TMath.Sqrt((totalMCUnc.GetBinError(xbin)/val)**2+self.mcUnc**2)
                        totalUncShape=ROOT.TMath.Sqrt((totalMCUncShape.GetBinError(xbin)/val)**2+self.mcUnc**2)
                        ratioframeshape.SetBinContent(xbin,totalMCUncShape.GetBinContent(xbin)/val)
                        ratioframeshape.SetBinError(xbin,totalUncShape)
                        limitToRange(ratioframeshape, self.ratiorange)
                    else: totalUnc=ROOT.TMath.Sqrt((totalMC.GetBinError(xbin)/val)**2+self.mcUnc**2)
                    ratioframe.SetBinContent(xbin,totalMCnoUnc.GetBinContent(xbin)/val)
                    ratioframe.SetBinError(xbin,totalUnc)
                    limitToRange(ratioframe, self.ratiorange)
                ratioframe.Draw('e2')
                if (len(self.mcsyst)>0): ratioframeshape.Draw('e2 same')
                #try:
                ratio=self.dataH.Clone('ratio')
                ratio.SetDirectory(0)
                self._garbageList.append(ratio)
                for xbin in xrange(1,ratio.GetNbinsX()+1):
                    if ratio.GetBinContent(xbin) > 0.:
                        ratio.SetBinError(xbin, ratio.GetBinError(xbin)/ratio.GetBinContent(xbin))
                        ratio.SetBinContent(xbin, 1.)
                    else:
                        ratio.SetBinError  (xbin, 0.)
                        ratio.SetBinContent(xbin, 0.)
                ratio.SetMarkerSize(0)
                ratio.SetFillStyle(3001)
                ratio.SetFillColor(ROOT.kBlack)
                ratio.SetLineColor(ROOT.kBlack)
                ratio.Draw('e2 same')
                ratio.Draw('same l')
                #gr=ROOT.TGraphAsymmErrors(ratio)
                #gr.SetMarkerStyle(self.data.GetMarkerStyle())
                #gr.SetMarkerSize(self.data.GetMarkerSize())
                #gr.SetMarkerColor(self.data.GetMarkerColor())
                #gr.SetLineColor(self.data.GetLineColor())
                #gr.SetLineWidth(self.data.GetLineWidth())
                #gr.Draw('p')
                #except:
                #    pass
                
                # Let's compute a chi2 from the ratios
                #if (len(self.mcsyst)>0):
                #    c.cd()
                #    txt.DrawLatex(inix,iniy-0.02,'#scale[0.5]{#chi^{2}/ndf = %0.2f}'%(ratio.Chi2Test(ratioframeshape, 'UW,CHI2/NDF')))

        #all done
        if p1: p1.RedrawAxis()
        c.cd()
        c.Modified()
        c.Update()

        #save
        for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'.'+ext))
        if self.savelog:
            p1.cd()
            frame.GetYaxis().SetRangeUser(1,maxY*50)
            p1.SetLogy()
            c.cd()
            c.Modified()
            c.Update()
            for ext in self.plotformats : c.SaveAs(os.path.join(outDir, self.name+'_log.'+ext))

        if saveTeX : self.convertToTeX(outDir=outDir)


    def convertToTeX(self, outDir):
        if len(self.mc)==0:
            print '%s is empty' % self.name
            return

        f = open(outDir+'/'+self.name+'.dat','w')
        f.write('------------------------------------------\n')
        f.write("Process".ljust(20),)
        f.write("Events after each cut\n")
        f.write('------------------------------------------\n')

        tot ={}
        err = {}
        f.write(' '.ljust(20),)
        try:
            for xbin in xrange(1,self.mc.values()[0].GetXaxis().GetNbins()+1):
                pcut=self.mc.values()[0].GetXaxis().GetBinLabel(xbin)
                f.write(pcut.ljust(40),)
                tot[xbin]=0
                err[xbin]=0
        except:
            pass
        f.write('\n')
        f.write('------------------------------------------\n')

        for pname in self.mc:
            h = self.mc[pname]
            f.write(pname.ljust(20),)

            for xbin in xrange(1,h.GetXaxis().GetNbins()+1):
                itot=h.GetBinContent(xbin)
                ierr=h.GetBinError(xbin)
                pval=' & %s'%toLatexRounded(itot,ierr)
                f.write(pval.ljust(40),)
                tot[xbin] = tot[xbin]+itot
                err[xbin] = err[xbin]+ierr*ierr
            f.write('\n')

        f.write('------------------------------------------\n')
        f.write('Total'.ljust(20),)
        for xbin in tot:
            pval=' & %s'%toLatexRounded(tot[xbin],math.sqrt(err[xbin]))
            f.write(pval.ljust(40),)
        f.write('\n')

        if self.dataH is None: return
        f.write('------------------------------------------\n')
        f.write('Data'.ljust(20),)
        for xbin in xrange(1,self.dataH.GetXaxis().GetNbins()+1):
            itot=self.dataH.GetBinContent(xbin)
            pval=' & %d'%itot
            f.write(pval.ljust(40))
        f.write('\n')
        f.write('------------------------------------------\n')
        f.close()



"""
converts a histogram to a graph with Poisson error bars
"""
def convertToPoissonErrorGr(h):
    try:
        htype=h.ClassName()
    except:
        return None
    if htype.find('TH1')<0 : return None

    #check https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
    alpha = 1 - 0.6827;
    grpois = ROOT.TGraphAsymmErrors(h);
    for i in xrange(0,grpois.GetN()+1) :
        N = grpois.GetY()[i]
        if N<200 :
            L = 0
            if N>0 : L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
            U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1)
            grpois.SetPointEYlow(i, N-L)
            grpois.SetPointEYhigh(i, U-N)
        else:
            grpois.SetPointEYlow(i, math.sqrt(N))
            grpois.SetPointEYhigh(i,math.sqrt(N))
    return grpois
