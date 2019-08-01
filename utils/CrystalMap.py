import ROOT
from ROOT import TH1F
from ROOT import gROOT
import csv
import json
import array as array
from ROOT import gSystem,gPad

class CrystalMap:
        
    nbins = 30
    xmin = -15
    xmax = 15
    ymin = -15
    ymax = 15

    def __init__(self,data):
        self.data = data ## data should be TChain already
        self.grmax = 0
        self.grmax_testBeam = 0
        
    def set_crystal(self,crystal):
        self.crystal = crystal
    def set_energy(self,energy):
        self.energy = energy   
    def set_selection(self):
        self.selection = "n_tracks==1"


    def plot(self,dict_crystal_centers):
        round_energy = round(float(self.energy),-1)
        if round_energy ==240 : round_energy = 250
        
        self.data.Draw("fit_ampl[%s]:Y:X>>hhnew2_%s_%s(%d,%d,%d,%d,%d,%d,0,10000)"%(self.crystal,self.crystal,self.energy,self.nbins,self.xmin,self.xmax,self.nbins,self.ymin,self.ymax),self.selection,"PROFCOLZ")
        self.htemp = ROOT.gPad.GetPrimitive("hhnew2_%s_%s"%(self.crystal,self.energy))
        self.htemp.GetXaxis().SetRange(self.htemp.GetXaxis().FindBin(self.xmin),self.htemp.GetXaxis().FindBin(self.xmax))
        self.htemp.GetYaxis().SetRange(self.htemp.GetYaxis().FindBin(self.ymin),self.htemp.GetYaxis().FindBin(self.ymax))
        xx, yy, zz = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
        self.htemp.GetBinXYZ(self.htemp.GetMaximumBin(), xx,yy,zz);
        xmean =  self.htemp.GetXaxis().GetBinCenter(xx)
        ymean = self.htemp.GetYaxis().GetBinCenter(yy)
        
        
        self.data.Draw("fit_ampl[%s]:Y:X>>map_%s_%s(30,-15,15,30,-15,15,0,10000)"%(self.crystal,self.crystal,self.energy),self.selection,"PROFCOLZ")
        self.hist2d = ROOT.gPad.GetPrimitive("map_%s_%s"%(self.crystal,self.energy))
        self.hist2d.SetTitle("%s %d GeV"%(self.crystal,round_energy))
        self.hist2d.GetXaxis().SetTitle('X (mm)')
        self.hist2d.GetYaxis().SetTitle('Y (mm)')
        self.hist2d.Draw("PROFCOLZ")   

        self.grmax = TH1F('gr_%s_%s'%(self.crystal,self.energy),'gr_%s_%s'%(self.crystal,self.energy),1,xmean-1,xmean+1)
        self.grmax.SetBinContent(1,ymean)
        self.grmax.SetMarkerColor(ROOT.kMagenta+2)
        self.grmax.SetMarkerStyle(29)
        self.grmax.Draw("LPsame")
        ROOT.gDirectory.Append(self.grmax)

        if bool(dict_crystal_centers):
            crystal_meanx = dict_crystal_centers[self.crystal][0]
            crystal_meany = dict_crystal_centers[self.crystal][1]

            self.grmax_testBeam = TH1F('grTestBeam_%s_%s'%(self.crystal,self.energy),'grTestBeam_%s_%s'%(self.crystal,self.energy),1,crystal_meanx-1,crystal_meanx+1)
            self.grmax_testBeam.SetBinContent(1,crystal_meany)
            self.grmax_testBeam.SetMarkerColor(2)
            self.grmax_testBeam.SetMarkerStyle(29)
            self.grmax_testBeam.Draw("LPsame")
            ROOT.gDirectory.Append(self.grmax_testBeam)

        

