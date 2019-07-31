import ROOT
from ROOT import gROOT
import csv
import json
import array as array
from ROOT import gSystem

class CrystalMap:
        
    nbins = 6
    xmin = -15
    xmax = 15
    ymin = -15
    ymax = 15

    def __init__(self,data):
        self.data = data ## data should be TChain already
        
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
        self.htemp2 = ROOT.gPad.GetPrimitive("hhnew2_%s_%s"%(self.crystal,self.energy))
        xmean =  self.htemp2.GetMean(1)
        ymean = self.htemp2.GetMean(2)
        print xmean,ymean
        
        self.data.Draw("fit_ampl[%s]:Y:X>>map_%s_%s(15,-15,15,15,-15,15,0,10000)"%(self.crystal,self.crystal,self.energy),self.selection,"PROFCOLZ")
        self.hist2d = ROOT.gPad.GetPrimitive("map_%s_%s"%(self.crystal,self.energy))
        self.hist2d.SetTitle("%s %d GeV"%(self.crystal,round_energy))
        self.hist2d.GetXaxis().SetTitle('X (mm)')
        self.hist2d.GetYaxis().SetTitle('Y (mm)')
        self.hist2d.Draw("PROFCOLZsame")   

        self.grmax = ROOT.TGraph(1,array.array( 'd' ,[xmean]),array.array( 'd' ,[ymean]) )
        self.grmax.SetName('gr_%s_%s'%(self.crystal,self.energy))
        self.grmax.SetMarkerColor(ROOT.kMagenta+2)
        self.grmax.SetMarkerStyle(29)
        self.grmax.Draw("LPsame")
        
        self.grmax_testBeam = ROOT.TGraph(1,array.array( 'd' ,[dict_crystal_centers[self.crystal][0]]),array.array( 'd' ,[dict_crystal_centers[self.crystal][1]]))
        self.grmax_testBeam.SetName('gr_testBeam_%s_%s'%(self.crystal,self.energy))
        self.grmax_testBeam.SetMarkerColor(2)
        self.grmax_testBeam.SetMarkerStyle(29)
        self.grmax_testBeam.Draw("LPsame")

