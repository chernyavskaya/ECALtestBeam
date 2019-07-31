import ROOT
from ROOT import RooRealVar,RooCBShape,RooDataHist,RooArgList,RooFit
import csv
import json
import array as array

class CBfunction:
    
    nbins = 600
    xmin = 0
    xmax = 12000
    xaxis_scale = 0.2
    a_initial = 0.5
    n_initial = 7

    def __init__(self,data):
        self.data = data ## data should be TChain already
      
    def set_crystal(self,crystal):
        self.crystal = crystal
    def set_energy(self,energy):
        self.energy = energy
    def set_position(self,x=0,y=0,window=20):
        self.xcenter = x
        self.ycenter = y
        self.window = window
        
    def set_selection(self):
        self.selection = "n_tracks==1 && fabs(X-(%.2f))<%.2f && fabs(Y-(%.2f))<%.2f"%(self.xcenter,self.window,self.ycenter,self.window)

    def prepare_sumhistogram(self,dict_crystals_calibration,matrix):
        self.set_selection()
        self.hist = ROOT.TH1F("ampl_%s_%s"%(self.crystal,self.energy),"ampl_%s_%s"%(self.crystal,self.energy),self.nbins,self.xmin,self.xmax)

        draw_function = '('
        for enum,cryst in enumerate(matrix):
            #draw_function +='fit_ampl[%s]*%.4f'%(cryst,dict_crystals_calibration[cryst])
            draw_function +='((fit_ampl[%s] > 0) ? fit_ampl[%s] : amp_max[%s])*%.4f'%(cryst,cryst,cryst,dict_crystals_calibration[cryst])
            if enum!=len(matrix)-1 : draw_function+='+'
          #  else : draw_function+=")*%.4f>>ampl_%s_%s"%(dict_crystals_calibration['conversion_factor'],self.crystal,self.energy)
            else : draw_function+=")>>ampl_%s_%s"%(self.crystal,self.energy)
   
    
        self.data.Draw(draw_function,self.selection,"goff")
        self.peak_position = self.hist.GetXaxis().GetBinCenter(self.hist.GetMaximumBin())
        self.ymax_value = self.hist.GetMaximum()        
        
        
    def prepare_histogram(self):
        self.set_selection()
        self.hist = ROOT.TH1F("ampl_%s_%s"%(self.crystal,self.energy),"ampl_%s_%s"%(self.crystal,self.energy),self.nbins,self.xmin,self.xmax)        
        self.data.Draw("fit_ampl[%s]>>ampl_%s_%s"%(self.crystal,self.crystal,self.energy),self.selection,"goff")
        self.peak_position = self.hist.GetXaxis().GetBinCenter(self.hist.GetMaximumBin())
        self.ymax_value = self.hist.GetMaximum()
        
    def CBintialization(self):
        round_energy = round(float(self.energy),-1)
        self.x = RooRealVar("signal_%s_%dGeV"%(self.crystal,round_energy),"signal_%s_%dGeV"%(self.crystal,round_energy),max(0.,self.peak_position*(1-self.xaxis_scale)),self.peak_position*(1+self.xaxis_scale))
        self.roohist = RooDataHist("roohist_fit_%s_%s"%(self.crystal,self.energy),"roohist_fit_%s_%s"%(self.crystal,self.energy),RooArgList(self.x),self.hist)
        self.m = RooRealVar("mean_%s_%s"%(self.crystal,self.energy),"mean_%s_%s"%(self.crystal,self.energy),self.peak_position,max(0.,self.peak_position*(1-self.xaxis_scale)),self.peak_position*(1+self.xaxis_scale))
        self.s = RooRealVar("sigma_%s_%s"%(self.crystal,self.energy),"sigma_%s_%s"%(self.crystal,self.energy),60,0,500)
        self.a = RooRealVar("alpha_%s_%s"%(self.crystal,self.energy),"alpha_%s_%s"%(self.crystal,self.energy),self.a_initial,-20.,20)
        self.n = RooRealVar("exp_%s_%s"%(self.crystal,self.energy),"exp_%s_%s"%(self.crystal,self.energy),self.n_initial,0.,150)
        self.sig = RooCBShape("signal_%s_%s"%(self.crystal,self.energy),"signal_%s_%s"%(self.crystal,self.energy),self.x,self.m,self.s,self.a,self.n)

    def fitToData(self):
        self.res = self.sig.fitTo(self.roohist)    
        
    def fitResults(self):
        self.dict_fit_results = {}
        self.dict_fit_results['CBmean'] = [self.m.getVal(),self.m.getError()]
        self.dict_fit_results['CBsigma'] = [self.s.getVal(),self.s.getError()]
        self.dict_fit_results['CBalpha'] = [self.a.getVal(),self.a.getError()]
        self.dict_fit_results['CBexp'] = [self.n.getVal(),self.n.getError()]
        self.dict_fit_results['chi2'] = self.chi2
        return self.dict_fit_results


    def plot(self):
        self.frame = self.x.frame()
        self.roohist.plotOn(self.frame,RooFit.Name("roohist_chi2_%s_%s"%(self.crystal,self.energy)))
        self.sig.plotOn(self.frame,RooFit.Name("signal_chi2_%s_%s"%(self.crystal,self.energy)))
        self.chi2 = self.frame.chiSquare("signal_chi2_%s_%s"%(self.crystal,self.energy),"roohist_chi2_%s_%s"%(self.crystal,self.energy),4) # 4 = nFitParameters from CB
        self.sig.paramOn(self.frame,RooFit.Layout(0.65,0.99,0.8))
        self.frame.getAttText().SetTextSize(0.03)

        txt_chi2 = ROOT.TText(self.peak_position*0.83,self.ymax_value*0.95,"Chi2 = %.1f"%self.chi2)
        txt_chi2.SetTextSize(0.04)
        txt_chi2.SetTextColor(ROOT.kRed)
        self.frame.addObject(txt_chi2)
        self.frame.Draw()
