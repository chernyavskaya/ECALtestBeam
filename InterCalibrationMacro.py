
# coding: utf-8

# In[1]:

import ROOT
from ROOT import RooRealVar,RooCBShape,RooDataHist,RooArgList,RooFit
import csv
import json
import array as array


# In[2]:

reader = csv.reader(open('data/intercalibration19C.csv', 'r'))
dict_run_cryst = {}
for row in reader:
    run,energy,crystal,_ = row
    if energy=='149.12' :
        dict_run_cryst[crystal] = run
    
##Now add runs for C2,C3,C4 
reader = csv.reader(open('data/energyscans19C.csv', 'r'))
for row in reader:
    run,energy,crystal,_,_,_ = row
    if energy=='149.12' and  (crystal=='C3' or crystal=='C2' or crystal=='C4' ) :
        dict_run_cryst[crystal] = run
        
## Read crystal centers prepared by Simone 
reader = csv.reader(open('data/crystalscenters.csv', 'r'))
dict_crystal_centers={}
for row in reader:
    crystal,xpos,ypos = row
    if crystal!='Crystal' :
        dict_crystal_centers[crystal] = [int(xpos),int(ypos)]



# In[3]:

naming = []
for i in range(5,0,-1):
    naming.append('A%d'%i)
    naming.append('B%d'%i)
    naming.append('C%d'%i)
    naming.append('D%d'%i)
    naming.append('E%d'%i)
missing = 'A5,A4,E5,E4'.split(',')
working = [item for item in naming]
for item in missing:
    working.remove(item)


# In[ ]:




# In[ ]:




# In[ ]:

#crystal_maps=[]
#c = ROOT.TCanvas("c","c",900,900)
#c.Divide(5,5)
#canvas_num=0
#files=[]
#grmax=[]
#grmaxSimone=[]
#file_num = 0
#for crystal in naming:
#    if crystal in missing:
#        canvas_num+=1
#    else :
#        c.cd(canvas_num+1)
#    
#        run = dict_run_cryst[crystal]
#        files.append(ROOT.TFile("/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2018/ntuples_v5/ECAL_H4_October2018_%s.root"%run))
#        tree = files[file_num].Get("h4")
#        
#        tree.Draw("fit_ampl[%s]:Y:X>>hhnew2_%s(6,-15,15,6,-15,15,0,10000)"%(crystal,crystal),"n_tracks==1","PROFCOLZ")
#        htemp2 = ROOT.gPad.GetPrimitive("hhnew2_%s"%crystal)
#        xmax =  htemp2.GetMean(1)
#        ymax = htemp2.GetMean(2)
#        
#        tree.Draw("fit_ampl[%s]:Y:X>>hhnew_%s(15,-15,15,15,-15,15,0,10000)"%(crystal,crystal),"n_tracks==1","PROFCOLZ")
#        hist2d = ROOT.gPad.GetPrimitive("hhnew_%s"%crystal)
#        hist2d.SetTitle(crystal)
#        hist2d.GetXaxis().SetTitle('X (mm)')
#        hist2d.GetYaxis().SetTitle('Y (mm)')
#        hist2d.Draw("PROFCOLZsame")
#        crystal_maps.append(hist2d)
#    
#        
        
 #       hist2d.GetXaxis().SetRange(hist2d.GetXaxis().FindBin(-10),hist2d.GetXaxis().FindBin(10))
#        hist2d.GetYaxis().SetRange(hist2d.GetYaxis().FindBin(-10),hist2d.GetYaxis().FindBin(10))
        
#        hist2d.GetXaxis().SetRange(0,hist2d.GetNbinsX())
#        hist2d.GetYaxis().SetRange(0,hist2d.GetNbinsY())
#
#        grmax.append(ROOT.TGraph(1,array.array( 'd' ,[xmax]),array.array( 'd' ,[ymax])))
#        grmax[file_num].SetMarkerColor(ROOT.kMagenta+2);
#        grmax[file_num].SetMarkerStyle(29);
#        grmax[file_num].Draw("LPsame")
#        
#        grmaxSimone.append(ROOT.TGraph(1,array.array( 'd' ,[dict_crystal_centers[crystal][0]]),array.array( 'd' ,[dict_crystal_centers[crystal][1]])))
#        grmaxSimone[file_num].SetMarkerColor(2);
#        grmaxSimone[file_num].SetMarkerStyle(29);
#        grmaxSimone[file_num].Draw("LPsame")
#        
#        file_num+=1
#        canvas_num+=1
#c.Draw()
#c.SaveAs('plots/crystal_map.pdf')
#c.SaveAs('plots/crystal_map.png')


# In[ ]:

pos_cut = 5

c = ROOT.TCanvas("c","c",900,900)
c.Divide(5,5)
canvas_num=0
files=[]
file_num=0
chi2s = []
results = []
hists_fits=[]
trees=[]
x,roohist,m,s,a,n,sig,res,frame = [],[],[],[],[],[],[],[],[]
dict_crystals_means = {}
#for crystal_num in range(0,len(naming)):
#for crystal_num in range(0,3):
for crystal in naming:
#for crystal in ['A3','C3','B3']:
    if crystal in missing:
        canvas_num+=1
    else :
        c.cd(canvas_num+1) 
        print canvas_num
        
        run = dict_run_cryst[crystal]
        files.append(ROOT.TFile("/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2018/ntuples_v5/ECAL_H4_October2018_%s.root"%run))
        trees.append(files[file_num].Get("h4"))
    
        hist_ampl = ROOT.TH1F("hist_ampl_%s"%crystal,"hist_ampl_%s"%crystal,300,2000,8000)
        trees[file_num].Draw("fit_ampl[%s]>>hist_ampl_%s"%(crystal,crystal),"n_tracks==1 && fabs(X-(%.2f))<%.2f && fabs(Y-(%.2f))<%.2f"%(dict_crystal_centers[crystal][0],dict_crystal_centers[crystal][1],pos_cut,pos_cut),"")
        peak_position = hist_ampl.GetXaxis().GetBinCenter(hist_ampl.GetMaximumBin())
        ymax_value = hist_ampl.GetMaximum()

    
        x.append(RooRealVar("fit_ampl_%s"%crystal,"fit_ampl_%s"%crystal,peak_position*0.8,peak_position*1.2))
        roohist.append(RooDataHist("roohist_fit_%s"%crystal,"roohist_fit_%s"%crystal,RooArgList(x[file_num]),hist_ampl))

        m.append(RooRealVar("mean_%s"%crystal,"mean_%s"%crystal,peak_position,peak_position*0.8,peak_position*1.2))
        s.append(RooRealVar("sigma_%s"%crystal,"sigma_%s"%crystal,60,0,500))
        a.append(RooRealVar("alpha_%s"%crystal,"alpha_%s"%crystal,0.5,0,10))
        n.append(RooRealVar("exp_%s"%crystal,"exp_%s"%crystal,7,0,150))
        sig.append(RooCBShape("signal_%s"%crystal,"signal_%s"%crystal,x[file_num],m[file_num],s[file_num],a[file_num],n[file_num]))

        res.append(sig[file_num].fitTo(roohist[file_num]))#RooFit.Save())
        #res.Print()

        frame.append(x[file_num].frame())
        roohist[file_num].plotOn(frame[file_num],RooFit.Name("roohist_%s"%crystal))
        sig[file_num].plotOn(frame[file_num],RooFit.Name("signal_%s"%crystal))
        hists_fits.append(hist_ampl)
        chi2s.append(frame[file_num].chiSquare("signal_%s"%crystal,"roohist_%s"%crystal,4)) # 4 = nFitParameters from CB

        sig[file_num].paramOn(frame[file_num],RooFit.Layout(0.65,0.99,0.8)) 
        frame[file_num].getAttText().SetTextSize(0.03)

        txt_chi2 = ROOT.TText(peak_position*0.83,ymax_value*0.95,"Chi2 = %.1f"%chi2s[file_num])
        txt_chi2.SetTextSize(0.04) 
        txt_chi2.SetTextColor(ROOT.kRed)
        frame[file_num].addObject(txt_chi2)
 
        frame[file_num].Draw()

        dict_crystals_means[crystal] = ['CBmean',[m[file_num].getVal(),m[file_num].getError()],['CBsigma',s[file_num].getVal(),s[file_num].getError()]]

        canvas_num+=1
        file_num+=1
c.Draw()
c.SaveAs('plots/fits.pdf')
c.SaveAs('plots/fits.png')

with open('output/intercarlibration_results.json', 'w') as fp:
    json.dump(dict_crystals_means, fp)
