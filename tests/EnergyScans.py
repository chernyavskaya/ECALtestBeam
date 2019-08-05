import ROOT
from ROOT import RooRealVar,RooCBShape,RooDataHist,RooArgList,RooFit
import csv
import json
import array as array

with open('output/intercarlibration_results.json') as json_file:
    fit_results = json.load(json_file)

dict_crystals_calibration = {}
for key in fit_results:
    dict_crystals_calibration[key] = fit_results['C3'][1:2][0][0]/fit_results[key][1:2][0][0]

## Read crystal centers prepared by Simone 
reader = csv.reader(open('data/crystalscenters.csv', 'r'))
dict_crystal_centers={}
for row in reader:
    crystal,xpos,ypos = row
    if crystal!='Crystal' :
        dict_crystal_centers[crystal] = [int(xpos),int(ypos)]



##Read C3 energy scans
dict_C3_energy_scan = {}
reader = csv.reader(open('data/energyscans19C.csv', 'r'))
for row in reader:
    run,energy,crystal,_,_,_ = row
    if crystal=='C3' :
        if energy in dict_C3_energy_scan.keys() :
            dict_C3_energy_scan[energy].append(run)
        else : 
            dict_C3_energy_scan[energy] = []
            dict_C3_energy_scan[energy].append(run)


pos_cut = 1

crystal_maps=[]
c = ROOT.TCanvas("c","c",2000,400)
c.Divide(5,1)
canvas_num=0
files=[]
file_num = 0
chi2s = []
results = []
hists_fits=[]
trees=[]
x,roohist,m,s,a,n,sig,res,frame = [],[],[],[],[],[],[],[],[]
dict_energy_means={}
crystal='C3'
energies = sorted([float(item) for item in dict_C3_energy_scan.keys()])
energies = [str(item) for item in energies]
for energy in energies :
        c.cd(canvas_num+1)
    
        runs = dict_C3_energy_scan[energy]
        trees.append(ROOT.TChain("h4"))
        for run in runs:
            trees[file_num].Add("/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2018/ntuples_v5/ECAL_H4_October2018_%s.root"%run)

        hist_ampl = ROOT.TH1F("hist_ampl_%s"%energy,"hist_ampl_%s"%energy,500,0,10000)
        trees[file_num].Draw("fit_ampl[%s]>>hist_ampl_%s"%(crystal,energy),"n_tracks==1 && fabs(X-(%.2f))<%.2f && fabs(Y-(%.2f))<%.2f"%(dict_crystal_centers[crystal][0],dict_crystal_centers[crystal][1],pos_cut,pos_cut),"")
        peak_position = hist_ampl.GetXaxis().GetBinCenter(hist_ampl.GetMaximumBin())
        ymax_value = hist_ampl.GetMaximum()


        x.append(RooRealVar("signal_%dGeV"%round(float(energy),-1),"signal_%dGeV"%round(float(energy),-1),peak_position*0.8,peak_position*1.2))
        roohist.append(RooDataHist("roohist_fit_%s"%energy,"roohist_fit_%s"%energy,RooArgList(x[file_num]),hist_ampl))

        m.append(RooRealVar("mean_%s"%energy,"mean_%s"%energy,peak_position,peak_position*0.8,peak_position*1.2))
        s.append(RooRealVar("sigma_%s"%energy,"sigma_%s"%energy,60,0,500))
        a.append(RooRealVar("alpha_%s"%energy,"alpha_%s"%energy,0.5,0,10))
        n.append(RooRealVar("exp_%s"%energy,"exp_%s"%energy,7,0,150))
        sig.append(RooCBShape("signal_%s"%energy,"signal_%s"%energy,x[file_num],m[file_num],s[file_num],a[file_num],n[file_num]))

        res.append(sig[file_num].fitTo(roohist[file_num]))#RooFit.Save())
        #res.Print()

        frame.append(x[file_num].frame())
        roohist[file_num].plotOn(frame[file_num],RooFit.Name("roohist_%s"%energy))
        sig[file_num].plotOn(frame[file_num],RooFit.Name("signal_%s"%energy))
        hists_fits.append(hist_ampl)
        chi2s.append(frame[file_num].chiSquare("signal_%s"%energy,"roohist_%s"%energy,4)) # 4 = nFitParameters from CB

        sig[file_num].paramOn(frame[file_num],RooFit.Layout(0.65,0.99,0.8))
        frame[file_num].getAttText().SetTextSize(0.027)

        txt_chi2 = ROOT.TText(peak_position*0.83,ymax_value*0.95,"Chi2 = %.1f"%chi2s[file_num])
        txt_chi2.SetTextSize(0.04)
        txt_chi2.SetTextColor(ROOT.kRed)
        frame[file_num].addObject(txt_chi2)

        frame[file_num].Draw()

        dict_energy_means[energy] = ['CBmean',[m[file_num].getVal(),m[file_num].getError()],['CBsigma',s[file_num].getVal(),s[file_num].getError()]]

        canvas_num+=1
        file_num+=1
c.Draw()
c.SaveAs('plots/C3_fits_cut%dmm.pdf'%pos_cut)
c.SaveAs('plots/C3_fits_cut%dmm.png'%pos_cut)


with open('output/energyscans_C3_results.json', 'w') as fp:
    json.dump(dict_energy_means, fp)
