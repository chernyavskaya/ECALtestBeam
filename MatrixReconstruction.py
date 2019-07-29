import ROOT
from ROOT import RooRealVar,RooCBShape,RooDataHist,RooArgList,RooFit
import csv
import json
import array as array

with open('output/conversion_factor.json') as json_file:
  	conversion_results = json.load(json_file)

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
matrix_3 = 'B4,C4,D4,B3,C3,D3,B2,C2,D2'.split(',')
matrix_5 = working 

pos_cut = 1
c = ROOT.TCanvas("c","c",900,900)
c.Divide(5,5)
canvas_num=0
files=[]
file_num = 0
chi2s = []
results = []
hists_fits=[]
trees=[]
x,roohist,m,s,a,n,sig,res,frame = [],[],[],[],[],[],[],[],[]
dict_energy_means={}
energies = sorted([float(item) for item in dict_C3_energy_scan.keys()])
energies = [str(item) for item in energies]
sum_energies_3 = {}
sum_energies_5 = {}
##Read C3 energy scans
for energy in energies:
     round_energy = round(float(energy),-1)
     sum_energies_3[energy] = 0.
     sum_energies_5[energy] = 0.
     for crystal in naming :
       if crystal in missing:
        canvas_num+=1
       else :
        c.cd(canvas_num+1)
    
        runs = dict_C3_energy_scan[energy]
        trees.append(ROOT.TChain("h4"))
        for run in runs:
            trees[file_num].Add("/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2018/ntuples_v5/ECAL_H4_October2018_%s.root"%run)

        nbins=500
        if 'crystal'!='C3' : nbins=1000
        hist_ampl = ROOT.TH1F("hist_ampl_%s_%s"%(crystal,energy),"hist_ampl_%s_%s"%(crystal,energy),nbins,0,10000)
        trees[file_num].Draw("fit_ampl[%s]>>hist_ampl_%s_%s"%(crystal,crystal,energy),"n_tracks==1 && fabs(X-(%.2f))<%.2f && fabs(Y-(%.2f))<%.2f"%(dict_crystal_centers[crystal][0],dict_crystal_centers[crystal][1],pos_cut,pos_cut),"")
        peak_position = hist_ampl.GetXaxis().GetBinCenter(hist_ampl.GetMaximumBin())
        mean_position = hist_ampl.GetMean()
        ymax_value = hist_ampl.GetMaximum()

        scale = 0.2
        if crystal!='C3' : scale  = 1.3
        x.append(RooRealVar("signal_%s_%dGeV"%(crystal,round_energy),"signal_%s_%dGeV"%(crystal,round_energy),max(0.,peak_position*(1-scale)),peak_position*(1+scale)))
        roohist.append(RooDataHist("roohist_fit_%s_%s"%(crystal,energy),"roohist_fit_%s_%s"%(crystal,energy),RooArgList(x[file_num]),hist_ampl))

        m.append(RooRealVar("mean_%s_%s"%(crystal,energy),"mean_%s_%s"%(crystal,energy),peak_position,peak_position*0.8,peak_position*1.2))
        s.append(RooRealVar("sigma_%s_%s"%(crystal,energy),"sigma_%s_%s"%(crystal,energy),60,0,500))
        alpha_initial = 0.5
        n_initial = 7
        if crystal!='C3' : 
            alpha_initial = -1.
            n_initial = 0.5
        a.append(RooRealVar("alpha_%s_%s"%(crystal,energy),"alpha_%s_%s"%(crystal,energy),alpha_initial,-10.,10))
        n.append(RooRealVar("exp_%s_%s"%(crystal,energy),"exp_%s_%s"%(crystal,energy),n_initial,0.,150))
        sig.append(RooCBShape("signal_%s_%s"%(crystal,energy),"signal_%s_%s"%(crystal,energy),x[file_num],m[file_num],s[file_num],a[file_num],n[file_num]))

        res.append(sig[file_num].fitTo(roohist[file_num]))#RooFit.Save())
        #res.Print()

        frame.append(x[file_num].frame())
        roohist[file_num].plotOn(frame[file_num],RooFit.Name("roohist_%s_%s"%(crystal,energy)))
        sig[file_num].plotOn(frame[file_num],RooFit.Name("signal_%s_%s"%(crystal,energy)))
        hists_fits.append(hist_ampl)
        chi2s.append(frame[file_num].chiSquare("signal_%s_%s"%(crystal,energy),"roohist_%s_%s"%(crystal,energy),4)) # 4 = nFitParameters from CB

        sig[file_num].paramOn(frame[file_num],RooFit.Layout(0.65,0.99,0.8))
        frame[file_num].getAttText().SetTextSize(0.027)

        txt_chi2 = ROOT.TText(peak_position*((1-scale)+0.3),ymax_value*0.95,"Chi2 = %.1f"%chi2s[file_num])
        txt_chi2.SetTextSize(0.04)
        txt_chi2.SetTextColor(ROOT.kRed)
        frame[file_num].addObject(txt_chi2)

        frame[file_num].Draw()

        dict_energy_means[energy] = ['CBmean',[m[file_num].getVal(),m[file_num].getError()],['CBsigma',s[file_num].getVal(),s[file_num].getError()]]

        if crystal in matrix_3:
            sum_energies_3[energy]+=m[file_num].getVal()*conversion_results[crystal]*conversion_results['conversion_factor']
            #if crystal=='C3' : sum_energies_3[energy]+=m[file_num].getVal()*conversion_results[crystal]*conversion_results['conversion_factor']
            #else : sum_energies_3[energy]+=peak_position*conversion_results[crystal]*conversion_results['conversion_factor']
        if crystal in matrix_5:
            sum_energies_5[energy]+=m[file_num].getVal()*conversion_results[crystal]*conversion_results['conversion_factor']
            #if crystal=='C3' : sum_energies_5[energy]+=m[file_num].getVal()*conversion_results[crystal]*conversion_results['conversion_factor']
            #else : sum_energies_5[energy]+=peak_position*conversion_results[crystal]*conversion_results['conversion_factor']

        canvas_num+=1
        file_num+=1
     c.Draw()
     c.SaveAs('plots/sum_energy_%s_cut%dmm.pdf'%(round_energy,pos_cut))
     c.SaveAs('plots/sum_energy_%s_cut%dmm.png'%(round_energy,pos_cut))


dict_sum_energy = {}
dict_sum_energy['matrix_3_3'] = sum_energies_3
dict_sum_energy['matrix_5_5'] = sum_energies_5

with open('output/energy_sums.json', 'w') as fp:
    json.dump(dict_sum_energy, fp)
