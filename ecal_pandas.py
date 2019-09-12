import ROOT
from ROOT import RooRealVar,RooCBShape,RooDataHist,RooArgList,RooFit
from ROOT import gROOT,gStyle,gPad
from ROOT import std
import csv
import json
from array import array
import os

import numpy as np


import sys
sys.path.insert(0, 'utils/')
import pandas as pd

outstr = 'timeuniformity_C3_MCP2'
#trees_path = '/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Oct2018/ntuples_v5/'
trees_path = '/eos/user/n/nchernya/ntuples/ECAL/Upgrade/ECAL_TB_Oct2018/ntuples_v7/'

dir_str = '/11_09_2019/'
plot_folder = 'plots/'+dir_str
output_folder = 'output/'+dir_str
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

dict_C3_energy_scan = {}
##Now add runs for C2,C3,C4 
reader = csv.reader(open('data/energyscans19C_and3x3.csv', 'r'))
#reader = csv.reader(open('data/energyscans19C_3x3.csv', 'r'))
for row in reader:
    run,energy,crystal,_,_,_ = row
    if ('C3' in crystal) or ('C2' in crystal) or ('C4' in crystal):
        if energy in dict_C3_energy_scan.keys() :
            dict_C3_energy_scan[energy].append(run)
        else : 
            dict_C3_energy_scan[energy] = []
            dict_C3_energy_scan[energy].append(run)        

            
## Read crystal centers prepared by Simone 
reader = csv.reader(open('data/crystalscenters.csv', 'r'))
dict_crystal_centers={}
for row in reader:
    crystal,xpos,ypos = row
    if crystal!='Crystal' :
        dict_crystal_centers[crystal] = [int(xpos),int(ypos)]
        
dict_crystal_centers['C3_3x3'] =[4,5]       
        
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
matrix_5 = 'B5,C5,D5,B4,C4,D4,B3,C3,D3,B2,C2,D2,A1,B1,C1,D1,E1'.split(',')        
matrix_3_C4 = 'B5,C5,D5,B4,C4,D4,B3,C3,D3'.split(',')        
matrix_3_C2 = 'B3,C3,D3,B2,C2,D2,B1,C1,D1'.split(',')        


energies = sorted([float(item) for item in dict_C3_energy_scan.keys()])
energies = [str(item) for item in energies]
round_energies = [round(float(energy),-1) for energy in energies]
if round_energies[-1] ==240 : round_energies[-1] = 250.


dict_df_energy = {}

for energy in energies:
#for energy in [energies[-1]]:
    
    runs = dict_C3_energy_scan[energy]
    tree = ROOT.TChain("h4")
    for run in runs:
        tree.Add("%s/ECAL_H4_October2018_%s.root"%(trees_path,run))
    pos_cut = 4
    data = []
    for evt in tree:
        ############CHange e3x3 cut per energy!!!!!########
       # if evt.n_tracks==1 and evt.e3x3<10000. and evt.e3x3>8300. and evt.X[0] > -9. and evt.X[0] < 1. and evt.Y[0]>-1. and evt.Y[0]<8. : 
        if evt.n_tracks==1 and evt.X[0] > -9. and evt.X[0] < 1. and evt.Y[0]>-1. and evt.Y[0]<8. : 

            evt_dict = {}
            e = round(float(evt.Energy),-1)
            if e ==240 : e = 250.
            evt_dict['Energy'] = e
            evt_dict['run'] = run
            evt_dict['e3x3'] = evt.e3x3
            evt_dict['seed'] = evt.seed
            evt_dict['fit_ampl_MCP1'] = evt.fit_ampl[evt.MCP1]/evt.b_rms[evt.MCP1]
            evt_dict['fit_ampl_MCP2'] = evt.fit_ampl[evt.MCP2]/evt.b_rms[evt.MCP2]
            evt_dict['fit_time_MCP1'] = evt.fit_time[evt.MCP1]
            evt_dict['fit_time_MCP2'] = evt.fit_time[evt.MCP2]
            evt_dict['fit_time_VFE_CLK'] = evt.fit_time[evt.VFE_CLK]            
            for xstal in working:
                evt_dict['fit_ampl_'+xstal] = evt.fit_ampl[getattr(evt, xstal)]
                evt_dict['amp_max_'+xstal] = evt.amp_max[getattr(evt, xstal)]
                evt_dict['noise_'+xstal] = evt.b_rms[getattr(evt, xstal)]
                evt_dict['fit_time_'+xstal] = evt.fit_time[getattr(evt, xstal)]
                evt_dict['fit_terr_'+xstal] = evt.fit_terr[getattr(evt, xstal)]
                evt_dict['%s'%xstal] = getattr(evt,xstal)
            data.append(evt_dict)

    df = pd.DataFrame(data) 
    dict_df_energy[energy] = df

    pd_name = '%s/ECAL_H4_October2018_pandas_C3beamscan_full5_%s.csv'%(trees_path,int(e))
   # pd_name = '%s/ECAL_H4_October2018_pandas_C3_3x3_%s.csv'%(trees_path,int(e))
    df.to_csv(pd_name)
