#ifndef LLMINIMIZER
#define LLMINIMIZER

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TError.h"
 
using namespace std;


class LLminimizer
{
private:
    vector<double>  m_mcp1_ampl_;
    vector<double>  m_mcp2_ampl_;
    vector<double>  m_dt_ecal_mcp1_;
    vector<double>  m_dt_ecal_mcp2_;
    vector<double>  m_dt_mcp2_mcp1_;
    long m_data_size_;
    TString m_crystal_;
    TString m_energy_;
 
public:
    LLminimizer(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1,  TString crystal, TString energy);
 
    void SetData(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1);
    void SetCrystal(TString crystal);
    void SetEnergy(TString energy);
    double NegLogLikelihood(const double* par=NULL);
 
};
 
#endif
