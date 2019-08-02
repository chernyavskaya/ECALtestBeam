#include "Riostream.h" 
#include "LLminimizer.h"
 
//  constructor
LLminimizer::LLminimizer(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1, TString crystal, TString energy)
{
    SetData(mcp1_ampl,mcp2_ampl,dt_ecal_mcp1,dt_ecal_mcp2,dt_mcp2_mcp1);
    SetCrystal(crystal);
    SetEnergy(energy);
}
 
// member functions
void LLminimizer::SetData(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1)
{
    m_mcp1_ampl_ = mcp1_ampl;
    m_mcp2_ampl_ = mcp2_ampl;
    m_dt_ecal_mcp1_ = dt_ecal_mcp1;
    m_dt_ecal_mcp2_ = dt_ecal_mcp2;
    m_dt_mcp2_mcp1_ = dt_mcp2_mcp1;
    m_data_size_ = mcp1_ampl.size();
}
void LLminimizer::SetCrystal(TString crystal)
{
    m_crystal_ = crystal;
}
void LLminimizer::SetEnergy(TString energy)
{
    m_energy_ = energy;
}

double LLminimizer::NegLogLikelihood(const double* par)
{
    double nll = 0;
    double delta = 0;
    double pi = TMath::Pi();
	
    for(int iSample=0; iSample<=m_data_size_; ++iSample)
    {	

		 double Cm = par[0];
		 double a1 = par[1];
		 double b1 = par[2];
		 double a2 = par[3];
		 double b2 = par[4];
		 double alpha = par[5];
		 double beta = par[6];

		 double sigma2_ecal_mcp1 = pow(Cm,2) + pow((a1/m_mcp1_ampl_[iSample]+b1) ,2) ;
		 double sigma2_ecal_mcp2 = pow(Cm,2) + pow((a2/m_mcp2_ampl_[iSample]+b2) ,2) ;
		 double sigma2_mcp2_mcp1 =  pow((a2/m_mcp2_ampl_[iSample]+b2) ,2) + pow((a1/m_mcp1_ampl_[iSample]+b1) ,2) ;	

       delta = log(2*pi*sigma2_ecal_mcp1) + log(2*pi*sigma2_ecal_mcp2) + log(2*pi*sigma2_mcp2_mcp1) + pow((m_dt_ecal_mcp1_[iSample]-(alpha+beta)),2)/sigma2_ecal_mcp1 +
					pow((m_dt_ecal_mcp2_[iSample]-alpha),2)/sigma2_ecal_mcp2 + pow((m_dt_mcp2_mcp1_[iSample]-beta),2)/sigma2_mcp2_mcp1 ;
       nll += delta;
    }

    return nll;
}





