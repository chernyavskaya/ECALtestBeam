#include "Riostream.h" 
#include "LLminimizer.h"
 
//  constructor
//LLminimizer::LLminimizer(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1, TString crystal, TString energy)
LLminimizer::LLminimizer(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double *dt_ecal_mcp2,double *dt_mcp2_mcp1, TString crystal, TString energy)
{
    SetData(events_num,mcp1_ampl,mcp2_ampl,dt_ecal_mcp1,dt_ecal_mcp2,dt_mcp2_mcp1);
    SetCrystal(crystal);
    SetEnergy(energy);
}
 
// member functions
//void LLminimizer::SetData(vector<double> mcp1_ampl,vector<double> mcp2_ampl,vector<double> dt_ecal_mcp1,vector<double> dt_ecal_mcp2,vector<double> dt_mcp2_mcp1)
void LLminimizer::SetData(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double *dt_ecal_mcp2,double *dt_mcp2_mcp1)
{
    m_mcp1_ampl_.assign(mcp1_ampl,mcp1_ampl+events_num);    
    m_mcp2_ampl_.assign(mcp2_ampl,mcp2_ampl+events_num);    
    m_dt_ecal_mcp1_.assign(dt_ecal_mcp1,dt_ecal_mcp1+events_num);    
    m_dt_ecal_mcp2_.assign(dt_ecal_mcp2,dt_ecal_mcp2+events_num);    
    m_dt_mcp2_mcp1_.assign(dt_mcp2_mcp1,dt_mcp2_mcp1+events_num);    
    m_data_size_ = events_num;

  /*  m_mcp1_ampl_ = mcp1_ampl;
    m_mcp2_ampl_ = mcp2_ampl;
    m_dt_ecal_mcp1_ = dt_ecal_mcp1;
    m_dt_ecal_mcp2_ = dt_ecal_mcp2;
    m_dt_mcp2_mcp1_ = dt_mcp2_mcp1;
    m_data_size_ = mcp1_ampl.size();
*/
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
    double delta_1 = 0;
    double delta_2 = 0;
    double delta = 0;
    double pi = TMath::Pi();
	
    for(int iSample=0; iSample<m_data_size_; ++iSample)
    {	

		 double Cm = par[0];
		 double a1 = par[1];
		 double b1 = par[2];
		 double a2 = par[3];
		 double b2 = par[4];
		 double alpha = par[5];
		 double beta = par[6];
		 double gamma = par[7];

		 double sigma2_ecal_mcp1 = pow(Cm,2) + pow((a1/m_mcp1_ampl_[iSample]+b1) ,2) ;
		 double sigma2_ecal_mcp2 = pow(Cm,2) + pow((a2/m_mcp2_ampl_[iSample]+b2) ,2) ;
		 double sigma2_mcp2_mcp1 =  pow((a2/m_mcp2_ampl_[iSample]+b2) ,2) + pow((a1/m_mcp1_ampl_[iSample]+b1) ,2) ;	

       delta_1 = log(2*pi*sigma2_ecal_mcp1) + log(2*pi*sigma2_ecal_mcp2) + log(2*pi*sigma2_mcp2_mcp1); 
		 delta_2 =  pow((m_dt_ecal_mcp1_[iSample]-alpha),2)/sigma2_ecal_mcp1 +	pow((m_dt_ecal_mcp2_[iSample]-beta),2)/sigma2_ecal_mcp2 + pow((m_dt_mcp2_mcp1_[iSample]-gamma),2)/sigma2_mcp2_mcp1 ;
       nll = nll + delta_1 + delta_2;
    //		std::cout<<delta_1<<"   "<<delta_2<<std::endl; //debug
    		
//		std::cout<<"MCP amps : "<<m_mcp1_ampl_[0]<<"   "<<m_mcp2_ampl_[0]<<std::endl; //debug
    }

//    std::cout<<nll<<std::endl; //debug
    return nll;
}


int LLminimizer::MinimizeNLL()
{

            ROOT::Math::Functor chi2(this, &LLminimizer::NegLogLikelihood, 8);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(5000);
            minimizer->SetTolerance(1e-7);
            minimizer->SetPrintLevel(4);
            minimizer->SetFunction(chi2);
            minimizer->SetLimitedVariable(0, "Cm",25,1.,10,40 );
            minimizer->SetLimitedVariable(1, "a1",2,1e-1,0.01,5 );
            minimizer->SetLimitedVariable(2, "b1",30,1.,10.,50 );
            minimizer->SetLimitedVariable(3, "a2",2,1e-1,0.01,5 );
            minimizer->SetLimitedVariable(4, "b2",30,1.,10.,50 );
            minimizer->SetLimitedVariable(5, "alpha",2500,10,2000,3000 );
            minimizer->SetLimitedVariable(6, "beta",1150,10,900.,1400 );
            minimizer->SetLimitedVariable(7, "gamma",-4900.,10,-5500.,-4000 );
            //---fit
            minimizer->Minimize();
            m_Cm_ = minimizer->X()[0];
            m_a1_ = minimizer->X()[1];
            m_b1_ = minimizer->X()[2];
            m_a2_ = minimizer->X()[3];
            m_b2_ = minimizer->X()[4];
            m_alpha_ = minimizer->X()[5];
            m_beta_ = minimizer->X()[6];
            m_gamma_ = minimizer->X()[7];
            
            m_Cm_e_ = minimizer->Errors()[0];
            m_a1_e_ = minimizer->Errors()[1];
            m_b1_e_ = minimizer->Errors()[2];
            m_a2_e_ = minimizer->Errors()[3];
            m_b2_e_ = minimizer->Errors()[4];
            m_alpha_e_ = minimizer->Errors()[5];
            m_beta_e_ = minimizer->Errors()[6];
            m_gamma_e_ = minimizer->Errors()[7];

	return 0;
}


