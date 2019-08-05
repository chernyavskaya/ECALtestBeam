#include "Riostream.h" 
#include "LLminimizer.h"
 
//  constructor
LLminimizer::LLminimizer(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double *dt_ecal_mcp2,double *dt_mcp2_mcp1, TString crystal, TString energy)
{
    SetData(events_num,mcp1_ampl,mcp2_ampl,dt_ecal_mcp1,dt_ecal_mcp2,dt_mcp2_mcp1);
    SetCrystal(crystal);
    SetEnergy(energy);
}
 
// member functions
void LLminimizer::SetData(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double *dt_ecal_mcp2,double *dt_mcp2_mcp1)
{
    m_mcp1_ampl_.assign(mcp1_ampl,mcp1_ampl+events_num);    
    m_mcp2_ampl_.assign(mcp2_ampl,mcp2_ampl+events_num);    
    m_dt_ecal_mcp1_.assign(dt_ecal_mcp1,dt_ecal_mcp1+events_num);    
    m_dt_ecal_mcp2_.assign(dt_ecal_mcp2,dt_ecal_mcp2+events_num);    
    m_dt_mcp2_mcp1_.assign(dt_mcp2_mcp1,dt_mcp2_mcp1+events_num);    
    m_data_size_ = events_num;

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
        double s1 = par[8];
        double s2 = par[9];
      

        double sigma2_ecal_mcp1 = pow(Cm,2) + pow(( a1/m_mcp1_ampl_[iSample] + b1 + s1/TMath::Sqrt(m_mcp1_ampl_[iSample])) ,2) ;
        double sigma2_ecal_mcp2 = pow(Cm,2) + pow(( a2/m_mcp2_ampl_[iSample] + b2 + s2/TMath::Sqrt(m_mcp2_ampl_[iSample])) ,2) ;
        double sigma2_mcp2_mcp1 =  pow((a2/m_mcp2_ampl_[iSample]+b2+s2/TMath::Sqrt(m_mcp2_ampl_[iSample])) ,2) + pow((a1/m_mcp1_ampl_[iSample]+b1+s1/TMath::Sqrt(m_mcp1_ampl_[iSample])) ,2) ;
   
 
        
        delta_1 = log(sigma2_ecal_mcp1) + log(sigma2_ecal_mcp2) + log(sigma2_mcp2_mcp1); 
        delta_2 =  pow((m_dt_ecal_mcp1_[iSample]-alpha),2)/sigma2_ecal_mcp1 + pow((m_dt_ecal_mcp2_[iSample]-beta),2)/sigma2_ecal_mcp2 + pow((m_dt_mcp2_mcp1_[iSample]-gamma),2)/sigma2_mcp2_mcp1 ;
   
  
        nll = nll + delta_1 + delta_2;
    }

//    std::cout<<nll<<std::endl; //debug 
        m_nll_.push_back(nll);

    return nll;
}


int LLminimizer::MinimizeNLL()
{

            ROOT::Math::Functor chi2(this, &LLminimizer::NegLogLikelihood, 9);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(1000000);
            //minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(0.001);
            //minimizer->SetPrintLevel(4);
            minimizer->SetPrintLevel(0);
            minimizer->SetFunction(chi2);
    
            minimizer->SetLimitedVariable(0, "Cm",25.,0.1,5,60. );
            minimizer->SetLimitedVariable(1, "a1",1000.,1,1.,2000 );
            minimizer->SetLimitedVariable(2, "b1",15.,0.1,1.,50 );
            minimizer->SetLimitedVariable(3, "a2",1000.,1,1,2000 );
            minimizer->SetLimitedVariable(4, "b2",35.,0.1,1.,50 );
            minimizer->SetLimitedVariable(5, "alpha",2500,10,2000,3000 );
            minimizer->SetLimitedVariable(6, "beta",1150,10,900.,1400 );
            minimizer->SetLimitedVariable(7, "gamma",-4900.,10,-5500.,-4000 ); 
            minimizer->SetLimitedVariable(8, "s1",100.,0.1,1.,400 );
            minimizer->SetLimitedVariable(9, "s2",100.,0.1,1.,400 );

    
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
            m_s1_ = minimizer->X()[8];
            m_s2_ = minimizer->X()[9];
            
            m_Cm_e_ = minimizer->Errors()[0];
            m_a1_e_ = minimizer->Errors()[1];
            m_b1_e_ = minimizer->Errors()[2];
            m_a2_e_ = minimizer->Errors()[3];
            m_b2_e_ = minimizer->Errors()[4];
            m_alpha_e_ = minimizer->Errors()[5];
            m_beta_e_ = minimizer->Errors()[6];
            m_gamma_e_ = minimizer->Errors()[7];
            m_s1_e_ = minimizer->Errors()[8];
            m_s2_e_ = minimizer->Errors()[9];


    
            delete minimizer;

    return 0;
}


