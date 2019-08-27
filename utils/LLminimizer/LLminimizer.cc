#include "Riostream.h" 
#include "LLminimizer.h"
 
//  constructor
LLminimizer::LLminimizer(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double* dt_ecal_mcp2,double *dt_mcp2_mcp1, TString crystal, TString energy, int print_level):
    minuit_print_level_(print_level)
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

    double Cm = par[0];
    double a1 = par[1];
    double b1 = par[2];
    double a2 = par[3];
    double b2 = par[4];
    double alpha = par[5];
    double beta = par[6];
    double gamma = par[7];
    double alpha_cb1 = par[8];
    double alpha_cb2 = par[9];
    //double alpha_cb3 = par[10];
    double n_cb1 = par[10];
    double n_cb2 = par[11];        
    //double n_cb3 = par[13];
                               
    for(int iSample=0; iSample<m_data_size_; ++iSample)
    {   
        double sigma2_ecal_mcp1 = pow(Cm, 2) + pow(a1/m_mcp1_ampl_[iSample], 2) + pow(b1, 2);
        double sigma2_ecal_mcp2 = pow(Cm, 2) + pow(a2/m_mcp2_ampl_[iSample], 2) + pow(b2, 2);
        double sigma2_mcp2_mcp1 = pow(a2/m_mcp2_ampl_[iSample], 2) + pow(b2, 2) + 
                                  pow(a1/m_mcp1_ampl_[iSample], 2) + pow(b1, 2) ;       
        
        auto prob_ecal_mcp1 = ROOT::Math::crystalball_pdf(m_dt_ecal_mcp1_[iSample], 
                                                          alpha_cb1, n_cb1,
                                                          sqrt(sigma2_ecal_mcp1),
                                                          alpha);
        auto prob_ecal_mcp2 = ROOT::Math::crystalball_pdf(m_dt_ecal_mcp2_[iSample], 
                                                          alpha_cb2, n_cb2,
                                                          sqrt(sigma2_ecal_mcp2),
                                                          beta);
        auto prob_mcp2_mcp1 = ROOT::Math::gaussian_pdf(m_dt_mcp2_mcp1_[iSample],
                                                       sqrt(sigma2_mcp2_mcp1),
                                                       gamma);

        //auto prob_mcp2_mcp1 = ROOT::Math::crystalball_pdf(m_dt_mcp2_mcp1_[iSample], 
        //                                                  alpha_cb3, n_cb3,
        //                                                  sqrt(sigma2_mcp2_mcp1),
        //                                                  gamma);        
        //auto prob_ecal_mcp1 = ROOT::Math::gaussian_pdf(m_dt_ecal_mcp1_[iSample],
        //                                               sqrt(sigma2_ecal_mcp1),
        //                                               alpha);
        //auto prob_ecal_mcp2 = ROOT::Math::gaussian_pdf(m_dt_ecal_mcp2_[iSample],
        //                                               sqrt(sigma2_ecal_mcp2),
        //                                               beta);
        
        nll += log(prob_ecal_mcp1 * prob_ecal_mcp2 * prob_mcp2_mcp1);            
    }

    m_nll_.push_back(nll);

    return -nll;
}


int LLminimizer::MinimizeNLL()
{

            ROOT::Math::Functor chi2(this, &LLminimizer::NegLogLikelihood, 12);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(1000);
            minimizer->SetTolerance(1.);
            minimizer->SetPrintLevel(minuit_print_level_);
            minimizer->SetFunction(chi2);
    
            minimizer->SetLimitedVariable(0, "Cm",25.,1.,5,50. );
            minimizer->SetLimitedVariable(1, "a1",800.,1,200.,2000 );
            minimizer->SetLimitedVariable(2, "b1",10.,1.,1.,20 );
            minimizer->SetLimitedVariable(3, "a2",800.,1,200.,2000 );
            minimizer->SetLimitedVariable(4, "b2",10.,1.,1.,20 );
            minimizer->SetLimitedVariable(5, "alpha",2500.,1.,2300.,2700. );
            minimizer->SetLimitedVariable(6, "beta",1226.,1.,900.,1300. );
            minimizer->SetLimitedVariable(7, "gamma",-4895.,1.,-5000.,-4000. ); 
            minimizer->SetLimitedVariable(8, "alpha_cb1",3,1.,0,10 );
            minimizer->SetLimitedVariable(9, "alpha_cb2",3,1.,0,10 );
            //minimizer->SetLimitedVariable(10, "alpha_cb3",10,1.,0,100);
            minimizer->SetLimitedVariable(10, "n_1",10,0.1,1.,100 );
            minimizer->SetLimitedVariable(11, "n_2",10,0.1,1.,100 );    
            //minimizer->SetLimitedVariable(13, "n_3",10,0.1,1.,20 );
    
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


    
            delete minimizer;

    return 0;
}


