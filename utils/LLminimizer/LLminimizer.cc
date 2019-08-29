#include "Riostream.h" 
#include "LLminimizer.h"
 
//  constructor
LLminimizer::LLminimizer(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double* dt_ecal_mcp2,double *dt_mcp2_mcp1, double *energy, TString crystal, int print_level):
    minuit_print_level_(print_level)
{
    SetData(events_num,mcp1_ampl,mcp2_ampl,dt_ecal_mcp1,dt_ecal_mcp2,dt_mcp2_mcp1,energy);
    SetCrystal(crystal);
}
 
// member functions
void LLminimizer::SetData(int events_num, double *mcp1_ampl,double *mcp2_ampl,double *dt_ecal_mcp1,double *dt_ecal_mcp2,double *dt_mcp2_mcp1,double *energy)
{
    m_mcp1_ampl_.assign(mcp1_ampl,mcp1_ampl+events_num);    
    m_mcp2_ampl_.assign(mcp2_ampl,mcp2_ampl+events_num);    
    m_dt_ecal_mcp1_.assign(dt_ecal_mcp1,dt_ecal_mcp1+events_num);    
    m_dt_ecal_mcp2_.assign(dt_ecal_mcp2,dt_ecal_mcp2+events_num);    
    m_dt_mcp2_mcp1_.assign(dt_mcp2_mcp1,dt_mcp2_mcp1+events_num);    
    m_data_size_ = events_num;
    m_energy_.assign(energy,energy+events_num);    


}
void LLminimizer::SetCrystal(TString crystal)
{
    m_crystal_ = crystal;
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
    double n_cb1 = par[10];
    double n_cb2 = par[11];  
    //double alpha_cb3 = par[10];
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

            m_Cm_.push_back(minimizer->X()[0]);
            m_Cm_e_.push_back(minimizer->Errors()[0]);

            m_alpha_.push_back(minimizer->X()[5]);
            m_alpha_e_.push_back(minimizer->Errors()[5]);
    
            m_beta_.push_back(minimizer->X()[6]);
            m_beta_e_.push_back(minimizer->Errors()[6]);
   
            m_a1_ = minimizer->X()[1];
            m_b1_ = minimizer->X()[2];
            m_a2_ = minimizer->X()[3];
            m_b2_ = minimizer->X()[4];
            m_gamma_ = minimizer->X()[7];
            
            m_a1_e_ = minimizer->Errors()[1];
            m_b1_e_ = minimizer->Errors()[2];
            m_a2_e_ = minimizer->Errors()[3];
            m_b2_e_ = minimizer->Errors()[4];
            m_gamma_e_ = minimizer->Errors()[7];
    
            delete minimizer;

    return 0;
}



double LLminimizer::NegLogLikelihoodSimultaneous(const double* par)
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
    double n_cb1 = par[10];
    double n_cb2 = par[11]; 
           
    for(int iSample=0; iSample<m_data_size_; ++iSample)
    {   
        int par_num = 12;
        if (m_energy_[iSample]==100) par_num=12;
        if (m_energy_[iSample]==150) par_num=15;       
        if (m_energy_[iSample]==200) par_num=18;
        if (m_energy_[iSample]==250) par_num=21;
           
        if (m_energy_[iSample]!=50){
            Cm = par[par_num];
            alpha = par[par_num+1];
            beta = par[par_num+2];        
        }   
        /*
        par_num = 24;
        if (m_energy_[iSample]==100) par_num=24;
        if (m_energy_[iSample]==150) par_num=28;       
        if (m_energy_[iSample]==200) par_num=32;
        if (m_energy_[iSample]==250) par_num=36;
           
        if (m_energy_[iSample]!=50){
            alpha_cb1 = par[par_num];
            alpha_cb2 = par[par_num+1];
        //    n_cb1 = par[par_num+2];
         //  n_cb2 = par[par_num+3];        
        }          
        */
         par_num = 24;
        if (m_energy_[iSample]==100) par_num=24;
        if (m_energy_[iSample]==150) par_num=26;       
        if (m_energy_[iSample]==200) par_num=28;
        if (m_energy_[iSample]==250) par_num=30;
           
        if (m_energy_[iSample]!=50){
            alpha_cb1 = par[par_num];
            alpha_cb2 = par[par_num+1];
        //    n_cb1 = par[par_num+2];
         //  n_cb2 = par[par_num+3];        
        }          
        
      
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

        nll += log(prob_ecal_mcp1 * prob_ecal_mcp2 * prob_mcp2_mcp1);            
    }

    m_nll_.push_back(nll);

    return -nll;
}


int LLminimizer::MinimizeNLLSimultaneous()
{

            ROOT::Math::Functor chi2(this, &LLminimizer::NegLogLikelihoodSimultaneous, 32);  //24  //40
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
            minimizer->SetLimitedVariable(10, "n_cb1",10,0.1,1.,100 );
            minimizer->SetLimitedVariable(11, "n_cb2",10,0.1,1.,100 );    
            minimizer->SetLimitedVariable(12, "Cm100",25.,1.,5,50. );
            minimizer->SetLimitedVariable(13, "alpha100",2500.,1.,2300.,2700. );
            minimizer->SetLimitedVariable(14, "beta100",1226.,1.,900.,1300. );
            minimizer->SetLimitedVariable(15, "Cm150",25.,1.,5,50. );
            minimizer->SetLimitedVariable(16, "alpha150",2500.,1.,2300.,2700. );
            minimizer->SetLimitedVariable(17, "beta150",1226.,1.,900.,1300. );
            minimizer->SetLimitedVariable(18, "Cm200",25.,1.,5,50. );
            minimizer->SetLimitedVariable(19, "alpha200",2500.,1.,2300.,2700. );
            minimizer->SetLimitedVariable(20, "beta200",1226.,1.,900.,1300. );
            minimizer->SetLimitedVariable(21, "Cm250",25.,1.,5,50. );
            minimizer->SetLimitedVariable(22, "alpha250",2500.,1.,2300.,2700. );
            minimizer->SetLimitedVariable(23, "beta250",1226.,1.,900.,1300. );  
    
            minimizer->SetLimitedVariable(24, "alpha_cb1_100",3,1.,0,10 );
            minimizer->SetLimitedVariable(25, "alpha_cb2_100",3,1.,0,10 );
         //   minimizer->SetLimitedVariable(26, "n_cb1_100",10,0.1,1.,100 );
        //    minimizer->SetLimitedVariable(27, "n_cb2_100",10,0.1,1.,100 );
    
            minimizer->SetLimitedVariable(28, "alpha_cb1_150",3,1.,0,10 );
            minimizer->SetLimitedVariable(29, "alpha_cb2_150",3,1.,0,10 );
          //  minimizer->SetLimitedVariable(30, "n_cb1_150",10,0.1,1.,100 );
          //  minimizer->SetLimitedVariable(31, "n_cb2_150",10,0.1,1.,100 );
    
            minimizer->SetLimitedVariable(32, "alpha_cb1_200",3,1.,0,10 );
            minimizer->SetLimitedVariable(33, "alpha_cb2_200",3,1.,0,10 );
       //     minimizer->SetLimitedVariable(34, "n_cb1_200",10,0.1,1.,100 );
         //   minimizer->SetLimitedVariable(35, "n_cb2_200",10,0.1,1.,100 );
    
            minimizer->SetLimitedVariable(36, "alpha_cb1_250",3,1.,0,10 );
            minimizer->SetLimitedVariable(37, "alpha_cb2_250",3,1.,0,10 );
         //   minimizer->SetLimitedVariable(38, "n_cb1_250",10,0.1,1.,100 );
         //   minimizer->SetLimitedVariable(39, "n_cb2_250",10,0.1,1.,100 );
    
            //---fit
            minimizer->Minimize();

            m_Cm_.push_back(minimizer->X()[0]);
            m_Cm_.push_back(minimizer->X()[12]);
            m_Cm_.push_back(minimizer->X()[15]);
            m_Cm_.push_back(minimizer->X()[18]);
            m_Cm_.push_back(minimizer->X()[21]);
    
            m_Cm_e_.push_back(minimizer->Errors()[0]);
            m_Cm_e_.push_back(minimizer->Errors()[12]);
            m_Cm_e_.push_back(minimizer->Errors()[15]);
            m_Cm_e_.push_back(minimizer->Errors()[18]);
            m_Cm_e_.push_back(minimizer->Errors()[21]);

            m_alpha_.push_back(minimizer->X()[5]);
            m_alpha_.push_back(minimizer->X()[13]);
            m_alpha_.push_back(minimizer->X()[16]);
            m_alpha_.push_back(minimizer->X()[19]);
            m_alpha_.push_back(minimizer->X()[22]);
    
            m_alpha_e_.push_back(minimizer->Errors()[5]);
            m_alpha_e_.push_back(minimizer->Errors()[13]);
            m_alpha_e_.push_back(minimizer->Errors()[16]);
            m_alpha_e_.push_back(minimizer->Errors()[19]);
            m_alpha_e_.push_back(minimizer->Errors()[22]);
    
            m_beta_.push_back(minimizer->X()[6]);
            m_beta_.push_back(minimizer->X()[14]);
            m_beta_.push_back(minimizer->X()[17]);
            m_beta_.push_back(minimizer->X()[20]);
            m_beta_.push_back(minimizer->X()[23]);    
    
            m_beta_e_.push_back(minimizer->Errors()[6]);
            m_beta_e_.push_back(minimizer->Errors()[14]);
            m_beta_e_.push_back(minimizer->Errors()[17]);
            m_beta_e_.push_back(minimizer->Errors()[20]);
            m_beta_e_.push_back(minimizer->Errors()[23]);   
    
            m_alpha_cb1_.push_back(minimizer->X()[8]);
            m_alpha_cb1_.push_back(minimizer->X()[24]);
            m_alpha_cb1_.push_back(minimizer->X()[26]);
            m_alpha_cb1_.push_back(minimizer->X()[28]);
            m_alpha_cb1_.push_back(minimizer->X()[30]);

            m_alpha_cb2_.push_back(minimizer->X()[9]);
            m_alpha_cb2_.push_back(minimizer->X()[25]);
            m_alpha_cb2_.push_back(minimizer->X()[27]);
            m_alpha_cb2_.push_back(minimizer->X()[29]);
            m_alpha_cb2_.push_back(minimizer->X()[31]);
    
            m_n_cb1_.push_back(minimizer->X()[10]);
            m_n_cb2_.push_back(minimizer->X()[11]);


    /*
            m_alpha_cb1_.push_back(minimizer->X()[8]);
            m_alpha_cb1_.push_back(minimizer->X()[24]);
            m_alpha_cb1_.push_back(minimizer->X()[28]);
            m_alpha_cb1_.push_back(minimizer->X()[32]);
            m_alpha_cb1_.push_back(minimizer->X()[36]);

            m_alpha_cb2_.push_back(minimizer->X()[9]);
            m_alpha_cb2_.push_back(minimizer->X()[25]);
            m_alpha_cb2_.push_back(minimizer->X()[29]);
            m_alpha_cb2_.push_back(minimizer->X()[33]);
            m_alpha_cb2_.push_back(minimizer->X()[37]);

            m_n_cb1_.push_back(minimizer->X()[10]);
            m_n_cb1_.push_back(minimizer->X()[26]);
            m_n_cb1_.push_back(minimizer->X()[30]);
            m_n_cb1_.push_back(minimizer->X()[34]);
            m_n_cb1_.push_back(minimizer->X()[38]);
    
            m_n_cb2_.push_back(minimizer->X()[11]);
            m_n_cb2_.push_back(minimizer->X()[27]);
            m_n_cb2_.push_back(minimizer->X()[31]);
            m_n_cb2_.push_back(minimizer->X()[35]);
            m_n_cb2_.push_back(minimizer->X()[39]);
            */
    
            m_alpha_cb1_e_.push_back(minimizer->Errors()[8]);
            m_alpha_cb2_e_.push_back(minimizer->Errors()[9]);     
    
            m_n_cb1_e_.push_back(minimizer->Errors()[10]);
            m_n_cb2_e_.push_back(minimizer->Errors()[11]);  
       
            m_a1_ = minimizer->X()[1];
            m_b1_ = minimizer->X()[2];
            m_a2_ = minimizer->X()[3];
            m_b2_ = minimizer->X()[4];
            m_gamma_ = minimizer->X()[7];
            
            m_a1_e_ = minimizer->Errors()[1];
            m_b1_e_ = minimizer->Errors()[2];
            m_a2_e_ = minimizer->Errors()[3];
            m_b2_e_ = minimizer->Errors()[4];
            m_gamma_e_ = minimizer->Errors()[7];
    
            delete minimizer;

    return 0;
}


