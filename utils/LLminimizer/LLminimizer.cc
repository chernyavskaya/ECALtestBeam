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

    auto cb_ecal_mcp1 = new TF1("fitfunc_ecal_mcp1", "crystalball", -10000.,10000.);
    auto cb_ecal_mcp2 = new TF1("fitfunc_ecal_mcp2", "crystalball", -10000.,10000.);
    auto cb_mcp2_mcp1 = new TF1("fitfunc_mcp2_mcp1", "crystalball", -10000.,10000.);

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
        double alpha_cb1 = par[8];
        double alpha_cb2 = par[9];
        double alpha_cb3 = par[10];
    //    double n_cb1 = par[11];
   //     double n_cb2 = par[12];        
   //     double n_cb3 = par[13];
        double n_cb1 = 10.;
        double n_cb2 = 10.;        
        double n_cb3 = 10.;
      
      // when adding stochastic term
     //   double -> s1/TMath::Sqrt(m_mcp1_ampl_[iSample])) ,2) ;
        double sigma2_ecal_mcp1 = pow(Cm,2) + pow(a1/m_mcp1_ampl_[iSample], 2) + pow(b1,2);
        double sigma2_ecal_mcp2 = pow(Cm,2) + pow(a2/m_mcp2_ampl_[iSample],2) + pow(b2,2);
        double sigma2_mcp2_mcp1 =  pow(a2/m_mcp2_ampl_[iSample],2) + pow(b2,2) + pow(a1/m_mcp1_ampl_[iSample],2) + pow(b1,2) ;
   
 
        ///gaussian distibutions
   //     delta_1 = log(sigma2_ecal_mcp1) + log(sigma2_ecal_mcp2) + log(sigma2_mcp2_mcp1); 
   //     delta_2 =  pow((m_dt_ecal_mcp1_[iSample]-alpha),2)/sigma2_ecal_mcp1 + pow((m_dt_ecal_mcp2_[iSample]-beta),2)/sigma2_ecal_mcp2 + pow((m_dt_mcp2_mcp1_[iSample]-gamma),2)/sigma2_mcp2_mcp1 ;
   
    ///    nll = nll + delta_1 + delta_2;
        
        ///// crystal ball (TF1 parameters in order : norm, mu, sigma, alpha, N
        
        cb_ecal_mcp1->SetParameters(1 ,alpha,sqrt(sigma2_ecal_mcp1),alpha_cb1,n_cb1)  ; 
        cb_ecal_mcp2->SetParameters(1,beta,sqrt(sigma2_ecal_mcp2),alpha_cb2,n_cb2) ; 
        cb_mcp2_mcp1->SetParameters(1,gamma,sqrt(sigma2_mcp2_mcp1),alpha_cb3,n_cb3)  ;
            
        //double function_value = cb_ecal_mcp1->Eval(m_dt_ecal_mcp1_[iSample])*cb_ecal_mcp2->Eval(m_dt_ecal_mcp2_[iSample])*cb_mcp2_mcp1->Eval(m_dt_mcp2_mcp1_[iSample]);
            
        //nll = nll -2*log(function_value);
        if ( (cb_ecal_mcp1->Eval(m_dt_ecal_mcp1_[iSample])<1e-8) || (cb_ecal_mcp2->Eval(m_dt_ecal_mcp2_[iSample])<1e-8) || (cb_mcp2_mcp1->Eval(m_dt_mcp2_mcp1_[iSample])<1e-8)) {
            //std::cout << cb_ecal_mcp1->Eval(m_dt_ecal_mcp1_[iSample]) << " " << cb_ecal_mcp2->Eval(m_dt_ecal_mcp2_[iSample]) << " " << cb_mcp2_mcp1->Eval(m_dt_mcp2_mcp1_[iSample]) << std::endl;
            //std::cout << m_dt_ecal_mcp1_[iSample] << " " << m_dt_ecal_mcp2_[iSample] << " " << m_dt_mcp2_mcp1_[iSample] << std::endl;
            continue;
        }
        nll += (log(cb_ecal_mcp1->Eval(m_dt_ecal_mcp1_[iSample])) + log(cb_ecal_mcp2->Eval(m_dt_ecal_mcp2_[iSample])) + log(cb_mcp2_mcp1->Eval(m_dt_mcp2_mcp1_[iSample])));

        //std::cout << cb_ecal_mcp1->Eval(m_dt_ecal_mcp1_[iSample]) << " " << cb_ecal_mcp2->Eval(m_dt_ecal_mcp2_[iSample]) << " " << cb_mcp2_mcp1->Eval(m_dt_mcp2_mcp1_[iSample]) << std::endl;        
        
  ///std::cout << m_dt_ecal_mcp1_[iSample] << " " << m_dt_ecal_mcp2_[iSample] << " " << m_dt_mcp2_mcp1_[iSample] << std::endl; //debug
     
    }

  ///  std::cout<<nll<< std::endl ;
    m_nll_.push_back(nll);

    return nll;
}


int LLminimizer::MinimizeNLL()
{

            ROOT::Math::Functor chi2(this, &LLminimizer::NegLogLikelihood, 11);
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(100000);
            minimizer->SetMaxIterations(1000);
            minimizer->SetTolerance(1.);
            //minimizer->SetMaxFunctionCalls(1000000);
            //minimizer->SetMaxIterations(100000);
            //minimizer->SetTolerance(0.001);
            minimizer->SetPrintLevel(4);
            minimizer->SetFunction(chi2);
    
            minimizer->SetLimitedVariable(0, "Cm",25.,0.1,5,60. );
            minimizer->SetLimitedVariable(1, "a1",804.0,1,500.,2000 );
            minimizer->SetLimitedVariable(2, "b1",6.285,0.1,1.,20 );
            minimizer->SetLimitedVariable(3, "a2",900.364,1,500.,2000 );
            minimizer->SetLimitedVariable(4, "b2",8.237,0.1,1.,20 );
            minimizer->SetLimitedVariable(5, "alpha",2500,10,2000,3000 );
            minimizer->SetLimitedVariable(6, "beta",1226.,10,900.,1400 );
            minimizer->SetLimitedVariable(7, "gamma",-4895.,10,-5500.,-4000 ); 
            minimizer->SetLimitedVariable(8, "alpha_cb1",3,0.1,0.3,10 );
            minimizer->SetLimitedVariable(9, "alpha_cb2",3,0.1,.3,10 );
            minimizer->SetLimitedVariable(10, "alpha_cb3",3,0.1,.3,10 );
      //      minimizer->SetLimitedVariable(11, "n_1",10,0.1,1.,20 );
      //      minimizer->SetLimitedVariable(12, "n_2",10,0.1,1.,20 );    
      //      minimizer->SetLimitedVariable(13, "n_3",10,0.1,1.,20 );
    
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


