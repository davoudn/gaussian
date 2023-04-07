#include <cmath>
#include <vector>
#include <armadillo>
#include <memory>


template <typename PDF>
struct move_base{
 using t_pdf = std::shared_ptr<PDF>;
 t_pdf m_pdf;
 bool m_acc;
 move_base(t_pdf& _pdf):m_pdf(_pdf),m_acc(false){}
 virtual void move (){std::cout << "Parent move. No move! \n  ";}
};
/* some famous generic moves !! */
template <typename PDF>
struct langevin:public move_base< PDF >{
  using t_pdf = std::shared_ptr<PDF>;
  double m_D, m_dt;
  drawGaussian<multi> m_grand;
  public:
   langevin (t_pdf& _pdf, double _D, double _dt):move_base<PDF>(_pdf), m_D(_D), m_dt(_dt), m_grand(_pdf->size())
   {}

   void move()
   {
    this->m_pdf->m_x = this->m_pdf->m_x + m_dt * this->m_pdf->force() + std::pow (2.0 * m_D * m_dt, 0.5) * m_grand();
    this->m_acc = true;
    return;
   }

};

template <typename PDF>
struct move_local:public move_base<PDF>{
        using t_pdf = std::shared_ptr<PDF>;
        unif <int>   m_u1;
        unif<double> m_u2;
  move_local(t_pdf &_pdf):move_base<PDF>(_pdf), m_u1(0,_pdf->size()), m_u2(0.0,1.0){}

  void move()
  {
    this->m_acc = false;
    int site = m_u1();
    double x_new = this->m_pdf->m_x[site] + (m_u2()-0.5)*1.0;
 //   std :: cout << site << " gg " << x_new << "\n";
    if (  this->m_pdf->ratio(site,x_new) > m_u2() ){
       this->m_acc = true;
       this->m_pdf->m_x[site] = x_new;
    }
    return;
  }
};



struct move_eig_decompose_gaussian: public move_base<gaussian> {
       using t_pdf = std::shared_ptr<gaussian>;
       drawGaussian<single> m_grand;
       unif<int> m_urand;
       int m_nchanged;
       arma::vec m_eta;
       move_eig_decompose_gaussian(t_pdf& _pdf, int _nchanged):move_base<gaussian>(_pdf), m_urand (0,_pdf->size()), m_nchanged(_nchanged), m_eta(_pdf->size())
       {
               this->m_pdf->set_eig_decompose();
       }
       void move()
       {
         int l_site {0};
	 m_eta = this->m_pdf->m_U.t() * this->m_pdf->m_x;
	 for (int i=0; i<m_nchanged; i++){
		 l_site = m_urand();
// std::cout << " ggg " << l_site << "\n";
		 m_eta(l_site) = m_grand() * std::pow ( 1.0 / this->m_pdf->m_lambda(l_site),0.5 );
	 }
	 this->m_pdf->m_x = this->m_pdf->m_U * m_eta;
	 this->m_acc = true;
       }
};


struct move_chol_decompose_gaussian: public move_base<gaussian> {
       using t_pdf = std::shared_ptr<gaussian>;
       drawGaussian<single> m_grand;
       unif<int> m_urand;
       int m_nchanged;
       arma::vec m_eta;
       move_chol_decompose_gaussian(t_pdf& _pdf, int _nchanged):move_base<gaussian>(_pdf), m_urand (0,_pdf->size()), m_nchanged(_nchanged), m_eta(_pdf->size())
       {
         this->m_pdf->set_chol_decompose();
       }
       void move()
       {
         int l_site {0};
         m_eta = this->m_pdf->m_R * this->m_pdf->m_x;
         for (int i=0; i<m_nchanged; i++){
                 l_site = m_urand();
                 m_eta(l_site) = m_grand();
         }
         this->m_pdf->m_x = this->m_pdf->m_R.t() * m_eta;
	 std::cout << this->m_pdf->m_x << "\n";
         this->m_acc = true;
       }
};


/* the algorithem presented in: Probability in the Engineering and Informational Sciences, 4, 1990, 369-389 
 * title : IMPROVING STOCHASTIC RELAXATION FOR GAUSSIAN RANDOM FIELDS 
 */

template <typename PDF>
struct move_stochastic_relaxation:public move_base< PDF >{
  private :
   using t_pdf = std::shared_ptr<PDF>;
   double m_omega, m_dt;
   int m_site;
   unif<int> m_urand;
   arma::mat m_force;
   drawGaussian<single> m_grand;
  public:
   move_stochastic_relaxation (t_pdf& _pdf, double _omega ,double _dt):move_base<PDF>(_pdf), m_dt(_dt), m_force (_pdf->size(),_pdf->size()),
	                                                               m_omega(_omega), m_urand (0,_pdf->size())
   {}

   void increment (){ m_site++; if (m_site == this->m_pdf->size()) m_site=0;  }
   void move()
   {
    double l_force {0.0};
    increment ();
    for (int i=0; i < this->m_pdf->size(); i++)
        if ( i != m_site )
	       l_force += this->m_pdf->m_Q(m_site,i) * this->m_pdf->m_x(i) * m_dt * m_omega;
    this->m_pdf->m_x(m_site) = -l_force + (1.0-m_omega) * m_dt * this->m_pdf->m_x(m_site ) + std::pow(m_dt,0.5) * m_grand();
    
    this->m_acc = true;
    return;
   }

};




