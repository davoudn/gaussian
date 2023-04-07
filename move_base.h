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



struct move_eigen_decompose_gaussian: public move_base<gaussian> {
       using t_pdf = std::shared_ptr<gaussian>;
       drawGaussian<single> m_grand;
       unif<int> m_urand;
       int m_nchanged;
       arma::vec m_eta;
       move_eigen_decompose_gaussian(t_pdf& _pdf, int _nchanged):move_base<gaussian>(_pdf), m_urand (0,_pdf->size()-1), m_nchanged(_nchanged), m_eta(_pdf->size())
       {
	 this->m_pdf->get_eig_decompose();
       }
       void move()
       {
         int l_site {0};
	 m_eta = this->m_pdf->m_U * this->m_pdf->m_x;
	 for (int i=0; i<m_nchanged; i++){
		 l_site = m_urand();
		 m_eta(l_site) = m_grand() * 2.0 * this->m_pdf->m_lambda(i);
	 }
	 this->m_pdf->m_x = this->m_pdf->m_U.t() * m_eta;
	 this->m_acc = true;
       }
};


struct move_chol_decompose_gaussian: public move_base<gaussian> {
       using t_pdf = std::shared_ptr<gaussian>;
       drawGaussian<single> m_grand;
       unif<int> m_urand;
       int m_nchanged;
       arma::vec m_eta;
       move_chol_decompose_gaussian(t_pdf& _pdf, int _nchanged):move_base<gaussian>(_pdf), m_urand (0,_pdf->size()-1), m_nchanged(_nchanged), m_eta(_pdf->size())
       {
         this->m_pdf->get_chol_decompose();
       }
       void move()
       {
         int l_site {0};
         m_eta = this->m_pdf->m_R.t() * this->m_pdf->m_x;
         for (int i=0; i<m_nchanged; i++){
                 l_site = m_urand();
                 m_eta(l_site) = m_grand() * 2.0 ;
         }
         this->m_pdf->m_x = this->m_pdf->m_R * m_eta;
         this->m_acc = true;
       }
};




