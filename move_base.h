#include <cmath>
#include <vector>
#include <armadillo>
#include <memory>


template <typename PDF>
struct move_base{
 using t_pdf = std::shared_ptr<PDF>;
 t_pdf m_pdf;
 move_base(t_pdf& _pdf):m_pdf(_pdf){}
};
/* some famous generic moves !! */
template <typename PDF>
struct langevin:public move_base< PDF >
{
  using t_pdf = std::shared_ptr<PDF>;
  double m_D, m_dt;
  bool m_acc;
  drawGaussian<multi> m_grand;
  public:
   langevin (t_pdf& _pdf, double _D, double _dt):move_base<PDF>(_pdf), m_D(_D), m_dt(_dt), m_grand(_pdf->size()),m_acc(false)
   {}

   void move()
   {
    this->m_pdf->m_x = this->m_pdf->m_x + m_dt * this->m_pdf->force() + std::pow (2.0 * m_D * m_dt, 0.5) * m_grand();
    m_acc = true;
    return;
   }

};

template <typename PDF>
struct move_local:public move_base<PDF>
{
private:
        using t_pdf = std::shared_ptr<PDF>;
        unif <int>   m_u1;
        unif<double> m_u2;
public:
 bool m_acc;
  move_local(t_pdf &_pdf):move_base<PDF>(_pdf), m_u1(0,_pdf->size()-1), m_u2(0.0,1.0),m_acc(false){}

  void move()
  {
    m_acc = false;
    int site = m_u1();
    double x_new = this->m_pdf->m_x[site] + (m_u2()-0.5);
    if (  this->m_pdf->ratio(site,x_new) > m_u2() ){
       m_acc = true;
       this->m_pdf->m_x[site] = x_new;
    }
    return;
  }
};

