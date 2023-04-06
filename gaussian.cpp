#include <cstdlib>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <armadillo>
#include <map>
#include <memory>

struct single{};
struct multi{};

template <typename T>
struct unif;

template <>
struct unif<double> {
	double m_a1,m_a2;
	unif (double _a1, double _a2):m_a1(_a1), m_a2(_a2) {
		std::srand(std::time(nullptr));
	}
	double operator ()(){
		return ((double)std::rand()/RAND_MAX)*(m_a2-m_a1);
	}
};

template <>
struct unif<int> {
	double m_a1,m_a2;
        unif (int _a1, int _a2) {
                std::srand(std::time(nullptr));
        }
        double operator () (){
		int dummy = ((double)std::rand()/RAND_MAX)*(m_a2-m_a1);
                return dummy;
        }
};

template <typename V>
struct drawGaussian;


template <>
struct drawGaussian<single> {
	drawGaussian()
	{
            std::srand(std::time(nullptr)); // use current time as seed for random generator
	}
 
	double box_muller(){
		double u1{0.0}, u2{0.0}, eps{0.0001};
		while (u1 < eps || u2 < eps){
			u1 = (double)std::rand()/RAND_MAX;
			u2 = (double)std::rand()/RAND_MAX;
		}
		double x = std::pow(-2 * std::log(u1), 0.5) * std::cos(2 * 3.1415 * u2);
                return x;
        	}
};

template <>
struct drawGaussian<multi> {
       private:
         drawGaussian<single> m_grand;
         int m_dim;
         arma::Col<double> m_eta;
	public:
         drawGaussian(int _dim):m_eta(_dim), m_dim(_dim),m_grand(){}

     const arma::Col<double>& operator ()() {
	       std::for_each(m_eta.begin(), m_eta.end(), [this](double& x) {x=this->m_grand.box_muller();});
	       return m_eta;
       }
};


template <typename T>
struct pdf {
        std::vector <arma::Col<double>> m_configs;
        pdf(){}
        virtual T energy();
        virtual arma::Col<T> force ();
        virtual T delta_energy(int n,T _xnew);
        virtual int size ();
	virtual double ratio(int n, T _xnew);
	void record (){}
};


struct gaussian:public pdf<double> {
	arma::Mat<double> m_Q, m_U;
	arma::Col<double> m_lambda, m_x;
        int m_dim;
        gaussian (arma::Mat<double> _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
						 ,m_lambda(_dim),m_x(_dim)
		                                 ,m_dim(_dim){}    

	double energy()	{
// return (m_x.t() * ( m_Q * m_x ))*0.5;
	}

	arma::Col<double> force (){
          return  m_Q * m_x;
	}

        double delta_energy(int n, double _xnew){
	  double delta { _xnew - m_x[n]}, dummy {0.0};	  
	  for (int i=0;i<m_dim; i++)
              dummy+= m_x[i]*delta;
	  dummy+= 0.5 * std::pow(delta,2.0);
          return dummy;
        }


        double ratio(int n, double _xnew){
          double delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++)
              dummy+= m_x[i]*delta;
          dummy+= 0.5 * std::pow(delta,2.0);
          return std::exp(-dummy);
        }

	int size (){return  m_dim;}
};

template <typename PDF>
struct move_base{
 using t_pdf = typename std::shared_ptr<PDF>;
 t_pdf m_pdf;
 move_base(t_pdf& _pdf):m_pdf(_pdf){}
};
/* some famous generic moves !! */
template <typename PDF>
struct langevin:public move_base< PDF >
{
  using t_pdf = typename std::shared_ptr<PDF>;
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
	using t_pdf = typename std::shared_ptr<PDF>;
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


template <typename PDF>
struct meassure_base{
 using t_pdf = typename std::shared_ptr<PDF>;
 t_pdf m_pdf;
 meassure_base(t_pdf& _pdf ):m_pdf(_pdf){}


};

template <typename PDF>
struct meassure_2_point_correlation: public meassure_base <PDF>{
        using t_pdf = typename std::shared_ptr<PDF>;
        
	std::vector<arma::Mat<double>> m_results;
        arma::Mat<double> m_correlations;
	int m_size;
        meassure_2_point_correlation (t_pdf& _pdf):meassure_base<PDF>(_pdf),m_correlations(_pdf->size(),_pdf->size()), 
						   m_size(_pdf->size()) {}
        
        void meassure () {
             for (int i=0; i<m_size; i++)
		     for (int j=0; j<m_size; j++){
			     m_correlations(i,j) = this->m_pdf->m_x(i) * this->m_pdf->m_x(j);
		     }
	     m_results.push_back(m_correlations);
	     return;
	}
};


template <typename PDF>
struct drawSample
{
        using t_pdf = typename std::shared_ptr<PDF>;
	t_pdf m_pdf;
	std::vector <arma::Col<double>> m_configs;
	std::vector<move_base<PDF>> m_moves;
	std::vector<meassure_base<PDF>> m_meassures;
        drawSample (t_pdf _pdf):m_pdf(_pdf), m_configs(1,arma::Col<double>(_pdf->size())) {}

        void meassure(){
               for (typename std::vector<meassure_base<PDF>>::iterator it = m_meassures.begin(); it!=m_meassures.end(); it++){
                       it->meassure();
                }
        return;	
	}

	void move(){
		for (typename std::vector< move_base<PDF> >::iterator it = m_moves.begin(); it!=m_moves.end(); it++){
		       it->move();
		       if (it->m_acc == true ){
	                  meassure ();		      
			  record_configs(); 
		       }
		}
		return;
	}

	void record_configs () {
		m_configs.push_back(m_pdf.m_x);
		return;
	}

        void do_it(int n_iter){
             for (int it=0;it < n_iter; it++) {
                 move();
             }
	     return;
        }

	void add_move (move_base<PDF> _move){
              m_moves.push_back(_move);
	}

	void add_meassure ( meassure_base<PDF> _meassure){
	      m_meassures.push_back (_meassure);
	}

	const arma::Col<double>& get_current_config () {
		return m_pdf->m_x;
	}
};


struct simulateGaussian: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian (t_pdf _pdf):drawSample<gaussian>(_pdf) {
                    add_move(move_local<gaussian>(_pdf));
		    add_meassure (meassure_2_point_correlation<gaussian>(_pdf));
       }
};


/* ------------------------------------ */
template <typename V>
struct histogram {};

template <>
struct histogram<single> {
	std::vector<double> v;
        std::vector<double> h;
	int mesh_no;
	double field_length;

        histogram (std::vector<double>& _v):v(_v), h(1,0)
        {}

        void make(int _mesh_no, double _field_length) 
	{
		double _delta = _field_length/(_mesh_no+1);
		double _norm{0};
                h.resize(_mesh_no,0);
                std::for_each (v.begin(), v.end(), [this,_field_length,_mesh_no,_delta](double x) {
				for(int i=0;i<_mesh_no; ++i) 
                                   if ( x > ( i*_delta - _field_length/2.0) && x < ((i+1)*_delta - _field_length/2.0) ) {
				      this->h[i]+=1.0;
				      }
				   } );
                for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
			_norm+= *it;

//		std::for_each ( h.begin(), h.end(), [this,_norm](double& x) { x=x/_norm; });
  		for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
			*it = *it/_norm;
                std::fstream fout;

	        fout.open ("hist.dat",std::fstream::out);
		for (int i=0; i < _mesh_no; i++)
			fout << i<< " "  << i*_delta - _field_length/2.0 << " " << h[i] << "\n";
	}
};

int main ()
{

	return 0;
}


/*
{
	std::vector<double> v(1,0.0);
        gaussian<single> g = gaussian<single>();

	int N{100}, N_draw{1000000};
	double field {10.0};

	for (int i=0; i < N_draw; i++)
            v.push_back( g.box_muller());

        histogram<single> hist = histogram<single> (v);
        hist.make(N,field);
}
*/
