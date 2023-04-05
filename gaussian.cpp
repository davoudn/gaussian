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

struct single{};
struct multi{};

template <typename V>
struct drawGaussian;

template <typename T>
struct unif;

struct unif<double> {
	unif (double a1, double a2) {
		std::srand(std::time(nullptr));
	}
	double operator (){
		return ((double)std::rand()/RAND_MAX)*(a2-a1);
	}
};

struct unif<int> {
        unif (int a1, int a2) {
                std::srand(std::time(nullptr));
        }
        double operator () {
		int dummy = ((double)std::rand()/RAND_MAX)*(a2-a1);
                return dummy;
        }
};

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
       drawGaussian<single> m_grand;
       int m_dim;
       arma::Col<double> m_x;
       drawGaussian(int _dim):m_x(_dim), m_dim(_dim),m_grand(){}

       const arma::Col<double>& box_muller() {
	       std::for_each(m_x.begin(), m_x.end(), [this](double& x) {x=this->m_grand.box_muller()});
	       return m_x;
       }
};


template <typename T>
struct pdf {
        std::vector <arma::Col<double>> m_configs;
        pdf(){}
        virtual T energy();
        virtual arma::Col<T> force ();
        virtual T delta_energy(int n,T _xnew);
        virtual int get_size ();
	virtual double ratio(int n, T _xnew);
	void record (){}
};


template <>
struct gaussian::public pdf<double> {
	arma::Mat<double> m_Q, m_U;
	arma::Col<double> m_lambda, m_x;
        int dim;
        pdf (arma::Mat<double> _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
						 ,lambda(_dim),x(_dim)
		                                 ,m_dim(_dim){}    

	double energy()	{
	  return (m_x.t() * ( m_Q * m_x ))/2.0;
	}

	arma::Col<double> force (){
          return  m_Q * m_x;
	}

        double delta_energy(int n, double _xnew){
	  double delta { _xnew - m_x[n]}, dummy {0.0};	  
	  for (int i=0;i<dim; i++)
              dummy+= m_x[i]*delta;
	  dummy+= 0.5 * std::pow(delta,2.0);
          return dummy;
        }


        double ratio(int n, double _xnew){
          double delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<dim; i++)
              dummy+= m_x[i]*delta;
          dummy+= 0.5 * std::pow(delta,2.0);
          return std::exp(-dummy);
        }

	int get_size (){return  dim;}
};

template <typename PDF>
struct move();

/* some famous generic moves !! */

template <typename PDF>
struct langevin::public move< PDF >
{
  double m_D, m_dt;
  drawGaussian<multi> m_grand();

  langevin (PDF& _pdf, double _D, double _dt):drawSample< PDF >(_pdf), m_D(_D), m_dt(_dt), m_grand(_pdf.get_size())
  {}

  void ()
  {
    m_pdf.m_x = m_pdf.m_x + m_dt * m_pdf.force() + std::pow (2.0 * m_D * dt, 0.5) * m_grand();
    acc = true;
    return;
  }

};

template <typename PDF>
struct move_local::public move<PDF>
{
        unif <int>   m_u1;
        unif<double> m_u2;
  metropolice(PDF _pdf):drawSample<PDF>(_pdf),m_u1(0,_pdf.get_size()-1),m_u2(0.0,1.0){}

  void move()
  {
    acc = false;	  
    int site = m_u1();
    double x_new = m_pdf.m_x[site] + (m_u2()-0.5);
    if (  m_pdf.ratio(site,x_new) > m_u2() ){
       acc = true;
       m_pdf.m_x[site] = x_new;
    }
    return;
  }
};


template <typename PDF>
struct meassure;

template <typename PDF>
struct drawSample
{
	PDF& m_pdf;
	std::vector <arma::Col<double>> m_configs;
	std::vector<move<PDF>> m_moves;
	std::vector<meassures<PDF>> m_meassures;
        drawSample (PDF& _pdf):m_pdf(_pdf), m_configs(1,arma::Col<double>(_pdf.get_size())) {}

        void meassure(){
               for (std::vector<meassure<PDF>>::iterator it = m_meassures.begin(), it!=m_meassures.end(), it++){
                       it->meassure();
                }
        return;	
	}

	void move(){
		for (std::vector<move<PDF>>::iterator it = m_moves.begin(), it!=m_moves.end(), it++){
		       it->move();
		       if (it->acc == true ){
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

	void add_move (move<PDF>& _move){
              m_moves.push_back(_move);
	}

	void add_meassure ( meassure<PDF>& _meassure){
	      m_meassures.push_back (_meassure);
	}

	const arma::Col<double>& get_current_config () {
		return m_pdf.m_x;
	}
};


/* smapling gaussian with only local move */
template <>
struct meassure_correlation:: public meassure <gaussian>{
	gassuain m_gaussian;
	meassure_correlation () {
		
	}
};

template <>
struct drawGaussian::drawSample<gaussian> {
       local_move<gaussian> m_local_move;
       std::shared_ptr<gaussain> m_gaussian;
       drawGaussian (gaussian _gaussian):drawSample<gaussian>(_gaussian) {
                    add_move(m_local_move);
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
	std::vector<double> v(1,0.0);
        gaussian<single> g = gaussian<single>();

	int N{100}, N_draw{1000000};
	double field {10.0};

	for (int i=0; i < N_draw; i++)
            v.push_back( g.box_muller());

        histogram<single> hist = histogram<single> (v);
        hist.make(N,field);
}
