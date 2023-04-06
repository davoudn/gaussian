#include <cmath>
#include <vector>
#include <algorithm>
#include <armadillo>


template <typename T>
struct pdf {
        std::vector <arma::Col<T>> m_configs;
//	arma::Col<T> m_x;
        pdf(){}
        virtual T energy(){ return T(1);}
        virtual arma::Col<T> force (){}
        virtual T delta_energy(int n,T _xnew){}
        virtual int size (){ return 0;}
        virtual double ratio(int n, T _xnew){}
        void record (){}
	void set_config () {}
};


struct gaussian:public pdf<double> {
        arma::Mat<double> m_Q, m_U, m_stupid_fix;
        arma::Col<double> m_lambda, m_x;
        int m_dim;
        gaussian (arma::Mat<double> _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
                                                 ,m_lambda(_dim),m_x(_dim),m_stupid_fix(1,1)
                                                 ,m_dim(_dim){}

        gaussian (arma::Mat<double>& _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
                                                 ,m_lambda(_dim),m_x(_dim),m_stupid_fix(1,1)
                                                 ,m_dim(_dim){}


        double energy() {
	  m_stupid_fix = m_x.t() * ( m_Q * m_x );
          return m_stupid_fix(0,0)*0.5;
//        return 0.0;
        }

        double energy(arma::Col<double>& _x) {
          m_stupid_fix = _x.t() * ( m_Q * _x );
          return m_stupid_fix(0,0);
        }
 
        double energy(arma::Col<double> _x) {
          m_stupid_fix = _x.t() * ( m_Q * _x );
          return m_stupid_fix(0,0)*0.5;
        }
 
        arma::Col<double> force (){
          return  m_Q * m_x;
        }

        double delta_energy(int n, double _xnew){
          double delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++) {
		  if ( i!=n ) dummy+= m_x[i]*delta;
          }
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

        void set_config (arma::Col<double> _x){
		m_x = _x; 
	}

        int size (){return  m_dim;}
};

