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
        virtual const arma::Col<T>& force (){ std::cout << "This is parrent pdf force (). No force!";}
        virtual T delta_energy(int n,T _xnew){}
        virtual int size (){ return 0;}
        virtual double ratio(int n, T _xnew){}
        void record (){}
	void set_config () {}
};


struct gaussian:public pdf<double> {
        arma::Mat<double> m_Q, m_U, m_stupid_fix,m_R, m_Q_inv;
        arma::Col<double> m_lambda, m_x, m_force;
        int m_dim;
        gaussian (arma::Mat<double> _Q, int _dim):m_Q(_Q),m_U(_dim,_dim), m_R(_dim,_dim),m_force(_dim), m_Q_inv(_dim,_dim)
                                                 ,m_lambda(_dim),m_x(arma::randu(_dim)),m_stupid_fix(1,1)
                                                 ,m_dim(_dim){}

      //  gaussian (arma::Mat<double>& _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
      //                                           ,m_lambda(_dim),m_x(_dim),m_stupid_fix(1,1)
                                                 //,m_dim(_dim){}


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
 
        const arma::Col<double>& force (){
		m_force = -m_Q * m_x;
	//	std::cout << m_force << "ff\n";
          return m_force;
        }

        double delta_energy(int n, double _xnew){
          double delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++) {
		  dummy+= m_Q(n,i)*m_x(i)*delta;
          }
          dummy+= 0.5 * std::pow(delta,2.0)*m_Q(n,n);
          return dummy;
        }


        double ratio(int n, double _xnew){
          double delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++){
               dummy+= m_Q(n,i)*m_x(i)*delta;
	  }
          dummy+= 0.5 * std::pow(delta,2.0)*m_Q(n,n);
          return std::exp(-dummy);
        }

        void set_config (arma::Col<double> _x){
		m_x = _x; 
	}

        int size (){return  m_dim;}

	const arma::Mat<double>& get_covariance ()
	{
		return m_Q;
	}
        void set_covariance_inv ( )
	{
		m_Q_inv = m_Q.i();
		return;
	}

        arma::Mat<double>& get_covariance_inv ()
        {
                return m_Q_inv;
        }

	
	void set_eig_decompose()
	{
		eig_sym(m_lambda, m_U, m_Q);
		std :: cout << m_lambda << "\n";
	}

	void set_chol_decompose()
        {
		m_R = arma::chol(m_Q);
                std :: cout << m_R << "\n";
        }


};

